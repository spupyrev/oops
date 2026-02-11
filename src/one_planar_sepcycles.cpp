#include "graph_algorithms.h"
#include "logging.h"
#include "one_planar.h"

/// TODO: merge with similar impl
struct SepCycleTraversal {
  static constexpr size_t MAX_CROSSINGS_FOR_3EQ = 2048;
  static constexpr size_t MAX_CROSSINGS_FOR_3GT = 1280;

  SepCycleTraversal(SATModel& model, const InputGraph& graph, const Params& params, ForbiddenTuples& tuples)
    : model(model), graph(graph), params(params), n(graph.n), m((int)graph.edges.size()), 
      verbose(params.verbose), forbiddenTuples(tuples)
    {}

  bool init(int numPairs) {
    for (int e1 = 0; e1 < m; e1++) {
      for (int e2 = e1 + 1; e2 < m; e2++) {
        if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
          continue;
        possibleCrossings.push_back({e1, e2});
      }
    }

    LOG_IF(verbose, "sep-cycles constraints for %d pairs; |possible_crossings| = %d", 
           numPairs, possibleCrossings.size());

    isEdge = std::vector<std::vector<bool>> (n, std::vector<bool>(n, false));
    adjList = std::vector<std::vector<int>>(n, std::vector<int>());
    isCrossEdge = std::vector<std::vector<bool>>(n, std::vector<bool>(n, false));
    isKiteEdge = std::vector<std::vector<int>>(n, std::vector<int>(n, 0));

    // init free edges
    for (const auto& [u, v] : graph.edges) {
      setEdge(isEdge, u, v, true);
      adjList[u].push_back(v);
      adjList[v].push_back(u);
    }

    return true;
  }

  /// Find pairs of crossings that cannot happen in a 1-planar drawing
  size_t build2Clauses() {
    size_t numClauses2 = 0;
    size_t numEqFlowClauses2 = 0;
    size_t numEqFlowClauses3 = 0;

    // the main crossing
    for (size_t cross0 = 0; cross0 < possibleCrossings.size(); cross0++) {
      const auto [e1_cross0, e2_cross0] = possibleCrossings[cross0];
      const auto [u, v] = graph.edges[e1_cross0];
      const auto [x, y] = graph.edges[e2_cross0];
      CHECK(all_unique({x, y, u, v}));

      // circular order: u--x--v--y
      markCrossing(u, v, x, y);

      std::vector<std::pair<int, int>> takenCrossings0 = {{e1_cross0, e2_cross0}};
      numEqFlowClauses2 += findEqualFlow(x, y, u, v, takenCrossings0);
      numEqFlowClauses2 += findEqualFlow(u, v, x, y, takenCrossings0);

      // secondary crossings
      for (size_t cross1 = 0; cross1 < possibleCrossings.size(); cross1++) {
        const auto [e1_cross1, e2_cross1] = possibleCrossings[cross1];
        if (cross0 == cross1)
          continue;
        // circular order: u2--x2--v2--y2
        const auto [u2, v2] = graph.edges[e1_cross1];
        const auto [x2, y2] = graph.edges[e2_cross1];

        // check that the selected crossings are compatibe
        if (!canMarkCrossing(u2, v2, x2, y2))
          continue;
        CHECK(all_unique({e1_cross0, e2_cross0, e1_cross1, e2_cross1}));

        // stop early, if the crossing pair is already processed
        const CrossingPair crossPair(m, e1_cross0, e2_cross0, e1_cross1, e2_cross1);
        if (forbiddenTuples.contains(crossPair))
          continue;

        markCrossing(u2, v2, x2, y2);

        // need to search cycles both via (x, y) and via (u, v)
        std::vector<std::pair<int, int>> takenCrossings = {{e1_cross0, e2_cross0}, {e1_cross1, e2_cross1}};
        if (findSepCycle(x, y, u, v, takenCrossings) || findSepCycle(u, v, x, y, takenCrossings)) {
          // save the 2-clause
          numClauses2++;
          forbiddenTuples.insert(crossPair);
        }

        if (possibleCrossings.size() < MAX_CROSSINGS_FOR_3EQ) {
          numEqFlowClauses3 += findEqualFlow(x, y, u, v, takenCrossings);
          numEqFlowClauses3 += findEqualFlow(u, v, x, y, takenCrossings);
        }

        unmarkCrossing(u2, v2, x2, y2);
      }

      unmarkCrossing(u, v, x, y);
    }

    LOG_IF(verbose, "  added %'9d equal-flow 2-clauses", numEqFlowClauses2);
    LOG_IF(verbose, "  added %'9d equal-flow 3-clauses", numEqFlowClauses3);

    LOG_IF(verbose && clauseIndex, "num_extra_clauses: %d", clauseIndex);
    return numClauses2;
  }

  /// Find triples of crossings that cannot happen in a 1-planar drawing
  size_t build3Clauses() {
    if (possibleCrossings.size() > MAX_CROSSINGS_FOR_3GT) {
      LOG_IF(verbose, "  skipped building sep-cycles 3-clauses due to too many possible_crossings");
      return 0;
    }

    for (const auto& [u, v] : graph.edges) {
      CHECK(isEdge[u][v]);
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        CHECK(!isCrossEdge[i][j]);
        CHECK(!isKiteEdge[i][j]);
      }
    }

    size_t numClauses3 = 0;

    // the main crossing
    for (size_t idx_uv = 0; idx_uv < possibleCrossings.size(); idx_uv++) {
      const int e1 = possibleCrossings[idx_uv].first;
      const int e2 = possibleCrossings[idx_uv].second;
      const auto [u, v] = graph.edges[e1];
      const auto [x, y] = graph.edges[e2];
      CHECK(all_unique({x, y, u, v}));

      markCrossing(u, v, x, y);

      // secondary crossing-A
      for (size_t cr0 = 0; cr0 < possibleCrossings.size(); cr0++) {
        if (idx_uv == cr0) 
          continue;

        const int e1_cr0 = possibleCrossings[cr0].first;
        const int e2_cr0 = possibleCrossings[cr0].second;
        const auto [u1, v1] = graph.edges[e1_cr0];
        const auto [u2, v2] = graph.edges[e2_cr0];

        // skip the triple if this pair is already forbidden
        if (forbiddenTuples.contains(CrossingPair(m, e1, e2, e1_cr0, e2_cr0))) 
          continue;

        if (!canMarkCrossing(u1, v1, u2, v2))
          continue;

        markCrossing(u1, v1, u2, v2);

        // secondary crossing-B
        for (size_t cr1 = cr0 + 1; cr1 < possibleCrossings.size(); cr1++) {
          if (idx_uv == cr1) 
            continue;

          const int e1_cr1 = possibleCrossings[cr1].first;
          const int e2_cr1 = possibleCrossings[cr1].second;
          const auto [w1, z1] = graph.edges[e1_cr1];
          const auto [w2, z2] = graph.edges[e2_cr1];

          // Approximation: need to have at least two common vertices among the crossings
          if (unique_size({x, y, u, v, u1, v1, u2, v2, w1, z1, w2, z2}) >= 11)
            continue;

          // skip the triple if these pairs are already forbidden
          if (forbiddenTuples.contains(CrossingPair(m, e1, e2, e1_cr1, e2_cr1))) 
            continue;
          if (forbiddenTuples.contains(CrossingPair(m, e1_cr0, e2_cr0, e1_cr1, e2_cr1))) 
            continue;

          // compatibility checks with already marked stuff
          if (!canMarkCrossing(w1, z1, w2, z2))
            continue;

          const CrossingTriple triple(m, e1, e2, e1_cr0, e2_cr0, e1_cr1, e2_cr1);
          if (forbiddenTuples.contains(triple))
            continue;

          markCrossing(w1, z1, w2, z2);

          // need to search cycles both via (x, y) and via (u, v)
          CHECK(all_unique({e1, e2, e1_cr0, e2_cr0, e1_cr1, e2_cr1}));
          const std::vector<std::pair<int,int>> takenCrossings = {{e1,e2}, {e1_cr0,e2_cr0}, {e1_cr1,e2_cr1}};
          if (findSepCycle(x, y, u, v, takenCrossings) || findSepCycle(u, v, x, y, takenCrossings)) {
            // save the 3-clause
            numClauses3++;
            forbiddenTuples.insert(triple);
          }

          unmarkCrossing(w1, z1, w2, z2);
        }

        unmarkCrossing(u1, v1, u2, v2);
      }

      unmarkCrossing(u, v, x, y);
    }

    return numClauses3;
  }

private:
  void setEdge(std::vector<std::vector<bool>> &isEdge, int u, int v, bool value) {
    isEdge[u][v] = value;
    isEdge[v][u] = value;
  }

  void markKiteEdge(int u, int v) {
    isKiteEdge[u][v] += 1;
    isKiteEdge[v][u] += 1;

    if (isKiteEdge[u][v] == 1 && !isEdge[u][v]) {
      adjList[u].push_back(v);
      adjList[v].push_back(u);
    }
  }

  void unmarkKiteEdge(int u, int v) {
    isKiteEdge[u][v] -= 1;
    isKiteEdge[v][u] -= 1;

    if (isKiteEdge[u][v] == 0 && !isEdge[u][v]) {
      CHECK(adjList[u].back() == v);
      adjList[u].pop_back();
      CHECK(adjList[v].back() == u);
      adjList[v].pop_back();
    }
  }

  bool canMarkCrossing(int u, int v, int x, int y) const {
    // check that the selected crossings are compatibe
    if (isCrossEdge[u][v] || isCrossEdge[x][y])
      return false;
    if (isKiteEdge[u][v] || isKiteEdge[x][y])
      return false;
    if (isCrossEdge[u][x] || isCrossEdge[x][v] || isCrossEdge[v][y] || isCrossEdge[y][u])
      return false;

    return true;
  }

  void markCrossing(int u, int v, int x, int y) {
    CHECK(isEdge[u][v] && isEdge[x][y]);
    CHECK(!isCrossEdge[u][v] && !isCrossEdge[x][y]);
    CHECK(isKiteEdge[u][v] == 0 && isKiteEdge[x][y] == 0);
    // cross
    setEdge(isCrossEdge, u, v, true);
    setEdge(isCrossEdge, x, y, true);

    // kite: circular order is u--x--v--y
    CHECK(!isCrossEdge[u][x] && !isCrossEdge[x][v] && !isCrossEdge[v][y] && !isCrossEdge[y][u]);
    markKiteEdge(u, x);
    markKiteEdge(x, v);
    markKiteEdge(v, y);
    markKiteEdge(y, u);
  }

  void unmarkCrossing(int u, int v, int x, int y) {
    CHECK(isEdge[u][v] && isEdge[x][y]);
    CHECK(isCrossEdge[u][v] && isCrossEdge[x][y]);
    CHECK(isKiteEdge[u][v] == 0 && isKiteEdge[x][y] == 0);
    // uncross
    setEdge(isCrossEdge, u, v, false);
    setEdge(isCrossEdge, x, y, false);

    // kite: circular order is u--x--v--y
    unmarkKiteEdge(y, u);
    unmarkKiteEdge(v, y);
    unmarkKiteEdge(x, v);
    unmarkKiteEdge(u, x);
  }

  bool isFreeEdge(int u, int v) {
    return !isCrossEdge[u][v] && !isKiteEdge[u][v];
  };

  /// Find all separating cycles: paths from x to y avoiding vertices {x, y, u, v}
  /// A path may contain as many kite/cross edges as possible but only one free edge
  /// If a path contains one crossed edge, then disjoint paths cannot use the second one
  bool findSepCycle(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings) {
    const size_t numTakenCrossings = takenCrossings.size();
    CHECK(numTakenCrossings == 2 || numTakenCrossings == 3);
    bool foundSepCycle = false;

    auto processCycle = [&](const std::vector<int>& cycle) -> bool {
      CHECK(cycle.front() == x && cycle.back() == y && cycle.size() >= 3);
      size_t numFree = 0;
      size_t numCross = 0;
      for (size_t i = 0; i < cycle.size(); i++) {
        const int v1 = cycle[i];
        const int v2 = cycle[(i + 1) % cycle.size()];
        numFree += size_t(isFreeEdge(v1, v2));
        numCross += size_t(isCrossEdge[v1][v2]);
      }

      // Search for paths from crossing edge pairs (2 or 3 crossings)
      if ((numTakenCrossings == 2 || numTakenCrossings == 3) && numCross == numTakenCrossings) {
        if (numTakenCrossings == 2 || params.unsatLevel >= 2) {
          if (hasManyDisjointPaths(takenCrossings, cycle, numFree)) {
            foundSepCycle = true;
            return true; // stop processing
          }
        }
      }

      // Check if there are many edge-disjoint paths from u to v disjoint from the cycle
      if (numFree >= (size_t)std::min(graph.degree(u) - 1, graph.degree(v) - 1))
        return false; // do not stop processing
      const auto crossCycleEdges = getCrossCycleEdges(cycle, takenCrossings);
      const size_t numPaths = countEdgeDisjointPaths({u}, {v}, graph.adj, cycle, crossCycleEdges, numFree + 1);
      if (numPaths > numFree) {
        foundSepCycle = true;
        return true; // stop processing
      }
      return false;
    };

    const std::vector<int> avoidedVertices = {x, y, u, v};
    // First, check cycle (x -- c -- y)
    forEachCycle(adjList, 3, x, y, avoidedVertices, processCycle);
    if (foundSepCycle)
      return true;

    // Second, check cycle (x -- c1 -- c2 -- y)
    forEachCycle(adjList, 4, x, y, avoidedVertices, processCycle);
    if (foundSepCycle)
      return true;

    // Third, check larger cycles in certain cases
    if (numTakenCrossings == 2 || numTakenCrossings == 3 /* || true*/) {
      for (size_t cycleLen = 5; cycleLen <= 8; cycleLen++) {
        forEachCycle(adjList, cycleLen, x, y, avoidedVertices, processCycle);
        if (foundSepCycle)
          return true;
      }
    }

    return false;
  }

  /// Check if there are more than `numFree` edge-disjoint paths from u to v disjoint from the cycle.
  /// It is assumed that the cycle contains exactly one edge from each crossing.
  bool hasManyDisjointPaths(const std::vector<std::pair<int, int>> &takenCrossings,
                            const std::vector<int>& cycle, 
                            int numFree) const {
    const size_t numTakenCrossings = takenCrossings.size();
    CHECK(numTakenCrossings == 2 || numTakenCrossings == 3);
    // Extract each cross-cycle edge per crossing
    std::vector<EdgeTy> crossCycleEdges;
    crossCycleEdges.reserve(numTakenCrossings);
    for (size_t i = 0; i < numTakenCrossings; i++) {
      const auto [e1, e2] = takenCrossings[i];
      const EdgeTy edge1 = graph.edges[e1];
      const EdgeTy edge2 = graph.edges[e2];
      const bool edge1OnCycle = isEdgeOnCycle(cycle, edge1);
      const bool edge2OnCycle = isEdgeOnCycle(cycle, edge2);
      // Need exactly one edge from each crossing on the cycle
      if (edge1OnCycle == edge2OnCycle)
        return false;
      crossCycleEdges.push_back(edge1OnCycle ? edge2 : edge1);
    }
    CHECK(numTakenCrossings == crossCycleEdges.size());

    for (size_t i = 0; i < crossCycleEdges.size(); i++) {
      const auto [x, y] = crossCycleEdges[i];
      if (contains(cycle, x) || contains(cycle, y))
        return false;
    }

    int degreeBound = 0;
    for (const auto& [x, y] : crossCycleEdges) {
      degreeBound += std::min(graph.degree(x) - 1, graph.degree(y) - 1);
    }
    if (numFree >= degreeBound)
      return false;

    // Orientation masks are symmetric under flipping all pairs at once
    // (swap sources <-> targets globally), so fix pair 0 as (a -> b).
    // Only pairs 1..k-1 need flipping choices.
    const int numPairs = static_cast<int>(numTakenCrossings);
    const int numMasks = 1 << (numPairs - 1);
    bool hasValidOrientation = false;
    for (int mask = 0; mask < numMasks; mask++) {
      const auto [x0, y0] = crossCycleEdges[0];
      std::vector<int> sources = {x0};
      std::vector<int> targets = {y0};
      for (int i = 1; i < numPairs; i++) {
        const auto [x, y] = crossCycleEdges[i];
        const bool flip = ((mask >> (i - 1)) & 1) != 0;
        sources.push_back(flip ? y : x);
        targets.push_back(flip ? x : y);
      }
      // The mask is invalid
      if (overlap(sources, targets))
        continue;
      hasValidOrientation = true;

      const int numPaths = countEdgeDisjointPaths(sources, targets, graph.adj, cycle, crossCycleEdges, numFree + 1);
      if (numPaths <= numFree)
        return false;
    }

    if (!hasValidOrientation) {
      // This configuration cannot happen; hence, disable the clause
      return true;
    }

    CHECK(hasValidOrientation);
    return true;
  }

  int findEqualFlow(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings) {
    CHECK(takenCrossings.size() == 1 || takenCrossings.size() == 2);
    int numAddedClauses = 0;

    auto processCycle = [&](const std::vector<int>& cycle) -> bool {
      CHECK(cycle.front() == x && cycle.back() == y);
      int numFree = 0;
      bool hasCrossableEdge = false;
      for (size_t i = 0; i + 1 < cycle.size(); i++) {
        const int v1 = cycle[i];
        const int v2 = cycle[i + 1];
        numFree += int(isFreeEdge(v1, v2));
        if (isEdge[v1][v2] && isFreeEdge(v1, v2))
          hasCrossableEdge = true;
      }

      if (!hasCrossableEdge)
        return false; // do not stop processing
      if (numFree > std::min(graph.degree(u) - 1, graph.degree(v) - 1))
        return false; // do not stop processing

      const auto crossCycleEdges = getCrossCycleEdges(cycle, takenCrossings);
      // TOOD: might want to use "numFree"
      // const int flowLB = (takenCrossings.size() == 1 ? numFree : numFree + 1);
      const int numPaths = countEdgeDisjointPaths(u, v, graph.adj, cycle, crossCycleEdges, numFree);
      if (numPaths == numFree) {
        CHECK(numPaths > 0);
        // found a clause
        for (size_t i = 0; i + 1 < cycle.size(); i++) {
          const int v1 = cycle[i];
          const int v2 = cycle[i + 1];
          if (!isEdge[v1][v2] || !isFreeEdge(v1, v2))
            continue;
          // the edge must cross smth
          numAddedClauses++;
          if (takenCrossings.size() == 1) {
            addEqFlowClause(takenCrossings[0], v1, v2);
          } else {
            CHECK(takenCrossings.size() == 2);
            addEqFlowClause(takenCrossings[0], takenCrossings[1], v1, v2);
          }
        }
      } 
      return false; // do not stop processing
    };

    const std::vector<int> avoidedVertices = {x, y, u, v};
    // First, check cycle (x -- c -- y)
    forEachCycle(adjList, 3, x, y, avoidedVertices, processCycle);

    // Second, check cycle (x -- c1 -- c2 -- y)
    forEachCycle(adjList, 4, x, y, avoidedVertices, processCycle);

    // Thritd, check 5-cycles when takenCrossings.size() == 1
    // TODO: check always?
    if (takenCrossings.size() == 1/* || true*/) {
      forEachCycle(adjList, 5, x, y, avoidedVertices, processCycle);
    }

    return numAddedClauses;
  }

  /// Check if the edge is on the cycle
  bool isEdgeOnCycle(const std::vector<int>& cycle, const EdgeTy& edge) const {
    for (size_t i = 0; i < cycle.size(); i++) {
      const EdgeTy cycleEdge = make_edge(cycle[i], cycle[(i + 1) % cycle.size()]);
      if (edge == cycleEdge)
        return true;
    }
    return false;
  }

  /// Collect edges crossing the cycle
  std::vector<EdgeTy> getCrossCycleEdges(
      const std::vector<int>& cycle,
      const std::vector<std::pair<int, int>> &takenCrossings) const {
    std::vector<EdgeTy> crossCycleEdges;
    for (const auto& [e1, e2] : takenCrossings) {
      const bool e1_on_cycle = isEdgeOnCycle(cycle, graph.edges[e1]);
      const bool e2_on_cycle = isEdgeOnCycle(cycle, graph.edges[e2]);
      if (e1_on_cycle) {
        crossCycleEdges.push_back(graph.edges[e2]);
      } 
      if (e2_on_cycle) {
        crossCycleEdges.push_back(graph.edges[e1]);
      } 
    }
    return crossCycleEdges;
  }

  void addEqFlowClause(std::pair<int, int> cr1, int u, int v) {
    const auto [e1, e2] = cr1;
    const int divUV = graph.findDivIndex(u, v);
    model.addClause({
        model.getCross2Var(e1 + n, e2 + n, false), 
        model.getCross1Var(divUV, true)
    });
  }

  void addEqFlowClause(std::pair<int, int> cr1, std::pair<int, int> cr2, int u, int v) {
    const auto [e1, e2] = cr1;
    const auto [e3, e4] = cr2;
    const int divUV = graph.findDivIndex(u, v);
    model.addClause({
        model.getCross2Var(e1 + n, e2 + n, false), 
        model.getCross2Var(e3 + n, e4 + n, false), 
        model.getCross1Var(divUV, true)
    });
  }

private:
  SATModel& model;
  const InputGraph& graph;
  const Params& params;
  const int n;
  const int m;
  const int verbose;
  ForbiddenTuples& forbiddenTuples;

  std::vector<std::pair<int, int>> possibleCrossings;
  std::vector<std::vector<bool>> isEdge;
  std::vector<std::vector<int>> adjList;

  std::vector<std::vector<bool>> isCrossEdge;
  std::vector<std::vector<int>> isKiteEdge;

  // Debug
  int clauseIndex = 0;
};

/// Add constraints based on separating cycles
void encodeSepCyclesConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  const int verbose = params.verbose;
  const int numPairs = to_int(params.sepCycleConstraints);
  CHECK(2 <= numPairs && numPairs <= 3, "incorrect value of sepCycleConstraints");

  SepCycleTraversal sct(model, graph, params, model.getForbiddenTuples());
  if (sct.init(numPairs)) {
    // forbidden pairs (2-clauses)
    const size_t numClauses2 = sct.build2Clauses();
    LOG_IF(verbose, "  found %'9d 2-clauses", numClauses2);

    if (numPairs == 3 || params.unsatLevel >= 2) {
      // forbidden triples (3-clauses)
      const size_t numClauses3 = sct.build3Clauses();
      LOG_IF(verbose, "  found %'9d 3-clauses", numClauses3);
    }
  }
}
