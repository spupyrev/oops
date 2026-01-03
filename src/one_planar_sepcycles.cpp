#include "graph_algorithms.h"
#include "logging.h"
#include "one_planar.h"

#include <utility>

struct CrossingPair {
  int e1;
  int e2;
  int p1;
  int p2;
  int64_t key;

  explicit CrossingPair(int m, int e1_, int e2_, int p1_, int p2_)
      : e1(e1_), e2(e2_), p1(p1_), p2(p2_) {
    // Normalize ordering so that (e1,e2) < (p1,p2) lexicographically
    CHECK(0 <= e1 && e1 < e2 && e2 < m);
    CHECK(0 <= p1 && p1 < p2 && p2 < m);
    if (std::make_pair(e1, e2) >= std::make_pair(p1, p2)) {
      std::swap(e1, p1);
      std::swap(e2, p2);
    }
    CHECK(std::make_pair(e1, e2) < std::make_pair(p1, p2));

    // Construct the key
    key = (int64_t)e1 * (int64_t)m * (int64_t)m * (int64_t)m +
          (int64_t)e2 * (int64_t)m * (int64_t)m +
          (int64_t)p1 * (int64_t)m +
          (int64_t)p2;
  }

  void encode(SATModel& model, const InputGraph& graph) const {
    const int divE1 = e1 + graph.n;
    const int divE2 = e2 + graph.n;
    CHECK(canBeMerged(divE1, divE2, graph));
    const int divP1 = p1 + graph.n;
    const int divP2 = p2 + graph.n;
    CHECK(canBeMerged(divP1, divP2, graph));

    model.addClause({
        model.getCross2Var(divE1, divE2, false),
        model.getCross2Var(divP1, divP2, false)
    });
  }
};

struct CrossingTriple {
  int a1, a2;
  int b1, b2;
  int c1, c2;
  int64_t key;

  explicit CrossingTriple(int m,
                          int x1, int x2,
                          int y1, int y2,
                          int z1, int z2)
      : a1(x1), a2(x2), b1(y1), b2(y2), c1(z1), c2(z2) {
    // normalize each pair
    CHECK(0 <= a1 && a1 < a2 && a2 < m);
    CHECK(0 <= b1 && b1 < b2 && b2 < m);
    CHECK(0 <= c1 && c1 < c2 && c2 < m);

    // normalize order of the three pairs lexicographically: (a1,a2) < (b1,b2) < (c1,c2)
    auto P = [](int u1, int u2) { return std::make_pair(u1, u2); };
    if (P(a1, a2) > P(b1, b2)) { std::swap(a1, b1); std::swap(a2, b2); }
    if (P(b1, b2) > P(c1, c2)) { std::swap(b1, c1); std::swap(b2, c2); }
    if (P(a1, a2) > P(b1, b2)) { std::swap(a1, b1); std::swap(a2, b2); }
    CHECK(P(a1, a2) < P(b1, b2) && P(b1, b2) < P(c1, c2));

    // pack 6 edge-indices base-m
    key = (((((int64_t)a1 * m + a2) * m + b1) * m + b2) * m + c1) * m + c2;
  }

  void encode(SATModel& model, const InputGraph& graph) const {
    const int divA1 = a1 + graph.n;
    const int divA2 = a2 + graph.n;
    CHECK(canBeMerged(divA1, divA2, graph));
    const int divB1 = b1 + graph.n;
    const int divB2 = b2 + graph.n;
    CHECK(canBeMerged(divB1, divB2, graph));
    const int divC1 = c1 + graph.n;
    const int divC2 = c2 + graph.n;
    CHECK(canBeMerged(divC1, divC2, graph));

    model.addClause({
        model.getCross2Var(divA1, divA2, false),
        model.getCross2Var(divB1, divB2, false),
        model.getCross2Var(divC1, divC2, false),
    });
  }
};

/// TODO: merge with similar impl
struct SepCycleTraversal {
  SepCycleTraversal(const InputGraph& graph, const Params& params)
    : graph(graph), n(graph.n), m((int)graph.edges.size()), verbose(params.verbose)
    {}

  bool init(int numPairs) {
    for (int e1 = 0; e1 < m; e1++) {
      for (int e2 = e1 + 1; e2 < m; e2++) {
        if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
          continue;
        possibleCrossings.push_back({e1, e2});
      }
    }
    LOG_IF(verbose, "encoding new-sep constraints for %d pairs; |possible_crossings| = %d", 
           numPairs, possibleCrossings.size());
    // TODO: maybe skip if possibleCrossings is large?

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
  std::vector<CrossingPair> build2Clauses() {
    std::vector<CrossingPair> clauses2;

    // the main crossing
    for (size_t idx_uv = 0; idx_uv < possibleCrossings.size(); idx_uv++) {
      const int e1 = possibleCrossings[idx_uv].first;
      const int e2 = possibleCrossings[idx_uv].second;
      const auto& [u, v] = graph.edges[e1];
      const auto& [x, y] = graph.edges[e2];
      CHECK(u < v && x < y);
      CHECK(all_unique({x, y, u, v}));

      // circular order: u--x--v--y
      markCrossing(u, v, x, y);

      // secondary crossings
      for (size_t cr0 = 0; cr0 < possibleCrossings.size(); cr0++) {
        const int e1_cr0 = possibleCrossings[cr0].first;
        const int e2_cr0 = possibleCrossings[cr0].second;
        if (idx_uv == cr0)
          continue;
        // circular order: u1--u2--v1--v2
        const auto& [u1, v1] = graph.edges[e1_cr0];
        const auto& [u2, v2] = graph.edges[e2_cr0];
        CHECK(u1 < v1 && u2 < v2);

        // stop early, if the crossing pair is already processed
        const CrossingPair crossPair(m, e1, e2, e1_cr0, e2_cr0);
        if (isForbiddenCrossingPair(crossPair))
          continue;

        // check that the selected crossings are compatibe
        if (!canMarkCrossing(u1, v1, u2, v2))
          continue;

        markCrossing(u1, v1, u2, v2);

        // need to search cycles both via (x, y) and via (u, v)
        CHECK(all_unique({e1, e2, e1_cr0, e2_cr0}));
        std::vector<std::pair<int, int>> takenCrossings = {{e1, e2}, {e1_cr0, e2_cr0}};
        if (findSepCycle(x, y, u, v, takenCrossings) || findSepCycle(u, v, x, y, takenCrossings)) {
          // save the 2-clause
          clauses2.push_back(crossPair);
          addForbiddenCrossingPair(crossPair);
        }

        unmarkCrossing(u1, v1, u2, v2);
      }

      unmarkCrossing(u, v, x, y);
    }

    return clauses2;
  }

  /// Find triples of crossings that cannot happen in a 1-planar drawing
  std::vector<CrossingTriple> build3Clauses() {
    for (const auto& [u, v] : graph.edges) {
      CHECK(isEdge[u][v]);
    }
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        CHECK(!isCrossEdge[i][j]);
        CHECK(!isKiteEdge[i][j]);
      }
    }

    std::vector<CrossingTriple> clauses3;

    // the main crossing
    for (size_t idx_uv = 0; idx_uv < possibleCrossings.size(); idx_uv++) {
      const int e1 = possibleCrossings[idx_uv].first;
      const int e2 = possibleCrossings[idx_uv].second;
      const auto& [u, v] = graph.edges[e1];
      const auto& [x, y] = graph.edges[e2];
      CHECK(u < v && x < y);
      CHECK(all_unique({x, y, u, v}));

      markCrossing(u, v, x, y);

      // secondary crossing-A
      for (size_t cr0 = 0; cr0 < possibleCrossings.size(); cr0++) {
        const int e1_cr0 = possibleCrossings[cr0].first;
        const int e2_cr0 = possibleCrossings[cr0].second;
        if (idx_uv == cr0) 
          continue;

        const auto& [u1, v1] = graph.edges[e1_cr0];
        const auto& [u2, v2] = graph.edges[e2_cr0];
        CHECK(u1 < v1 && u2 < v2);

        // skip the triple if this pair is already forbidden
        if (isForbiddenCrossingPair(CrossingPair(m, e1, e2, e1_cr0, e2_cr0))) 
          continue;

        if (!canMarkCrossing(u1, v1, u2, v2))
          continue;

        markCrossing(u1, v1, u2, v2);

        // secondary crossing-B
        for (size_t cr1 = cr0 + 1; cr1 < possibleCrossings.size(); cr1++) {
          const int e1_cr1 = possibleCrossings[cr1].first;
          const int e2_cr1 = possibleCrossings[cr1].second;
          if (idx_uv == cr1) 
            continue;

          const auto& [w1, z1] = graph.edges[e1_cr1];
          const auto& [w2, z2] = graph.edges[e2_cr1];
          CHECK(w1 < z1 && w2 < z2);

          // skip the triple if these pairs are already forbidden
          if (isForbiddenCrossingPair(CrossingPair(m, e1, e2, e1_cr1, e2_cr1))) 
            continue;
          if (isForbiddenCrossingPair(CrossingPair(m, e1_cr0, e2_cr0, e1_cr1, e2_cr1))) 
            continue;

          // compatibility checks with already marked stuff
          if (!canMarkCrossing(w1, z1, w2, z2))
            continue;

          const CrossingTriple triple(m, e1, e2, e1_cr0, e2_cr0, e1_cr1, e2_cr1);
          if (isForbiddenCrossingTriple(triple))
            continue;

          markCrossing(w1, z1, w2, z2);

          // need to search cycles both via (x, y) and via (u, v)
          CHECK(all_unique({e1, e2, e1_cr0, e2_cr0, e1_cr1, e2_cr1}));
          const std::vector<std::pair<int,int>> takenCrossings = {{e1,e2}, {e1_cr0,e2_cr0}, {e1_cr1,e2_cr1}};
          if (findSepCycle(x, y, u, v, takenCrossings) || findSepCycle(u, v, x, y, takenCrossings)) {
            // save the 3-clause
            clauses3.push_back(triple);
            addForbiddenCrossingTriple(triple);
          }

          unmarkCrossing(w1, z1, w2, z2);
        }

        unmarkCrossing(u1, v1, u2, v2);
      }

      unmarkCrossing(u, v, x, y);
    }

    return clauses3;
  }


private:
  void addForbiddenCrossingPair(const CrossingPair& crossPair) {
    forbiddenCrossingPairs.insert(crossPair.key);
  }

  bool isForbiddenCrossingPair(const CrossingPair& crossPair) const {
    return forbiddenCrossingPairs.find(crossPair.key) != forbiddenCrossingPairs.end();
  }

  void addForbiddenCrossingTriple(const CrossingTriple& triple) {
    forbiddenCrossingTriples.insert(triple.key);
  }

  bool isForbiddenCrossingTriple(const CrossingTriple& triple) const {
    return forbiddenCrossingTriples.find(triple.key) != forbiddenCrossingTriples.end();
  }

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

  /// Find all separating cycles: paths from x to y avoiding vertices {x, y, u, v}
  /// A path may contain as many kite/cross edges as possible but only one free edge
  /// If a path contains one crossed edge, then disjoint paths cannot use the second one
  bool findSepCycle(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings, bool Debug=false) {
    // First, check cycle (x -- c -- y)
    if (findSepCycle1(x, y, u, v, takenCrossings, Debug))
      return true;

    // Second, check cycle (x -- c1 -- c2 -- y)
    if (findSepCycle2(x, y, u, v, takenCrossings, Debug)) {
      return true;
    }

    return false;
  }

  /// Find a separating cycle with one intermediate vertex
  bool findSepCycle1(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings, bool Debug=false) {
    // iterate over "neighbors" (including kite) of x
    for (const int c: adjList[x]) {
      if (contains({x, y, u, v}, c))
        continue;

      CHECK(isEdge[x][c] || isKiteEdge[x][c]);
      if (!isEdge[c][y] && !isKiteEdge[c][y])
        continue;

      const int numFree = int(!isCrossEdge[x][c] && !isKiteEdge[x][c]) + 
                          int(!isCrossEdge[c][y] && !isKiteEdge[c][y]);
      if (numFree >= std::min(graph.degree(u) - 1, graph.degree(v) - 1))
        continue;

      // Check if there are two edge-disjoint paths from u to v disjoint from (x, y, c); one extra path is via edge (u, v)
      const auto forbiddenEdges = getForbiddenCycleEdges({x, c, y}, takenCrossings);
      const int numPaths = countEdgeDisjointPaths(u, v, graph.adj, {x, y, c}, forbiddenEdges, numFree + 1);
      if (numPaths > numFree) {
        if (Debug) {
          LOG("found sep cycle: u = %d; v = %d; x = %d; y = %d; c = %d", u, v, x, y, c);
        }
        return true;
      }
    }

    return false;
  }

  /// Find a separating cycle with two intermediate vertices
  bool findSepCycle2(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings, bool Debug=false) {
    // iterate over "neighbors" (including kite) of x
    for (const int c1: adjList[x]) {
      if (contains({x, y, u, v}, c1))
        continue;
      CHECK(isEdge[x][c1] || isKiteEdge[x][c1]);

      // TODO: iterate over "neighbors" (including kite) of c1
      for (const int c2: adjList[c1]) {
        if (contains({x, y, u, v, c1}, c2))
          continue;
        CHECK(all_unique({x, y, u, v, c1, c2}));
        CHECK(isEdge[c1][c2] || isKiteEdge[c1][c2]);
        if (!isEdge[c2][y] && !isKiteEdge[c2][y])
          continue;

        const int numFree = int(!isCrossEdge[x][c1] && !isKiteEdge[x][c1]) + 
                            int(!isCrossEdge[c1][c2] && !isKiteEdge[c1][c2]) +
                            int(!isCrossEdge[c2][y] && !isKiteEdge[c2][y]);
        if (numFree >= std::min(graph.degree(u) - 1, graph.degree(v) - 1))
          continue;

        // Check if there are two edge-disjoint paths from u to v disjoint from (x, y, c); one extra path is via edge (u, v)
        const auto forbiddenEdges = getForbiddenCycleEdges({x, c1, c2, y}, takenCrossings);
        const int numPaths = countEdgeDisjointPaths(u, v, graph.adj, {x, y, c1, c2}, forbiddenEdges, numFree + 1);
        if (numPaths > numFree) {
          if (Debug) {
            LOG("found sep cycle: u = %d; v = %d; x = %d; y = %d; c1 = %d; c2 = %d", u, v, x, y, c1, c2);
          }
          return true;
        }
      }
    }

    return false;
  }

  std::vector<EdgeTy> getForbiddenCycleEdges(
      const std::vector<int>& cycle,
      const std::vector<std::pair<int, int>> &takenCrossings) {
    std::vector<EdgeTy> forbiddenEdges;
    for (size_t i = 0; i < cycle.size(); i++) {
      const EdgeTy cycleEdge = make_edge(cycle[i], cycle[(i + 1) % cycle.size()]);
      for (const auto& [e1, e2] : takenCrossings) {
        if (graph.edges[e1] == cycleEdge) {
          forbiddenEdges.push_back(graph.edges[e2]);
        }
        if (graph.edges[e2] == cycleEdge) {
          forbiddenEdges.push_back(graph.edges[e1]);
        }
      }
    }
    return forbiddenEdges;
  }

private:
  const InputGraph& graph;
  const int n;
  const int m;
  const int verbose;

  std::unordered_set<int64_t> forbiddenCrossingPairs;
  std::unordered_set<int64_t> forbiddenCrossingTriples;

  std::vector<std::pair<int, int>> possibleCrossings;
  std::vector<std::vector<bool>> isEdge;
  std::vector<std::vector<int>> adjList;

  std::vector<std::vector<bool>> isCrossEdge;
  std::vector<std::vector<int>> isKiteEdge;
};

/// Separating cycles
void encodeSepCyclesConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  const int verbose = params.verbose;
  const int numPairs = to_int(params.sepCycleConstraints);
  CHECK(2 <= numPairs && numPairs <= 3, "incorrect value of sepCycleConstraints");

  SepCycleTraversal sct(graph, params);
  if (sct.init(numPairs)) {
    // forbidden pairs (2-clauses)
    auto clauses2 = sct.build2Clauses();
    for (const auto& clause : clauses2) {
      clause.encode(model, graph);
    }
    LOG_IF(verbose, "  added %'7d unique sep-cycle 2-clauses", clauses2.size());

    if (numPairs == 3) {
      // forbidden triples (3-clauses)
      auto clauses3 = sct.build3Clauses();
      for (const auto& clause : clauses3) {
        clause.encode(model, graph);
      }
      LOG_IF(verbose, "  added %'7d unique sep-cycle 3-clauses", clauses3.size());
    }

    // CHECK(params.custom != "");
    // const int clauseIndex = to_int(params.custom);
    // CHECK(0 <= clauseIndex && clauseIndex < (int)newClauses.size(), "incorrect custom index");
    // LOG_IF(verbose, "  force custom-pair crossings: %d out of %d", clauseIndex, newClauses.size());
    // // should be UNSAT if we force the two pairs to cross
    // const auto& [d1, d2, d3, d4] = newClauses[clauseIndex];
    // LOG_IF(verbose, "    d1 = %d; d2 = %d; d3 = %d; d4 = %d", d1, d2, d3, d4);
    // model.addClause(model.getCross2Var(d1, d2, true));
    // model.addClause(model.getCross2Var(d3, d4, true));
  }

}
