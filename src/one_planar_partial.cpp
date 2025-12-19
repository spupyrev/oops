#include "graph_algorithms.h"
#include "logging.h"
#include "one_planar.h"

using namespace std;

/// TODO: merge with similar impl
struct PartialFinder {
  PartialFinder(const InputGraph& graph, const std::vector<std::pair<int, int>>& possibleCrossings, const int verbose) 
    : graph(graph),
      possibleCrossings(possibleCrossings),
      verbose(verbose)
    {}

  std::vector<std::vector<std::pair<int, int>>> search(int maxNumCrossings) {
    takenCrossings.clear();
    foundConstraints.clear();
    isFreeEdge = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, 0));
    isCrossedEdge = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, 0));
    takenVertexDegree = std::vector<int>(graph.n, 0);
    index = std::vector<int>(graph.n, 0);

    searchRec(0, maxNumCrossings);

    return foundConstraints;
  }

private:
  void searchRec(size_t curIdx, int remainingCrossings) {
    if (remainingCrossings >= 0) {
      if (!takenCrossings.empty() && !isPartialPlanar()) {
        foundConstraints.push_back(takenCrossings);
        LOG_IF(verbose && foundConstraints.size() % 5000 == 0, "  found %d partial constraints so far...", foundConstraints.size());
        // print 
        LOG_IF(verbose >= 2, "found %d-th partial constraint:", foundConstraints.size());
        for (const auto& [e1, e2] : takenCrossings) {
          LOG_IF(verbose >= 2, "  (%d, %d) -- (%d, %d);   e1=%d, e2=%d", 
                 graph.edges[e1].first, graph.edges[e1].second, graph.edges[e2].first, graph.edges[e2].second, e1, e2);
        }
        return;
      }
    }
    if (remainingCrossings == 0) {
      return;
    }
    CHECK(possibleCrossings.size() >= curIdx);
    if (curIdx == possibleCrossings.size()) {
      return;
    }

    // try to get the current
    const int e1 = possibleCrossings[curIdx].first;
    const int e2 = possibleCrossings[curIdx].second;
    const auto& [u1, v1] = graph.edges[e1];
    const auto& [u2, v2] = graph.edges[e2];
    // circular order: u1--u2--v1--v2

    if (isCrossedEdge[u1][v1] == 0 && isCrossedEdge[u2][v2] == 0) {
      if (isFreeEdge[u1][v1] == 0 && isFreeEdge[u2][v2] == 0) {
        if (isCrossedEdge[u1][u2] == 0 && isCrossedEdge[u2][v1] == 0 && isCrossedEdge[v1][v2] == 0 && isCrossedEdge[v2][u1] == 0) {
          if (takenCrossings.empty() || takenVertexDegree[u1] || takenVertexDegree[v1] || takenVertexDegree[u2] || takenVertexDegree[v2]) {
            // register the crossings
            addCrossedEdge(u1, v1, 1);
            addCrossedEdge(u2, v2, 1);

            addFreeEdge(u1, u2, 1);
            addFreeEdge(u2, v1, 1);
            addFreeEdge(v1, v2, 1);
            addFreeEdge(v2, u1, 1);

            takenVertexDegree[u1]++;
            takenVertexDegree[v1]++;
            takenVertexDegree[u2]++;
            takenVertexDegree[v2]++;

            takenCrossings.push_back(possibleCrossings[curIdx]);

            // recurse
            searchRec(curIdx + 1, remainingCrossings - 1);

            // rollback
            addCrossedEdge(u1, v1, -1);
            addCrossedEdge(u2, v2, -1);

            addFreeEdge(u1, u2, -1);
            addFreeEdge(u2, v1, -1);
            addFreeEdge(v1, v2, -1);
            addFreeEdge(v2, u1, -1);

            takenVertexDegree[u1]--;
            takenVertexDegree[v1]--;
            takenVertexDegree[u2]--;
            takenVertexDegree[v2]--;

            takenCrossings.pop_back();
          }
        }
      }
    }

    // skip the current
    searchRec(curIdx + 1, remainingCrossings);
  }

  void addFreeEdge(const int u, const int v, const int delta) {
    isFreeEdge[u][v] += delta; 
    isFreeEdge[v][u] += delta;
  }

  void addCrossedEdge(const int u, const int v, const int delta) {
    isCrossedEdge[u][v] += delta; 
    isCrossedEdge[v][u] += delta;
  }

  // check if takenCrossings can be extended to a planar graph
  bool isPartialPlanar() const {
    CHECK(takenCrossings.size() >= 1);
    const int n = graph.n;

    // Stop early: there cannot be (orig) vertices that belong to a single crossing
    for (int i = 0; i < n; i++) {
      if (takenVertexDegree[i] == 1)
        return true;
    }

    // Create a partial graph
    std::fill(index.begin(), index.end(), -1);
    int N = 0;
    for (int i = 0; i < n; i++) {
      if (takenVertexDegree[i]) {
        index[i] = N;
        N++;
      }
    }
    takenEdges.clear();
    for (const auto& [e1, e2] : takenCrossings) {
      const auto& [u1, v1] = graph.edges[e1];
      const auto& [u2, v2] = graph.edges[e2];
      CHECK(index[u1] != -1 && index[v1] != -1 && index[u2] != -1 && index[v2] != -1);
      const int d = N;
      N++;

      takenEdges.push_back(make_edge(index[u1], d));
      takenEdges.push_back(make_edge(index[v1], d));
      takenEdges.push_back(make_edge(index[u2], d));
      takenEdges.push_back(make_edge(index[v2], d));

      takenEdges.push_back(make_edge(index[u1], index[u2]));
      takenEdges.push_back(make_edge(index[u2], index[v1]));
      takenEdges.push_back(make_edge(index[v1], index[v2]));
      takenEdges.push_back(make_edge(index[v2], index[u1]));
    }
    sort_unique(takenEdges);

    numCheckedPlanar++;
    const bool res = isPlanar(N, takenEdges, verbose - 1);
    // if (!res) {
    //   // check degrees of orig vertices
    //   std::vector<int> degree(N);
    //   for (const auto& [s, t] : takenEdges) {
    //     degree[s]++;
    //     degree[t]++;
    //   }
    //   for (int i = 0; i + takenCrossings.size() < N; i++) {
    //     CHECK(degree[i] >= 4);
    //   }
    // }
    // if (!res && N >= 10) {
    //   LOG("non-planar with n = %d, m = %d", N, takenEdges.size());
    // }
    return res;
  }

public:
  mutable int numCheckedPlanar = 0;  

private:
  const InputGraph& graph;
  const std::vector<std::pair<int, int>>& possibleCrossings;
  const int verbose;

  std::vector<std::pair<int, int>> takenCrossings;
  std::vector<std::vector<std::pair<int, int>>> foundConstraints;
  std::vector<std::vector<int>> isFreeEdge;
  std::vector<std::vector<int>> isCrossedEdge;
  std::vector<int> takenVertexDegree;

  mutable std::vector<EdgeTy> takenEdges;
  mutable std::vector<int> index;

};

/// Add constraints based on partial planarity
void encodePartialConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  CHECK(!graph.isDirected());

  const int numPairs = to_int(params.partialConstraints);
  CHECK(2 <= numPairs && numPairs <= 3, "incorrect value of partialConstraints");

  const int n = graph.n;
  const int m = (int)graph.edges.size();
  std::vector<std::pair<int, int>> possibleCrossings;
  for (int e1 = 0; e1 < m; e1++) {
    for (int e2 = e1 + 1; e2 < m; e2++) {
      if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
        continue;
      possibleCrossings.push_back({e1, e2});
    }
  }
  if (possibleCrossings.size() > 4096) {
    LOG_IF(params.verbose, "skipped partial constraints due to too many possible_crossings (%d)", possibleCrossings.size());
    return;
  }
  LOG_IF(params.verbose, "encoding partial constraints for %d pairs; |possible_crossings| = %d", numPairs, possibleCrossings.size());

  PartialFinder finder(graph, possibleCrossings, params.verbose);
  auto constraints = finder.search(numPairs);
  LOG_IF(params.verbose && constraints.empty(), "  no partial constraints; checked planar: %d", finder.numCheckedPlanar);
  LOG_IF(params.verbose && !constraints.empty(), "  found %d partial constraints; checked planar: %d", constraints.size(), finder.numCheckedPlanar);

  std::vector<int> numConstraints(numPairs + 1, 0);
  for (const auto& constraint: constraints) {
    CHECK((int)constraint.size() <= numPairs);
    numConstraints[constraint.size()]++;
    MClause clause;
    for (const auto& [e1, e2] : constraint) {
      clause.addVar(model.getCross2Var(e1 + n, e2 + n, false));
    }
    model.addClause(clause);
  }
  for (size_t i = 0; i < numConstraints.size(); i++) {
    LOG_IF(params.verbose >= 2 && numConstraints[i] > 0, "    %'d %d-partial constraints", numConstraints[i], i);
  }
}

/*
std::vector<std::tuple<int, int, int, int>> findCase3(const InputGraph& graph, const Params& params) {
  std::vector<std::tuple<int, int, int, int>> clauses;
  const int n = graph.n;
  const auto& adj = graph.adj;  

  auto isCase3 = [&](int x, int y, int u, int v, int a, int b, int c) -> bool {
    CHECK(all_unique({x, y, u, v, a, b, c}));
    CHECK(u < v);
    const int divXY = graph.findDivIndex(x, y);
    const int divUV = graph.findDivIndex(u, v);
    CHECK(canBeMerged(divXY, divUV, graph));
    const int divAC = graph.findDivIndex(a, c);
    const int divYB = graph.findDivIndex(y, b);
    CHECK(canBeMerged(divAC, divYB, graph));

    // check if there are two edge-disjoint paths from u to v disjoint from (x, y, c)
    std::vector<std::pair<int, int>> forbiddenEdges = {{u, v}};
    const int numPaths = countEdgeDisjointPaths(u, v, graph.adj, {x, y, c}, forbiddenEdges);
    // CHECK(numPaths >= 1);
    if (numPaths >= 2) {
      // LOG_IF(params.verbose >= 3, "case3: x = %d, y = %d, c = %d, a = %d, b = %d, u = %d, v = %d", x, y, c, a, b, u, v);
      // LOG_IF(params.verbose >= 3, "  num_paths = %d", numPaths);
      return true;
    }

    return false;
  };

  // pick x
  for (int x = 0; x < n; x++) {
    // pick y and c
    for (int y : adj[x]) {
      for (int c : adj[x]) {
        if (y == c) continue;
        // pick a and b
        for (int b : adj[y]) {
          for (int a : adj[c]) {
            if (b == x || b == c) continue;
            if (a == x || a == y || a == b) continue;
            CHECK(all_unique({x, y, c, a, b}));
            const int divAC = graph.findDivIndex(a, c);
            const int divYB = graph.findDivIndex(y, b);
            if (!canBeMerged(divYB, divAC, graph)) 
              continue;

            // pick u and v
            for (int u = 0; u < n; u++) {
              for (int v : adj[u]) {
                if (u > v) continue;
                if (contains({x, y, a, b, c}, u)) continue;
                if (contains({x, y, a, b, c}, v)) continue;
                CHECK(all_unique({x, y, c, a, b, u, v}));

                const int divXY = graph.findDivIndex(x, y);
                const int divUV = graph.findDivIndex(u, v);
                if (!canBeMerged(divXY, divUV, graph)) 
                  continue;

                if (isCase3(x, y, u, v, a, b, c)) {
                  clauses.push_back(std::make_tuple(divXY, divUV, divAC, divYB));
                }
              }
            }
          }
        }
      }
    }
  }

  return clauses;
}
*/

int encode2Clauses(SATModel& model, const std::vector<std::tuple<int, int, int, int>> &clauses) {
  for (const auto& [a, b, c, d] : clauses) {
    model.addClause({
        model.getCross2Var(a, b, false), 
        model.getCross2Var(c, d, false)
    });
  }
  return (int)clauses.size();
}

std::vector<std::tuple<int, int, int, int>> makeUnique(std::vector<std::tuple<int, int, int, int>> &clauses) {
  std::vector<std::tuple<int, int, int, int>> newClauses;
  for (const auto& [a, b, c, d] : clauses) {
    auto m1 = std::minmax(a, b);
    auto m2 = std::minmax(c, d);
    if (m1 < m2) {
      newClauses.push_back(std::make_tuple(m1.first, m1.second, m2.first, m2.second));
    } else {
      newClauses.push_back(std::make_tuple(m2.first, m2.second, m1.first, m1.second));
    }
  }

  sort_unique(newClauses);
  return newClauses;
}

/// TODO: merge with similar impl
struct SepCycleTraversal {
  SepCycleTraversal(const InputGraph& graph, const Params& params) 
    : graph(graph), n(graph.n), verbose(params.verbose)
    {}

  bool init() {
    const int m = (int)graph.edges.size();

    for (int e1 = 0; e1 < m; e1++) {
      for (int e2 = e1 + 1; e2 < m; e2++) {
        // if (e1 == e2)
        //   continue;
        if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
          continue;
        possibleCrossings.push_back({e1, e2});
      }
    }
    LOG_IF(verbose, "encoding new-sep constraints; |possible_crossings| = %d", possibleCrossings.size());
    // TODO: maybe skip if possibleCrossings is large?

    isFreeEdge = std::vector<std::vector<bool>> (n, std::vector<bool>(n, 0));
    isCrossEdge = std::vector<std::vector<bool>>(n, std::vector<bool>(n, 0));
    isKiteEdge = std::vector<std::vector<int>>(n, std::vector<int>(n, 0));

    // init free edges
    for (const auto& [u, v] : graph.edges) {
      setEdge(isFreeEdge, u, v, true);
    }

    return true;
  }

  std::vector<std::tuple<int, int, int, int>> build2Clauses() {
    // the main crossing
    for (size_t idx_uv = 0; idx_uv < possibleCrossings.size(); idx_uv++) {
      const int e1 = possibleCrossings[idx_uv].first;
      const int e2 = possibleCrossings[idx_uv].second;
      const auto& [u, v] = graph.edges[e1];
      const auto& [x, y] = graph.edges[e2];
      CHECK(u < v && x < y);

      // circular order: u--x--v--y
      markCrossing(u, v, x, y);

      // secondary crossings
      for (size_t cr0 = 0; cr0 < possibleCrossings.size(); cr0++) {
        const int e1_cr0 = possibleCrossings[cr0].first;
        const int e2_cr0 = possibleCrossings[cr0].second;
        // circular order: u1--u2--v1--v2
        const auto& [u1, v1] = graph.edges[e1_cr0];
        const auto& [u2, v2] = graph.edges[e2_cr0];
        CHECK(u1 < v1 && u2 < v2);

        // check that the selected crossings are compatibe
        if (isCrossEdge[u1][v1] || isCrossEdge[u2][v2])
          continue;
        if (isKiteEdge[u1][v1] || isKiteEdge[u2][v2])
          continue;
        if (isCrossEdge[u1][u2] || isCrossEdge[u2][v1] || isCrossEdge[v1][v2] || isCrossEdge[v2][u1])
          continue;

        markCrossing(u1, v1, u2, v2);

        // need to search cycles both via (x, y) and via (u, v)
        std::vector<std::pair<int, int>> takenCrossings = {{e1, e2}, {e1_cr0, e2_cr0}};
        if (findSepCycle(x, y, u, v, takenCrossings) || findSepCycle(u, v, x, y, takenCrossings)) {
          // add the clause
          CHECK(all_unique({x, y, u, v}));
          CHECK(u < v && x < y);
          const int divXY = graph.findDivIndex(x, y);
          const int divUV = graph.findDivIndex(u, v);
          CHECK(canBeMerged(divXY, divUV, graph));
          const int divAC = graph.findDivIndex(u1, v1);
          const int divYB = graph.findDivIndex(u2, v2);
          CHECK(canBeMerged(divAC, divYB, graph));

          // clauses.push_back(std::make_tuple(divXY, divUV, divAC, divYB));
          auto m1 = std::minmax(divXY, divUV);
          auto m2 = std::minmax(divAC, divYB);
          if (m1 < m2) {
            clauses.push_back(std::make_tuple(m1.first, m1.second, m2.first, m2.second));
          } else {
            clauses.push_back(std::make_tuple(m2.first, m2.second, m1.first, m1.second));
          }

          // const auto last = clauses.back();
          // if (std::get<0>(last) == 9 && std::get<1>(last) == 14 && std::get<2>(last) == 10 && std::get<3>(last) == 18) {
          //   LOG("x = %d; y = %d; u = %d; v = %d; u1 = %d; v1 = %d; u2 = %d; v2 = %d", x, y, u, v, u1, v1, u2, v2);
          //   findSepCycle(x, y, u, v, takenCrossings, true);
          //   findSepCycle(u, v, x, y, takenCrossings, true);
          // }

        }

        unmarkCrossing(u1, v1, u2, v2);
      }

      unmarkCrossing(u, v, x, y);
    }

    return clauses;
  }

private:
  void setEdge(std::vector<std::vector<bool>> &isEdge, int u, int v, bool value) {
    isEdge[u][v] = value;
    isEdge[v][u] = value;
  }

  void setEdge(std::vector<std::vector<int>> &isEdge, int u, int v, int delta) {
    isEdge[u][v] += delta;
    isEdge[v][u] += delta;
  }

  void markCrossing(int u, int v, int x, int y) {
    // cross
    CHECK(!isCrossEdge[u][v] && !isCrossEdge[x][y]);
    CHECK(isFreeEdge[u][v] && isFreeEdge[x][y]);
    CHECK(isKiteEdge[u][v] == 0 && isKiteEdge[x][y] == 0);
    setEdge(isCrossEdge, u, v, true);
    setEdge(isCrossEdge, x, y, true);
    setEdge(isFreeEdge, u, v, false);
    setEdge(isFreeEdge, x, y, false);

    // kite: circular order is u--x--v--y
    CHECK(!isCrossEdge[u][x] && !isCrossEdge[x][v] && !isCrossEdge[v][y] && !isCrossEdge[y][u]);
    setEdge(isKiteEdge, u, x, 1);
    setEdge(isKiteEdge, x, v, 1);
    setEdge(isKiteEdge, v, y, 1);
    setEdge(isKiteEdge, y, u, 1);
  };
  
  void unmarkCrossing(int u, int v, int x, int y) {
    // cross
    CHECK(isCrossEdge[u][v] && isCrossEdge[x][y]);
    CHECK(!isFreeEdge[u][v] && !isFreeEdge[x][y]);
    CHECK(isKiteEdge[u][v] == 0 && isKiteEdge[x][y] == 0);
    setEdge(isCrossEdge, u, v, false);
    setEdge(isCrossEdge, x, y, false);
    setEdge(isFreeEdge, u, v, true);
    setEdge(isFreeEdge, x, y, true);

    // kite: circular order is u--x--v--y
    setEdge(isKiteEdge, u, x, -1);
    setEdge(isKiteEdge, x, v, -1);
    setEdge(isKiteEdge, v, y, -1);
    setEdge(isKiteEdge, y, u, -1);
  };

  bool findSepCycle(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings, bool Debug=false) {
    /// Find all separating cycles: paths from x to y avoiding vertices {x, y, u, v}
    /// A path may contain as many kite/cross edges as possible but only one free edge
    /// If a path contains one crossed edge, then disjoint paths cannot use the second one

    // For now, check a cycle (x -- c -- y)
    for (int c = 0; c < n; c++) {
      if (contains({x, y, u, v}, c))
        continue;
      if (!isFreeEdge[x][c] && !isCrossEdge[x][c] && !isKiteEdge[x][c])
        continue;
      if (!isFreeEdge[c][y] && !isCrossEdge[c][y] && !isKiteEdge[c][y])
        continue;
      
      int numFree = 0;
      int numCrossKite = 0;
      if (isFreeEdge[x][c])
        numFree++;
      else
        numCrossKite++;
      if (isFreeEdge[c][y])
        numFree++;
      else
        numCrossKite++;

      if (numFree > 1)
        continue;

      CHECK(numCrossKite >= 0 && numFree <= 1);

      // Check if there are two edge-disjoint paths from u to v disjoint from (x, y, c); one extra path is via edge (u, v)
      const auto forbiddenEdges = getForbiddenEdges({x, y, c}, takenCrossings);
      const int numPaths = countEdgeDisjointPaths(u, v, graph.adj, {x, y, c}, forbiddenEdges);
      //CHECK(numPaths >= 1);
      if (numPaths >= 2) {
        if (Debug) {
          LOG("found sep cycle: u = %d; v = %d; x = %d; y = %d; c = %d", u, v, x, y, c);
        }
        return true;
      }
    }

    return false;
  }

  std::vector<EdgeTy> getForbiddenEdges(
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
  const int verbose;

  std::vector<std::pair<int, int>> possibleCrossings;
  std::vector<std::tuple<int, int, int, int>> clauses;

  std::vector<std::vector<bool>> isFreeEdge;
  std::vector<std::vector<bool>> isCrossEdge;
  std::vector<std::vector<int>> isKiteEdge;
};

/// Separating cycles
void encodeSepCyclesConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  const int verbose = params.verbose;

  // auto case3Clauses = findCase3(graph, params);
  // auto reducedCase3Clauses = makeUnique(case3Clauses);
  // LOG_IF(verbose, "reduced case-3 clauses from %d to %d", case3Clauses.size(), reducedCase3Clauses.size());
  // const int numCase3 = encode2Clauses(model, reducedCase3Clauses);
  // LOG_IF(verbose, "added %d sep-cycle constraints (%d case3)", numCase3, numCase3);

  SepCycleTraversal sct(graph, params);
  if (sct.init()) {
    auto clauses = sct.build2Clauses();
    LOG_IF(verbose, "  found %d sep-cycle constraints", clauses.size());
    auto reducedClauses = makeUnique(clauses);
    LOG_IF(verbose, "  reduced sep-cycle clauses from %d to %d", clauses.size(), reducedClauses.size());
    encode2Clauses(model, reducedClauses);
    LOG_IF(verbose, "  added %'d separating-cycle clauses", reducedClauses.size());

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
