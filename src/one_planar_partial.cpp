#include "graph_algorithms.h"
#include "logging.h"
#include "one_planar.h"

using namespace std;

/// TODO
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

int encodeCase3(SATModel& model, const InputGraph& graph, const Params& params) {
  const int verbose = params.verbose;
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
    const int numPaths = countEdgeDisjointPaths(u, v, graph.adj, {x, y, c});
    CHECK(numPaths >= 1);
    if (numPaths >= 3) {
      LOG_IF(verbose >= 2, "case3: x = %d, y = %d, c = %d, a = %d, b = %d, u = %d, v = %d", x, y, c, a, b, u, v);
      LOG_IF(verbose >= 2, "  num_paths = %d", numPaths);
      return true;
    }

    return false;
  };

  int numCase3 = 0;
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
                  model.addClause({
                      model.getCross2Var(divXY, divUV, false), 
                      model.getCross2Var(divAC, divYB, false)
                  });
                  numCase3++;
                  // auto m1 = std::minmax(divXY, divUV);
                  // auto m2 = std::minmax(divAC, divYB);
                  // LOG_IF(params.verbose >= 1, "case3: (%d, %d) -- (%d, %d)", m1.first, m1.second, m2.first, m2.second);
                }
              }
            }
          }
        }
      }
    }
  }
  return numCase3;
}

/// Separating cycles
void encodeSepCyclesConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  const int verbose = params.verbose;
  // const int n = graph.n;
  // const auto& adj = graph.adj;  

  int numCase3 = encodeCase3(model, graph, params);

  // ad-hoc

  LOG_IF(verbose, "added %d sep-cycle constraints (%d case3)", numCase3, numCase3);
}
