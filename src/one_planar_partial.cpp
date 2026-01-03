#include "graph_algorithms.h"
#include "logging.h"
#include "one_planar.h"

#include <utility>

using namespace std;

/// TODO: merge with similar impl
struct PartialTraversal {
  PartialTraversal(const InputGraph& graph, const Params& params, ForbiddenTuples& tuples)
    : graph(graph), n(graph.n), m((int)graph.edges.size()), verbose(params.verbose), forbiddenTuples(tuples)
    {}

  bool init(int numPairs) {
    for (int e1 = 0; e1 < m; e1++) {
      for (int e2 = e1 + 1; e2 < m; e2++) {
        if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
          continue;
        possibleCrossings.push_back({e1, e2});
      }
    }
    if (possibleCrossings.size() > 4096) {
      LOG_IF(verbose, "skipped partial constraints due to too many possible_crossings (%d)", possibleCrossings.size());
      return false;
    }
    LOG_IF(verbose, "partial constraints for %d pairs; |possible_crossings| = %d", numPairs, possibleCrossings.size());

    takenCrossings.clear();
    isFreeEdge = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, 0));
    isCrossedEdge = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, 0));
    takenVertexDegree = std::vector<int>(graph.n, 0);
    index = std::vector<int>(graph.n, 0);

    return true;
  }

  std::pair<int, int> search(int numPairs) {
    searchRec(0, numPairs);
    return {numConstraints2, numConstraints3};
  }

private:
  void searchRec(size_t curIdx, int remainingCrossings) {
    if (remainingCrossings >= 0) {
      if (!takenCrossings.empty() && !isPartialPlanar()) {
        processConstraint();
        return;
      }
    }
    if (remainingCrossings == 0)
      return;
    CHECK(possibleCrossings.size() >= curIdx);
    if (curIdx == possibleCrossings.size())
      return;

    // try to get the current
    const int e1 = possibleCrossings[curIdx].first;
    const int e2 = possibleCrossings[curIdx].second;
    const auto [u1, v1] = graph.edges[e1];
    const auto [u2, v2] = graph.edges[e2];
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

  void processConstraint() {
    CHECK(takenCrossings.size() == 2 || takenCrossings.size() == 3);
    bool added = false;
    if (takenCrossings.size() == 2) {
      const CrossingPair crossPair(m, takenCrossings[0].first, takenCrossings[0].second, 
                                      takenCrossings[1].first, takenCrossings[1].second);
      if (!forbiddenTuples.contains(crossPair)) {
        forbiddenTuples.insert(crossPair);
        numConstraints2++;
        added = true;
      }
    } else {
      const CrossingTriple crossTriple(m, takenCrossings[0].first, takenCrossings[0].second, 
                                          takenCrossings[1].first, takenCrossings[1].second,
                                          takenCrossings[2].first, takenCrossings[2].second);
      if (!forbiddenTuples.contains(crossTriple)) {
        forbiddenTuples.insert(crossTriple);
        numConstraints3++;
        added = true;
      }
    }

    LOG_IF(verbose && added && (numConstraints2 + numConstraints3) % 5000 == 0, 
           "  found %d partial constraints so far...", numConstraints2 + numConstraints3);
  }

  void addFreeEdge(const int u, const int v, const int delta) {
    isFreeEdge[u][v] += delta;
    isFreeEdge[v][u] += delta;
  }

  void addCrossedEdge(const int u, const int v, const int delta) {
    isCrossedEdge[u][v] += delta;
    isCrossedEdge[v][u] += delta;
  }

  /// Check if takenCrossings can be extended to a planar graph
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
      const auto [u1, v1] = graph.edges[e1];
      const auto [u2, v2] = graph.edges[e2];
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
    return isPlanar(N, takenEdges, verbose - 1);
  }

public:
  mutable int numCheckedPlanar = 0;

private:
  const InputGraph& graph;
  const int n;
  const int m;
  const int verbose;
  ForbiddenTuples& forbiddenTuples;
  std::vector<std::pair<int, int>> possibleCrossings;

  std::vector<std::pair<int, int>> takenCrossings;
  size_t numConstraints2 = 0;
  size_t numConstraints3 = 0;
  std::vector<std::vector<int>> isFreeEdge;
  std::vector<std::vector<int>> isCrossedEdge;
  std::vector<int> takenVertexDegree;

  mutable std::vector<EdgeTy> takenEdges;
  mutable std::vector<int> index;
};

/// Add constraints based on partial planarity
void encodePartialConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  // TODO: need this??
  CHECK(!graph.isDirected());

  const int verbose = params.verbose;
  const int numPairs = to_int(params.partialConstraints);
  CHECK(2 <= numPairs && numPairs <= 3, "incorrect value of partialConstraints");

  PartialTraversal finder(graph, params, model.getForbiddenTuples());
  if (!finder.init(numPairs)) {
    return;
  }

  const auto [numClauses2, numClauses3] = finder.search(numPairs);
  LOG_IF(verbose && numClauses2 + numClauses3 == 0, "  no partial constraints; checked planar: %d", finder.numCheckedPlanar);
  LOG_IF(verbose && numClauses2 > 0, "  found %'9d 2-clauses", numClauses2);
  LOG_IF(verbose && numClauses3 > 0, "  found %'9d 3-clauses", numClauses3);
}
