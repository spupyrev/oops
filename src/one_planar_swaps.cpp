#include "one_planar.h"
#include "logging.h"

using namespace std;

/// TODO: merge with similar impl
struct SwapTraversal {
  SwapTraversal(const InputGraph& graph, const Params& params, ForbiddenTuples& tuples)
    : graph(graph), n(graph.n), m((int)graph.edges.size()), verbose(params.verbose), forbiddenTuples(tuples)
    {}

  bool init(int numPairs, int maxNumSwaps) {
    for (int e1 = 0; e1 < m; e1++) {
      for (int e2 = e1 + 1; e2 < m; e2++) {
        if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
          continue;
        possibleCrossings.push_back({e1, e2});
      }
    }
    LOG_IF(verbose, "swap constraints for %d pairs and %d swaps; |possible_crossings| = %d", 
           numPairs, maxNumSwaps, possibleCrossings.size());

    isFreeEdge = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, 0));
    isCrossedEdge = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, 0));
    takenCrossings.clear();
    cost = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, -1));
    numSwaps = maxNumSwaps;

    return true;
  }

  size_t search2() {
    searchRec(0, 2);
    return numConstraints2;
  }

  size_t search3() {
    if (possibleCrossings.size() > 1280) {
      LOG_IF(verbose, "skipped building swap 3-clauses due to too many possible_crossings (%d)", possibleCrossings.size());
      return 0;
    } 
      
    searchRec(0, 3);
    return numConstraints3;
  }

private:
  void searchRec(size_t curIdx, int remainingCrossings) {
    CHECK(remainingCrossings >= 0);
    CHECK(curIdx <= possibleCrossings.size());

    if (remainingCrossings == 0) {
      // Approximation: there should be some overlap in taken vertices
      if (takenCrossings.size() == 3) {
        const auto [v1, u1] = graph.edges[takenCrossings[0].first];
        const auto [v2, u2] = graph.edges[takenCrossings[0].second];
        const auto [v3, u3] = graph.edges[takenCrossings[1].first];
        const auto [v4, u4] = graph.edges[takenCrossings[1].second];
        const auto [v5, u5] = graph.edges[takenCrossings[2].first];
        const auto [v6, u6] = graph.edges[takenCrossings[2].second];
        if (unique_size({v1, u1, v2, u2, v3, u3, v4, u4, v5, u5, v6, u6}) >= 11)
          return;
      }
      if (takenCrossings.size() >= 2 && canSwap()) {
        processConstraint();
      }
      return;
    }

    if (curIdx == possibleCrossings.size())
      return;

    // try to get the current
    const int e1 = possibleCrossings[curIdx].first;
    const int e2 = possibleCrossings[curIdx].second;
    const auto [u1, v1] = graph.edges[e1];
    const auto [u2, v2] = graph.edges[e2];
    // circular order: u1--u2--v1--v2

    if (canMarkCrossing(u1, v1, u2, v2)) {
      takenCrossings.push_back(possibleCrossings[curIdx]);

      // check if this crossing tuple is forbidden
      const bool isForbidden = takenCrossings.size() > 1 && isTakenForbidden();

      // recurse
      if (!isForbidden) {
        markCrossing(u1, v1, u2, v2);
        searchRec(curIdx + 1, remainingCrossings - 1);
        unmarkCrossing(u1, v1, u2, v2);
      }

      takenCrossings.pop_back();
    }

    // skip the current
    searchRec(curIdx + 1, remainingCrossings);
  }

  void processConstraint() {
    CHECK(takenCrossings.size() == 2 || takenCrossings.size() == 3);
    if (takenCrossings.size() == 2) {
      const CrossingPair crossPair(m, takenCrossings[0].first, takenCrossings[0].second, 
                                      takenCrossings[1].first, takenCrossings[1].second);
      CHECK(!forbiddenTuples.contains(crossPair));
      forbiddenTuples.insert(crossPair);
      numConstraints2++;
    } else {
      const CrossingTriple crossTriple(m, takenCrossings[0].first, takenCrossings[0].second, 
                                          takenCrossings[1].first, takenCrossings[1].second,
                                          takenCrossings[2].first, takenCrossings[2].second);
      CHECK(!forbiddenTuples.contains(crossTriple));
      forbiddenTuples.insert(crossTriple);
      numConstraints3++;
    }

    LOG_IF(verbose && (numConstraints2 + numConstraints3) % 50000 == 0, 
           "  found %d swap constraints so far...", numConstraints2 + numConstraints3);
  }

  bool isTakenForbidden() const {
    CHECK(takenCrossings.size() == 2 || takenCrossings.size() == 3);
    if (takenCrossings.size() == 2) {
      const CrossingPair crossPair(m, takenCrossings[0].first, takenCrossings[0].second, 
                                      takenCrossings[1].first, takenCrossings[1].second);
      return forbiddenTuples.contains(crossPair);
    } 

    const CrossingTriple crossTriple(m, takenCrossings[0].first, takenCrossings[0].second, 
                                        takenCrossings[1].first, takenCrossings[1].second,
                                        takenCrossings[2].first, takenCrossings[2].second);
    return forbiddenTuples.contains(crossTriple);
  }

  bool canMarkCrossing(int u1, int v1, int u2, int v2) const {
    if (isCrossedEdge[u1][v1] == 0 && isCrossedEdge[u2][v2] == 0) {
      if (isFreeEdge[u1][v1] == 0 && isFreeEdge[u2][v2] == 0) {
        if (isCrossedEdge[u1][u2] == 0 && isCrossedEdge[u2][v1] == 0 && isCrossedEdge[v1][v2] == 0 && isCrossedEdge[v2][u1] == 0) {
          return true;
        }
      }
    }
    return false;
  }

  void markCrossing(int u1, int v1, int u2, int v2) {
    addCrossedEdge(u1, v1, 1);
    addCrossedEdge(u2, v2, 1);

    addFreeEdge(u1, u2, 1);
    addFreeEdge(u2, v1, 1);
    addFreeEdge(v1, v2, 1);
    addFreeEdge(v2, u1, 1);
  }

  void unmarkCrossing(int u1, int v1, int u2, int v2) {
    addCrossedEdge(u1, v1, -1);
    addCrossedEdge(u2, v2, -1);

    addFreeEdge(u1, u2, -1);
    addFreeEdge(u2, v1, -1);
    addFreeEdge(v1, v2, -1);
    addFreeEdge(v2, u1, -1);
  }

  void addFreeEdge(const int u, const int v, const int delta) {
    isFreeEdge[u][v] += delta; 
    isFreeEdge[v][u] += delta;
  }

  void addCrossedEdge(const int u, const int v, const int delta) {
    isCrossedEdge[u][v] += delta; 
    isCrossedEdge[v][u] += delta;
  }

  /// check if takenCrossings can be reordered
  bool canSwap() const {
    CHECK(takenCrossings.size() >= 1);
    const int n = graph.n;
    const auto& edges = graph.edges;
    const auto& adj = graph.adj;

    // set current costs
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        cost[i][j] = (i == j ? 0 : -1);
      }
    }
    for (const auto& [e1, e2] : takenCrossings) {
      const int u1 = edges[e1].first;
      const int v1 = edges[e1].second;
      const int u2 = edges[e2].first;
      const int v2 = edges[e2].second;
      CHECK(u1 < v1);
      CHECK(u2 < v2);

      if (cost[u1][v1] == -1)
        cost[u1][v1] = cost[v1][u1] = 1;
      if (cost[u2][v2] == -1)
        cost[u2][v2] = cost[v2][u2] = 1;

      cost[u1][u2] = cost[u2][u1] = 0;
      cost[u2][v1] = cost[v1][u2] = 0;
      cost[v1][v2] = cost[v2][v1] = 0;
      cost[v2][u1] = cost[u1][v2] = 0;
    }

    // collect taken vertices and edges
    std::vector<int> takenVertices;
    takenVertices.reserve(4 * takenCrossings.size());
    for (const auto& [e1, e2] : takenCrossings) {
      const auto [u1, v1] = edges[e1];
      const auto [u2, v2] = edges[e2];

      // this is a K4 edge, which should be taken care of by earlier constraints
      if (cost[u1][v1] == 0 || cost[u2][v2] == 0)
        return false;
      CHECK(cost[u1][v1] == 1 && cost[u2][v2] == 1);

      takenVertices.push_back(u1);
      takenVertices.push_back(v1);
      takenVertices.push_back(u2);
      takenVertices.push_back(v2);
    }
    sort_unique(takenVertices);

    for (const int v : takenVertices) {
      for (const int u : adj[v]) {
        CHECK(u != v);
        if (cost[u][v] == -1) {
          cost[u][v] = 2;
          cost[v][u] = 2;
        }
      }
    }

    // check if swapping takenVertices[i1] and takenVertices[i2] is beneficial
    CHECK(numSwaps >= 2);
    for (size_t i1 = 0; i1 < takenVertices.size(); i1++) {
      for (size_t i2 = i1 + 1; i2 < takenVertices.size(); i2++) {
        if (adj[takenVertices[i1]].size() != adj[takenVertices[i2]].size())
          continue;
        if (canSwap2(takenVertices[i1], takenVertices[i2], takenVertices))
          return true;
      }
    }

    if (numSwaps >= 3) {
      // swap takenVertices[i1] -> takenVertices[i2] -> takenVertices[i3]
      for (size_t i1 = 0; i1 < takenVertices.size(); i1++) {
        for (size_t i2 = i1 + 1; i2 < takenVertices.size(); i2++) {
          for (size_t i3 = i1 + 1; i3 < takenVertices.size(); i3++) {
            if (i2 == i3)
              continue;
            if (adj[takenVertices[i1]].size() != adj[takenVertices[i2]].size())
              continue;
            if (adj[takenVertices[i1]].size() != adj[takenVertices[i3]].size())
              continue;
            if (adj[takenVertices[i2]].size() != adj[takenVertices[i3]].size())
              continue;
            if (canSwap3(takenVertices[i1], takenVertices[i2], takenVertices[i3], takenVertices)) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  /// check if S and T can be swapped
  bool canSwap2(const int S, const int T, const std::vector<int>& takenVertices) const {
    const bool Debug = false;

    // collect new edges
    newEdges.clear();
    for (const int v : takenVertices) {
      for (const int u : graph.adj[v]) {
        CHECK(u != v);
        const int s = v != S && v != T ? v : (v == S ? T : S);
        const int t = u != S && u != T ? u : (u == S ? T : S);

        const int u1 = std::min(s, t);
        const int v1 = std::max(s, t);
        // stop early
        if (cost[u1][v1] == -1)
          return false;

        if (cost[u1][v1] == 1) {
          newEdges.push_back({u1, v1});
        }
      }
    }

    if (newEdges.empty()) {
      LOG_IF(Debug, "yay-1, swapped %d and %d", S, T);
      return true;
    }

    sort_unique(newEdges);

    if (canSwapX()) {
      LOG_IF(Debug, "yay-2, swapped %d and %d", S, T);
      return true;
    }

    return false;
  }

  bool canSwap3(const int V1, const int V2, const int V3, const std::vector<int>& takenVertices) const {
    const bool Debug = false;

    // collect new edges
    newEdges.clear();
    for (const int v : takenVertices) {
      for (const int u : graph.adj[v]) {
        CHECK(u != v);
        int s = v;
        if (s == V1)
          s = V2;
        else if (s == V2)
          s = V3;
        else if (s == V3)
          s = V1;

        int t = u;
        if (t == V1)
          t = V2;
        else if (t == V2)
          t = V3;
        else if (t == V3)
          t = V1;

        const int u1 = std::min(s, t);
        const int v1 = std::max(s, t);
        // stop early
        if (cost[u1][v1] == -1)
          return false;

        if (cost[u1][v1] == 1) {
          newEdges.push_back({u1, v1});
        }
      }
    }
    sort_unique(newEdges);

    if (newEdges.empty()) {
      LOG_IF(Debug, "yay, swapped %d -> %d -> %d", V1, V2, V3);
      return true;
    }

    if (canSwapX()) {
      LOG_IF(Debug, "yay, swapped %d -> %d -> %d", V1, V2, V3);
      return true;
    }

    return false;
  }

  bool canSwapX() const {
    // assumes new edges are collected and sorted
    CHECK(!newEdges.empty());
    const bool Debug = false;

    // count crossings
    const int oldCrossings = (int)takenCrossings.size();
    int newCrossings = 0;
    std::vector<bool> processedNewEdge(newEdges.size(), false);
    for (size_t i = 0; i < newEdges.size(); i++) {
      if (newCrossings >= oldCrossings)
        break;
      if (processedNewEdge[i])
        continue;
      
      const auto [u, v] = newEdges[i];

      // case 1: cost-0 => drop
      if (cost[u][v] == 0) {
        LOG_IF(Debug, " droppped 1: (%d, %d)", u, v);
        processedNewEdge[i] = true;
        CHECK(false, "should be processed earlier");
        continue;
      }

      // case 2: cross edge
      if (cost[u][v] == 1) {
        std::pair<int, int> crossEdge = findSiblingCrossEdge(newEdges[i]);
        auto it = std::find(newEdges.begin(), newEdges.end(), crossEdge);
        if (it != newEdges.end()) {
          // case 2a: two from cross => +1
          const size_t j = std::distance(newEdges.begin(), it);
          CHECK(i < j);
          LOG_IF(Debug, " droppped 2a: (%d, %d)", newEdges[j].first, newEdges[j].second);
          processedNewEdge[j] = true;
          LOG_IF(Debug, " droppped 2a: (%d, %d)", u, v);
          processedNewEdge[i] = true;
          newCrossings++;
          continue;
        } else {
          // case 2b: one from cross => +0
          LOG_IF(Debug, " droppped 2b: (%d, %d)", u, v);
          processedNewEdge[i] = true;
          continue;
        }
      }

      // case 3: from other => drop
      if (cost[u][v] == 2) {
        LOG_IF(Debug, " droppped 3: (%d, %d)", u, v);
        processedNewEdge[i] = true;
        CHECK(false, "should be processed earlier");
        continue;
      }

      // case 4: something else => +INF
      newCrossings = oldCrossings + 1;
      CHECK(false);
    }

    // cannot swap
    if (newCrossings >= oldCrossings) {
      return false;
    }

    return true;
  }

  std::pair<int, int> findSiblingCrossEdge(const std::pair<int, int>& edge) const {
    for (const auto& [e1, e2] : takenCrossings) {
      if (graph.edges[e1] == edge) {
        return graph.edges[e2];
      }
      if (graph.edges[e2] == edge) {
        return graph.edges[e1];
      }
    }
    ERROR("unreachanble");
    return {};
  }

private:
  const InputGraph& graph;
  const int n;
  const int m;
  const int verbose;
  ForbiddenTuples& forbiddenTuples;
  std::vector<std::pair<int, int>> possibleCrossings;

  std::vector<std::pair<int, int>> takenCrossings;
  // std::vector<std::vector<std::pair<int, int>>> foundConstraints;
  size_t numConstraints2 = 0;
  size_t numConstraints3 = 0;
  std::vector<std::vector<int>> isFreeEdge;
  std::vector<std::vector<int>> isCrossedEdge;
  int numSwaps;

  mutable std::vector<std::vector<int>> cost;
  mutable std::vector<std::pair<int, int>> newEdges;
};

/// Disable pairs/triplse of crossings that can be eliminated by swapping some vertices
void encodeSwapConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  // TODO: need this??
  CHECK(!graph.isDirected());

  const int verbose = params.verbose;
  const auto sc = SplitNotNullInt(params.swapConstraints, "/");
  CHECK(sc.size() == 2, "incorrect format for swap-constraints");
  const int numPairs = sc[0];
  const int numSwaps = sc[1];
  CHECK(numPairs >= 1, "incorrect format for swap-constraints");
  CHECK(2 <= numSwaps && numSwaps <= 3, "numSwaps=%d is not supported", numSwaps);

  SwapTraversal finder(graph, params, model.getForbiddenTuples());
  if (!finder.init(numPairs, numSwaps)) {
    return;
  }

  // forbidden pairs (2-clauses)
  const size_t numClauses2 = finder.search2();
  LOG_IF(verbose, "  found %'9d 2-clauses", numClauses2);

  if (numPairs == 3) {
    // forbidden triples (3-clauses)
    const size_t numClauses3 = finder.search3();
    LOG_IF(verbose, "  found %'9d 3-clauses", numClauses3);
  }
}
