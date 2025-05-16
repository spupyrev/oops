#include "one_planar.h"
#include "logging.h"

using namespace std;

/// TODO
struct SwapFinder {
  SwapFinder(const InputGraph& graph, const std::vector<std::pair<int, int>>& possibleCrossings, const int verbose) 
    : graph(graph),
      possibleCrossings(possibleCrossings),
      verbose(verbose)
    {}

  std::vector<std::vector<std::pair<int, int>>> search(int numCrossings, int numSwaps) {
    isCrossed = std::vector<bool>(graph.edges.size(), false);
    takenCrossings.clear();
    foundConstraints.clear();
    cost = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, -1));

    searchRec(0, numCrossings);

    return foundConstraints;
  }

private:
  void searchRec(size_t curIdx, int remainingCrossings) {
    if (remainingCrossings == 0) {
      if (canSwap()) {
        foundConstraints.push_back(takenCrossings);
        // print 
        LOG_IF(verbose >= 3, "found %d-th contraint:", foundConstraints.size());
        for (const auto& [e1, e2] : takenCrossings) {
          LOG_IF(verbose >= 3, "  (%d, %d) -- (%d, %d)", 
                 graph.edges[e1].first, graph.edges[e1].second, graph.edges[e2].first, graph.edges[e2].second);
        }
      }
      return;
    }
    CHECK(possibleCrossings.size() >= curIdx);
    const int leftCrossings = int(possibleCrossings.size() - curIdx);
    if (remainingCrossings > leftCrossings) {
      return;
    }
    if (curIdx == possibleCrossings.size()) {
      return;
    }

    // TODO: somehow stop if a subset is forbidden

    // try to get current
    const int e1 = possibleCrossings[curIdx].first;
    const int e2 = possibleCrossings[curIdx].second;
    if (!isCrossed[e1] && !isCrossed[e2]) {
      isCrossed[e1] = true;
      isCrossed[e2] = true;
      takenCrossings.push_back(possibleCrossings[curIdx]);
      searchRec(curIdx + 1, remainingCrossings - 1);
      isCrossed[e1] = false;
      isCrossed[e2] = false;
      takenCrossings.pop_back();
    }

    // skip current
    searchRec(curIdx + 1, remainingCrossings);
  }

  // check if takenCrossings can be reordered
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
    std::vector<int> takenVertices;
    std::vector<int> takenCrossEdges;
    for (const auto& [e1, e2] : takenCrossings) {
      takenCrossEdges.push_back(e1);
      takenCrossEdges.push_back(e2);
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

      cost[u1][u2] = 0;
      cost[u2][u1] = 0;

      cost[u2][v1] = 0;
      cost[v1][u2] = 0;

      cost[v1][v2] = 0;
      cost[v2][v1] = 0;

      cost[v2][u1] = 0;
      cost[u1][v2] = 0;

      takenVertices.push_back(u1);
      takenVertices.push_back(v1);
      takenVertices.push_back(u2);
      takenVertices.push_back(v2);
    }
    sort_unique(takenVertices);

    // verify that input is sane
    sort_unique(takenCrossEdges);
    // no overlapping crossing pairs
    CHECK(takenCrossEdges.size() == 2 * takenCrossings.size());
    // no trivially re-routable edges
    for (const int e : takenCrossEdges) {
      const int u = edges[e].first;
      const int v = edges[e].second;
      // hmm?
      if (cost[u][v] == 0)
        return false;
      CHECK(cost[u][v] == 1);
    }
    //return false;

    // collect current edges
    std::vector<std::pair<int, int>> curCrossEdges;
    for (const int e : takenCrossEdges) {
      CHECK(edges[e].first < edges[e].second);
      const int u = edges[e].first;
      const int v = edges[e].second;
      curCrossEdges.push_back({u, v});
    }
    std::vector<std::pair<int, int>> curOtherEdges;
    for (const int v : takenVertices) {
      for (const int u : adj[v]) {
        CHECK(u != v);
        std::pair<int, int> e = {std::min(u, v), std::max(u, v)};
        if (!contains(curCrossEdges, e) && !contains(curOtherEdges, e))
          curOtherEdges.push_back(e);
      }
    }
    auto findSiblingCrossEdge = [&](const std::vector<std::pair<int, int>>& newEdges, const std::pair<int, int>& edge) -> int {
      CHECK(edge.first < edge.second);
      for (const auto& [e1, e2] : takenCrossings) {
        if (edges[e1] == edge) {
          auto it = std::find(newEdges.begin(), newEdges.end(), edges[e2]);
          if (it == newEdges.end())
            return -1;
          return int(it - newEdges.begin());
        }
        if (edges[e2] == edge) {
          auto it = std::find(newEdges.begin(), newEdges.end(), edges[e1]);
          if (it == newEdges.end())
            return -1;
          return int(it - newEdges.begin());
        }
      }
      return -1;
    };

    // check if swapping takenVertices[i1] and takenVertices[i2] is beneficial 
    for (size_t i1 = 0; i1 < takenVertices.size(); i1++) {
      for (size_t i2 = i1 + 1; i2 < takenVertices.size(); i2++) {
        // collect new edges
        std::vector<std::pair<int, int>> newEdges;
        for (const int v : takenVertices) {
          for (const int u : adj[v]) {
            CHECK(u != v);
            int s, t;
            if (v != takenVertices[i1] && v != takenVertices[i2]) 
              s = v;
            else if (v == takenVertices[i1])
              s = takenVertices[i2];
            else
              s = takenVertices[i1];

            if (u != takenVertices[i1] && u != takenVertices[i2]) 
              t = u;
            else if (u == takenVertices[i1])
              t = takenVertices[i2];
            else
              t = takenVertices[i1];

            newEdges.push_back({std::min(s, t), std::max(s, t)});
          }
        }
        sort_unique(newEdges);
        CHECK(newEdges.size() == curCrossEdges.size() + curOtherEdges.size());

        // count crossings
        int oldCrossings = (int)takenCrossings.size();
        int newCrossings = 0;
        for (int i = 0; i < (int)newEdges.size(); i++) {
          // case 1: cost-0 => drop
          if (cost[newEdges[i].first][newEdges[i].second] == 0) {
            remove_value(newEdges, newEdges[i]);
            i--;
            continue;
          }

          // case 2: cross edge
          if (contains(curCrossEdges, newEdges[i])) {
            const int j = findSiblingCrossEdge(newEdges, newEdges[i]);
            CHECK(i != j);
            if (j != -1) {
              // case 2a: two from cross => +1
              CHECK(i < j);
              remove_value(newEdges, newEdges[i]);
              remove_value(newEdges, newEdges[j]);
              i--;
              newCrossings++;
              continue;
            } else {
              // case 2b: one from cross => +0
              remove_value(newEdges, newEdges[i]);
              i--;
              continue;
            }
          }

          // case 3: from other => drop
          if (contains(curOtherEdges, newEdges[i])) {
            remove_value(newEdges, newEdges[i]);
            i--;
            continue;
          }

          // case 4: something else => +INF
          newCrossings = oldCrossings + 1;
          break;
        }

        if (newCrossings >= oldCrossings) {
          // cannot reroute
          continue;
        }

        // yay!
        CHECK(newEdges.empty());
        return true;
      }
    }

    return false;
  }

private:
  const InputGraph& graph;
  const std::vector<std::pair<int, int>>& possibleCrossings;
  const int verbose;

  std::vector<bool> isCrossed;
  std::vector<std::pair<int, int>> takenCrossings;
  std::vector<std::vector<std::pair<int, int>>> foundConstraints;
  mutable std::vector<std::vector<int>> cost;
};

void encodeSwapConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  CHECK(!graph.isDirected());

  const auto sc = SplitNotNullInt(params.swapConstraints, "/");
  CHECK(sc.size() == 2, "incorrect format for swap-constraints");
  const int numPairs = sc[0];
  const int numSwaps = sc[1];
  CHECK(numPairs >= 1 && numSwaps >= 2, "incorrect format for swap-constraints");
  LOG_IF(params.verbose, "encoding swap constraints for %d pairs and %d swaps", numPairs, numSwaps);

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

  for (int c = 1; c <= numPairs; c++) {
    SwapFinder finder(graph, possibleCrossings, params.verbose);
    CHECK(numSwaps == 2, "numSwaps=%d is not implemeted yet", numSwaps);
    auto constraints = finder.search(c, numSwaps);
    if (c == 1) {
      CHECK(constraints.empty());
      continue;
    }

    LOG_IF(params.verbose, "  found %d %d-swap constraints", constraints.size(), numSwaps);
    for (const auto& constraint: constraints) {
      CHECK((int)constraint.size() == c);
      MClause clause;
      for (const auto& [e1, e2] : constraint) {
        clause.addVar(model.getCross2Var(e1 + n, e2 + n, false));
      }
      model.addClause(clause);
    }
  }
}