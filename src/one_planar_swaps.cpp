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

  std::vector<std::vector<std::pair<int, int>>> search(int maxNumCrossings, int maxNumSwaps) {
    isCrossed = std::vector<bool>(graph.edges.size(), false);
    takenCrossings.clear();
    foundConstraints.clear();
    cost = std::vector<std::vector<int>>(graph.n, std::vector<int>(graph.n, -1));
    numSwaps = maxNumSwaps;

    searchRec(0, maxNumCrossings);

    return foundConstraints;
  }

private:
  void searchRec(size_t curIdx, int remainingCrossings) {
    // Debug only
    // if (foundConstraints.size() >= 1000000)
    //   return;

    if (remainingCrossings >= 0) {
      if (!takenCrossings.empty() && canSwap()) {
        foundConstraints.push_back(takenCrossings);
        LOG_IF(verbose && foundConstraints.size() % 5000 == 0, "  found %d constraints so far...", foundConstraints.size());
        // print 
        LOG_IF(verbose >= 3, "found %d-th constraint:", foundConstraints.size());
        for (const auto& [e1, e2] : takenCrossings) {
          LOG_IF(verbose >= 3, "  (%d, %d) -- (%d, %d)", 
                 graph.edges[e1].first, graph.edges[e1].second, graph.edges[e2].first, graph.edges[e2].second);
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
      const int u1 = edges[e1].first;
      const int v1 = edges[e1].second;
      const int u2 = edges[e2].first;
      const int v2 = edges[e2].second;

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
              LOG("yay!");
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  // check if S and T can be swapped
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

    if (canSwapX(takenVertices)) {
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

    if (canSwapX(takenVertices)) {
      LOG_IF(Debug, "yay, swapped %d -> %d -> %d", V1, V2, V3);
      return true;
    }

    return false;
  }

  bool canSwapX(const std::vector<int>& takenVertices) const {
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
      
      const int u = newEdges[i].first;
      const int v = newEdges[i].second;

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
  const std::vector<std::pair<int, int>>& possibleCrossings;
  const int verbose;

  std::vector<bool> isCrossed;
  std::vector<std::pair<int, int>> takenCrossings;
  std::vector<std::vector<std::pair<int, int>>> foundConstraints;
  int numSwaps;

  mutable std::vector<std::vector<int>> cost;
  mutable std::vector<std::pair<int, int>> newEdges;
};

void encodeSwapConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  CHECK(!graph.isDirected());

  const auto sc = SplitNotNullInt(params.swapConstraints, "/");
  CHECK(sc.size() == 2, "incorrect format for swap-constraints");
  const int numPairs = sc[0];
  const int numSwaps = sc[1];
  CHECK(numPairs >= 1 && numSwaps >= 2, "incorrect format for swap-constraints");
  CHECK(numSwaps <= 3, "numSwaps=%d is not implemeted yet", numSwaps);
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

  SwapFinder finder(graph, possibleCrossings, params.verbose);
  auto constraints = finder.search(numPairs, numSwaps);

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
    LOG_IF(params.verbose && numConstraints[i] > 0, "  found %'d %d-swap constraints", numConstraints[i], i);
  }
}

/// Disable pairs of crossings that can be eliminated by swapping pairs of vertices
void encodeSwapConstraintsNaive(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const auto& adj = graph.adj;

  int numConstraints = 0;
  // try to eliminate crossings (x, b)--(y, a) and (x, v)--(y, u) by exchanging x and y
  for (int x = 0; x < n; x++) {
    for (int y = x + 1; y < n; y++) {
      for (size_t ib = 0; ib < adj[x].size(); ib++) {
        for (size_t iv = 0; iv < adj[x].size(); iv++) {
          if (ib == iv) continue;
          const int b = adj[x][ib];
          const int v = adj[x][iv];
          if (b == y || v == y) continue;
          CHECK(b != v);

          const int be = graph.findDivIndex(b, x);
          const int ve = graph.findDivIndex(v, x);
          for (size_t ia = 0; ia < adj[y].size(); ia++) {
            for (size_t iu = 0; iu < adj[y].size(); iu++) {
              if (ia == iu) continue;
              const int a = adj[y][ia];
              const int u = adj[y][iu];
              if (a == x || u == x) continue;
              CHECK(a != u);

              const int ae = graph.findDivIndex(a, y);
              const int ue = graph.findDivIndex(u, y);

              if (b == a || b == u || v == u || v == a)
                continue;

              CHECK(be != ae && ve != ue);
              if (!canBeMerged(be, ae, n, edges))
                continue;
              if (!canBeMerged(ve, ue, n, edges))
                continue;
              CHECK(b != a && b != u && v != u && v != a);

              // vertices that can be reached from x w/o violating 1-planarity
              std::vector<int> reachX = adj[x];
              CHECK(contains(reachX, b) && contains(reachX, b));
              reachX.push_back(x);
              reachX.push_back(y);
              reachX.push_back(a);
              reachX.push_back(u);
              // vertices that can be reached from y w/o violating 1-planarity
              std::vector<int> reachY = adj[y];
              CHECK(contains(reachY, a) && contains(reachY, u));
              reachY.push_back(y);
              reachY.push_back(x);
              reachY.push_back(b);
              reachY.push_back(v);

              if (!equal_unsorted(reachX, reachY))
                continue;

              if (std::minmax(ae, be) < std::minmax(ve, ue)) {
                // TODO: this isn't always correct, since we do not necessarily reduce the number of crossings
                LOG_IF(verbose >= 3, "swap-2 crossing (x = %d; y = %d)", x, y);
                LOG_IF(verbose >= 3, "  (%d, %d) -- (%d, %d)", min(x, b), max(x, b), min(y, a), max(y, a));
                LOG_IF(verbose >= 3, "  (%d, %d) -- (%d, %d)", min(x, v), max(x, v), min(y, u), max(y, u));

                model.addClause({
                    model.getCross2Var(ae, be, false), 
                    model.getCross2Var(ve, ue, false)
                });
                numConstraints++;
              }
            }
          }
        }
      }
    }
  }

  LOG_IF(verbose, "added %3d 2-swap constraints", numConstraints);

  int num3Constraints = 0;
  // forbid crossings between every pair of edges for two degree-3 vertices
  for (int u = 0; u < n; u++) {
    if (graph.degree(u) != 3) 
      continue;
    for (int v = u + 1; v < n; v++) {
      if (graph.degree(v) != 3) 
        continue;
      if (graph.hasEdge(u, v))
        continue;

      // vertex u is adjacent to (a, b, c) via edges (ua, ub, uc)
      const int a = adj[u][0];
      const int ua = graph.findDivIndex(u, a);
      const int b = adj[u][1];
      const int ub = graph.findDivIndex(u, b);
      const int c = adj[u][2];
      const int uc = graph.findDivIndex(u, c);
      vector<int> perm = identity(3);

      do {
        // vertex v is adjacent to (x, y, z) via edges (vx, vy, vz)
        const int x = adj[v][perm[0]];
        const int vx = graph.findDivIndex(v, x);
        const int y = adj[v][perm[1]];
        const int vy = graph.findDivIndex(v, y);
        const int z = adj[v][perm[2]];
        const int vz = graph.findDivIndex(v, z);

        if (!canBeMerged(ua, vx, n, edges))
          continue;
        if (!canBeMerged(ub, vy, n, edges))
          continue;
        if (!canBeMerged(uc, vz, n, edges))
          continue;

        LOG_IF(verbose >= 3, "swap-3 crossing: (%d--%d); (%d--%d); (%d--%d)", ua, vx, ub, vy, uc, vz);
        model.addClause({
            model.getCross2Var(ua, vx, false), 
            model.getCross2Var(ub, vy, false),
            model.getCross2Var(uc, vz, false)
        });
        num3Constraints++;
      } while (std::next_permutation(perm.begin(), perm.end()));
    }
  }

  LOG_IF(verbose, "added %3d 3-swap constraints", num3Constraints);
}
