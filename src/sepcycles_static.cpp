#include "graph_algorithms.h"
#include "logging.h"
#include "one_planar.h"

struct SepCyclesStatic {
  SepCyclesStatic(SATModel& model, const InputGraph& graph, const Params& params, ForbiddenTuples& tuples)
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
    crossingTwin = std::vector<std::vector<EdgeTy>>(n, std::vector<EdgeTy>(n, EdgeTy(-1, -1)));
    isVertexOnCycle = std::vector<char>(n, 0);
    isVertexBlocked = std::vector<int>(n, 0);
    cycleBuffer.reserve(10);

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
    if (possibleCrossings.size() > MAX_CROSSINGS_FOR_2GT) {
      LOG_IF(verbose, "  skipped building sep-cycles 2-clauses due to too many possible_crossings");
      return 0;
    }

    size_t numOverfullClauses2 = 0;
    size_t numTightClauses2 = 0;
    size_t numTightClauses3 = 0;

    // the main crossing
    for (size_t cross0 = 0; cross0 < possibleCrossings.size(); cross0++) {
      const auto [e1_cross0, e2_cross0] = possibleCrossings[cross0];
      const auto [u, v] = graph.edges[e1_cross0];
      const auto [x, y] = graph.edges[e2_cross0];
      CHECK(all_unique({x, y, u, v}));

      // circular order: u--x--v--y
      markCrossing(u, v, x, y);

      std::vector<std::pair<int, int>> takenCrossings0 = {{e1_cross0, e2_cross0}};
      numTightClauses2 += findTightSeparators(x, y, u, v, takenCrossings0);
      numTightClauses2 += findTightSeparators(u, v, x, y, takenCrossings0);

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
        if (findOverfullSeparator(x, y, u, v, takenCrossings) || findOverfullSeparator(u, v, x, y, takenCrossings)) {
          // save the 2-clause
          numOverfullClauses2++;
          forbiddenTuples.insert(crossPair);
        }

        unmarkCrossing(u2, v2, x2, y2);
      }

      unmarkCrossing(u, v, x, y);
    }

    // tight 3-clauses
    if (possibleCrossings.size() < MAX_CROSSINGS_FOR_3EQ) {
      for (size_t cross0 = 0; cross0 < possibleCrossings.size(); cross0++) {
        const auto [e1_cross0, e2_cross0] = possibleCrossings[cross0];
        const auto [u, v] = graph.edges[e1_cross0];
        const auto [x, y] = graph.edges[e2_cross0];
        CHECK(all_unique({x, y, u, v}));

        markCrossing(u, v, x, y);

        for (size_t cross1 = 0; cross1 < possibleCrossings.size(); cross1++) {
          if (cross0 == cross1)
            continue;
          const auto [e1_cross1, e2_cross1] = possibleCrossings[cross1];
          const auto [u2, v2] = graph.edges[e1_cross1];
          const auto [x2, y2] = graph.edges[e2_cross1];

          if (!canMarkCrossing(u2, v2, x2, y2))
            continue;
          CHECK(all_unique({e1_cross0, e2_cross0, e1_cross1, e2_cross1}));

          const CrossingPair crossPair(m, e1_cross0, e2_cross0, e1_cross1, e2_cross1);
          if (forbiddenTuples.contains(crossPair))
            continue;

          markCrossing(u2, v2, x2, y2);

          std::vector<std::pair<int, int>> takenCrossings = {{e1_cross0, e2_cross0}, {e1_cross1, e2_cross1}};
          numTightClauses3 += findTightSeparators(x, y, u, v, takenCrossings);
          numTightClauses3 += findTightSeparators(u, v, x, y, takenCrossings);

          unmarkCrossing(u2, v2, x2, y2);
        }

        unmarkCrossing(u, v, x, y);
      }
    }

    LOG_IF(verbose, "  added %'9d tight-separator 2-clauses", numTightClauses2);
    LOG_IF(verbose, "  added %'9d tight-separator 3-clauses", numTightClauses3);

    return numOverfullClauses2;
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
          if (findOverfullSeparator(x, y, u, v, takenCrossings) || findOverfullSeparator(u, v, x, y, takenCrossings)) {
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

  void setCrossingTwin(int u, int v, const EdgeTy& twin) {
    crossingTwin[u][v] = twin;
    crossingTwin[v][u] = twin;
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
    setCrossingTwin(u, v, make_edge(x, y));
    setCrossingTwin(x, y, make_edge(u, v));

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
    CHECK(crossingTwin[u][v] == make_edge(x, y));
    CHECK(crossingTwin[x][y] == make_edge(u, v));
    // uncross
    setEdge(isCrossEdge, u, v, false);
    setEdge(isCrossEdge, x, y, false);
    setCrossingTwin(u, v, EdgeTy(-1, -1));
    setCrossingTwin(x, y, EdgeTy(-1, -1));

    // kite: circular order is u--x--v--y
    unmarkKiteEdge(y, u);
    unmarkKiteEdge(v, y);
    unmarkKiteEdge(x, v);
    unmarkKiteEdge(u, x);
  }

  bool isFreeEdge(int u, int v) const {
    return !isCrossEdge[u][v] && !isKiteEdge[u][v];
  };

  // Upper bound on paths that can enter/leave an edge without using the edge itself.
  int edgeDegreeBound(int u, int v) const {
    return std::min(graph.degree(u) - 1, graph.degree(v) - 1);
  }

  int edgeDegreeBound(const EdgeTy& edge) const {
    return edgeDegreeBound(edge.first, edge.second);
  }

  struct SeparatorDfsResult {
    bool foundOverfull = false;
    int numTightClauses = 0;
  };

  // Constant data for one DFS from x to y. The edge (x,y) is already a marked
  // crossing edge; the DFS path from x to y, together with (y,x), forms the
  // candidate separator.
  struct SeparatorDfsContext {
    bool enableTestForFewCrossings;
    bool searchOverfull;
    bool searchTight;
    const std::vector<std::pair<int, int>>& takenCrossings;
    // Optimistic shortest-path distances to y in the currently marked graph.
    // Dynamic DFS blocks are ignored, so these distances can only make paths
    // look easier than they really are; that makes distance pruning safe.
    const std::vector<int>& distToY;
    int x;
    int y;
    int u;
    int v;
  };

  // Mutable data carried by recursive DFS calls.
  struct SeparatorDfsState {
    int remainingRecDepth;
    int numFreeEdgesOnCycle;
    // Upper bounds on how many disjoint paths may still cross the separator.
    // optimisticPathUB may still use not-yet-seen crossing edges; representedPathUB
    // counts only crossing edges already represented by the current cycle.
    int optimisticPathUB;
    int representedPathUB;
  };

  /// Find all separating cycles: paths from x to y avoiding vertices {x, y, u, v}
  /// A path may contain as many kite/cross edges as possible but only one free edge
  /// If a path contains one crossed edge, then disjoint paths cannot use the second one
  bool findOverfullSeparator(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings) {
    const int numTakenCrossings = (int)takenCrossings.size();
    CHECK(numTakenCrossings == 2 || numTakenCrossings == 3);
    const int maxCycleLen = (numTakenCrossings == 2 ? MAX_CYCLE_LEN_FOR_SEPCYCLES : 4);
    return findSeparators(x, y, u, v, takenCrossings, maxCycleLen, true, false).foundOverfull;
  }

  /// Find tight separators: cycles where the number of free cycle edges equals
  /// the number of disjoint paths crossing the cycle. Each free edge then gives
  /// one clause saying that this edge cannot stay free with the chosen crossings.
  int findTightSeparators(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings) {
    CHECK(takenCrossings.size() == 1 || takenCrossings.size() == 2);
    const int degreeBound = edgeDegreeBound(u, v);
    const int maxCycleLen = takenCrossings.size() == 1 ? std::min(MAX_CYCLE_LEN_FOR_SEPCYCLES, degreeBound + 1) : 4;
    return findSeparators(x, y, u, v, takenCrossings, maxCycleLen, false, true).numTightClauses;
  }

  SeparatorDfsResult findSeparators(int x, int y, int u, int v,
                                    const std::vector<std::pair<int, int>>& takenCrossings,
                                    int maxCycleLen,
                                    bool searchOverfull,
                                    bool searchTight) {
    CHECK(searchOverfull != searchTight);
    CHECK(isCrossEdge[x][y] && isKiteEdge[x][y] == 0);

    // The cycle already contains the represented crossing edge (x,y).
    // Vertices u and v are blocked because they are the endpoints of its twin edge.
    const int numTakenCrossings = (int)takenCrossings.size();
    const bool enableTestForFewCrossings = searchOverfull && (numTakenCrossings == 2 || params.unsatLevel >= 2);
    int optimisticPathUB = edgeDegreeBound(u, v);
    const int representedPathUB = optimisticPathUB;
    if (enableTestForFewCrossings) {
      optimisticPathUB = 0;
      for (const auto& [e1, e2] : takenCrossings) {
        optimisticPathUB += std::max(edgeDegreeBound(graph.edges[e1]), edgeDegreeBound(graph.edges[e2]));
      }
      const int primaryMaxPathUB = std::max(edgeDegreeBound(x, y), representedPathUB);
      optimisticPathUB += representedPathUB - primaryMaxPathUB;
    }
    CHECK(representedPathUB <= optimisticPathUB);

    isVertexOnCycle[x] = 1;
    isVertexBlocked[x]++;
    isVertexBlocked[u]++;
    isVertexBlocked[v]++;
    cycleBuffer.clear();
    cycleBuffer.push_back(x);
    const int maxDepth = maxCycleLen - 1;

    // Precompute how many DFS steps are needed to reach y from each vertex.
    // The graph here is the current adjList: original edges plus pseudo-kite
    // edges for the crossings already marked by the caller. We cap distances at
    // maxDepth + 1 because larger values are all equivalent for pruning.
    distToYBuffer.assign(n, maxDepth + 1);
    distToYBuffer[y] = 0;
    distQueue.clear();
    distQueue.push_back(y);
    for (size_t i = 0; i < distQueue.size(); i++) {
      const int cur = distQueue[i];
      if (distToYBuffer[cur] == maxDepth)
        continue;
      for (const int next : adjList[cur]) {
        if (isVertexBlocked[next] && next != y)
          continue;
        if (distToYBuffer[next] <= distToYBuffer[cur] + 1)
          continue;
        distToYBuffer[next] = distToYBuffer[cur] + 1;
        distQueue.push_back(next);
      }
    }
    const SeparatorDfsContext context = {
        enableTestForFewCrossings, searchOverfull, searchTight, takenCrossings, distToYBuffer, x, y, u, v};
    const SeparatorDfsState initialState = {
        maxDepth, 0, optimisticPathUB, representedPathUB};
    SeparatorDfsResult result;
    findSeparatorsRec(x, context, initialState, result);

    cycleBuffer.clear();
    isVertexOnCycle[x] = 0;
    isVertexBlocked[x]--;
    isVertexBlocked[u]--;
    isVertexBlocked[v]--;

    return result;
  }

  // Reject a partial path once the number of free edges is already too large.
  // For overfull separators, "too large" means there can no longer be enough
  // disjoint paths crossing the separator. For tight separators, we only need
  // cycles with at most the path upper bound many free edges.
  bool shouldPruneSeparatorDfs(const SeparatorDfsContext& context, const SeparatorDfsState& state) const {
    if (context.searchTight)
      return state.numFreeEdgesOnCycle > state.optimisticPathUB;

    CHECK(context.searchOverfull);
    if (!context.enableTestForFewCrossings)
      return state.numFreeEdgesOnCycle >= state.optimisticPathUB;

    CHECK(state.representedPathUB <= state.optimisticPathUB);
    const int pathUB = state.remainingRecDepth == 0 ? state.representedPathUB : state.optimisticPathUB;
    return state.numFreeEdgesOnCycle >= pathUB;
  }

  // An overfull separator is a cycle whose free edges are fewer than the number
  // of edge-disjoint paths forced to cross it.
  bool processOverfullSeparatorCycle(const SeparatorDfsContext& context, const SeparatorDfsState& state) {
    CHECK(cycleBuffer.front() == context.x && cycleBuffer.back() == context.y && cycleBuffer.size() >= 3);

    if (!context.enableTestForFewCrossings) {
      CHECK(state.numFreeEdgesOnCycle < state.optimisticPathUB);
      return hasManyDisjointPathsForOneCrossing(
          context.u, context.v, context.takenCrossings, cycleBuffer, state.numFreeEdgesOnCycle + 1);
    }

    CHECK(state.numFreeEdgesOnCycle < state.representedPathUB);
    return hasManyDisjointPathsForFewCrossings(context.takenCrossings, cycleBuffer, state.numFreeEdgesOnCycle);
  }

  // A tight separator has exactly as many free cycle edges as crossing paths.
  // Each free edge of such a cycle produces one clause forbidding that edge to stay free.
  int processTightSeparatorCycle(const SeparatorDfsContext& context, const SeparatorDfsState& state) {
    CHECK(cycleBuffer.front() == context.x && cycleBuffer.back() == context.y && cycleBuffer.size() >= 3);
    if (state.numFreeEdgesOnCycle == 0 || state.numFreeEdgesOnCycle > state.optimisticPathUB)
      return 0;
    if (!hasManyDisjointPathsForOneCrossing(
            context.u, context.v, context.takenCrossings, cycleBuffer, state.numFreeEdgesOnCycle))
      return 0;

    int numAddedClauses = 0;
    for (size_t i = 0; i + 1 < cycleBuffer.size(); i++) {
      const int v1 = cycleBuffer[i];
      const int v2 = cycleBuffer[i + 1];
      if (!isEdge[v1][v2] || !isFreeEdge(v1, v2))
        continue;
      numAddedClauses++;
      addTightSeparatorClause(context.takenCrossings, v1, v2);
    }
    return numAddedClauses;
  }

  // Extend the current x-to-y path. Edges in adjList are original graph edges plus
  // pseudo-kite edges introduced by currently marked crossings.
  bool findSeparatorsRec(int currentVertex,
                         const SeparatorDfsContext& context,
                         const SeparatorDfsState& state,
                         SeparatorDfsResult& result) {
    // Reaching y closes the separator; cycleBuffer is the path and the closing
    // edge (y,x) is the represented crossing edge.
    if (currentVertex == context.y) {
      CHECK(cycleBuffer.size() >= 3);
      CHECK(cycleBuffer.front() == context.x && cycleBuffer.back() == context.y);
      if (context.searchTight) {
        result.numTightClauses += processTightSeparatorCycle(context, state);
        return false;
      }

      CHECK(context.searchOverfull);
      result.foundOverfull = processOverfullSeparatorCycle(context, state);
      return result.foundOverfull;
    }

    if (state.remainingRecDepth == 0)
      return false;

    for (const int nextVertex : adjList[currentVertex]) {
      // Used vertices are on the current cycle; blocked vertices are endpoints
      // of twin crossing edges that the cycle is not allowed to touch.
      if (isVertexBlocked[nextVertex])
        continue;
      if (nextVertex == context.y && cycleBuffer.size() < 2)
        continue;

      // After taking currentVertex -> nextVertex, there must still be enough
      // depth left to reach y. distToY ignores dynamic DFS blocks, so failing
      // this test means no valid constrained continuation can exist.
      const int nextRemainingRecDepth = nextVertex == context.y ? 0 : state.remainingRecDepth - 1;
      if (context.distToY[nextVertex] > nextRemainingRecDepth)
        continue;

      SeparatorDfsState nextState = state;
      nextState.remainingRecDepth = nextRemainingRecDepth;
      nextState.numFreeEdgesOnCycle += isFreeEdge(currentVertex, nextVertex);

      int blockedU = -1;
      int blockedV = -1;
      if (isCrossEdge[currentVertex][nextVertex]) {
        const auto [u, v] = crossingTwin[currentVertex][nextVertex];
        CHECK(u != -1 && v != -1);
        blockedU = u;
        blockedV = v;
        // If the cycle uses one edge of a crossing, it may not also visit an
        // endpoint of the twin edge.
        if (isVertexOnCycle[blockedU] || isVertexOnCycle[blockedV])
          continue;
        if (context.enableTestForFewCrossings) {
          // The crossing is now represented by a specific edge on the cycle, so
          // replace the previous optimistic bound by the bound for the twin edge.
          const int edgePathUB = edgeDegreeBound(currentVertex, nextVertex);
          const int twinPathUB = edgeDegreeBound(blockedU, blockedV);
          nextState.optimisticPathUB += twinPathUB - std::max(edgePathUB, twinPathUB);
          nextState.representedPathUB += twinPathUB;
          CHECK(nextState.representedPathUB <= nextState.optimisticPathUB);
        }
      }

      if (shouldPruneSeparatorDfs(context, nextState))
        continue;

      // Commit the step only after pruning, so rejected states do not mutate
      // the global DFS buffers.
      if (blockedU != -1) {
        isVertexBlocked[blockedU]++;
        isVertexBlocked[blockedV]++;
      }

      isVertexOnCycle[nextVertex] = 1;
      isVertexBlocked[nextVertex]++;
      cycleBuffer.push_back(nextVertex);
      const bool stop = findSeparatorsRec(nextVertex, context, nextState, result);
      cycleBuffer.pop_back();
      isVertexBlocked[nextVertex]--;
      isVertexOnCycle[nextVertex] = 0;
      if (blockedU != -1) {
        isVertexBlocked[blockedU]--;
        isVertexBlocked[blockedV]--;
      }
      if (stop)
        return true;
    }
    return false;
  }

  /// Check if there are more than `numFree` edge-disjoint paths from u to v disjoint from the cycle.
  /// It is assumed that the cycle contains at most one edge from each crossing.
  bool hasManyDisjointPathsForFewCrossings(const std::vector<std::pair<int, int>> &takenCrossings,
                                           const std::vector<int>& cycle,
                                           const int numFree) const {
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
      // Need at most one edge from each crossing on the cycle
      CHECK(!edge1OnCycle || !edge2OnCycle);
      if (!edge1OnCycle && !edge2OnCycle)
        continue;
      crossCycleEdges.push_back(edge1OnCycle ? edge2 : edge1);
    }
    CHECK(!crossCycleEdges.empty());

    // The twin edge must be completely outside the separator cycle.
    for (size_t i = 0; i < crossCycleEdges.size(); i++) {
      const auto [x, y] = crossCycleEdges[i];
      if (contains(cycle, x) || contains(cycle, y))
        return false;
    }

    int degreeBound = 0;
    for (const auto& [x, y] : crossCycleEdges)
      degreeBound += edgeDegreeBound(x, y);
    if (numFree >= degreeBound)
      return false;

    // Orientation masks are symmetric under flipping all pairs at once
    // (swap sources <-> targets globally), so fix pair 0 as (a -> b).
    // Only pairs 1..k-1 need flipping choices.
    const int numPairs = static_cast<int>(crossCycleEdges.size());
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

  bool hasManyDisjointPathsForOneCrossing(int u, int v,
                                          const std::vector<std::pair<int, int>>& takenCrossings,
                                          const std::vector<int>& cycle,
                                          int numPathsLowerBound) const {
    const auto crossCycleEdges = getCrossCycleEdges(cycle, takenCrossings);
    const int numPaths = countEdgeDisjointPaths(
        {u}, {v}, graph.adj, cycle, crossCycleEdges, numPathsLowerBound);
    return numPaths >= numPathsLowerBound;
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

  void addTightSeparatorClause(const std::vector<std::pair<int, int>>& takenCrossings, int u, int v) {
    CHECK(takenCrossings.size() == 1 || takenCrossings.size() == 2);
    const int divUV = graph.findDivIndex(u, v);
    MClause clause;
    for (const auto& [e1, e2] : takenCrossings) {
      clause.addVar(model.getCross2Var(e1 + n, e2 + n, false));
    }
    clause.addVar(model.getCross1Var(divUV, true));
    model.addClause(std::move(clause));
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
  std::vector<std::vector<EdgeTy>> crossingTwin;
  std::vector<char> isVertexOnCycle;
  std::vector<int> isVertexBlocked;
  std::vector<int> cycleBuffer;
  std::vector<int> distToYBuffer;
  std::vector<int> distQueue;

};

/// Add constraints based on separating cycles
void encodeSepCyclesConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  const int verbose = params.verbose;
  const int numPairs = to_int(params.sepCycleConstraints);
  CHECK(2 <= numPairs && numPairs <= 3, "incorrect value of sepCycleConstraints");

  SepCyclesStatic sct(model, graph, params, model.getForbiddenTuples());
  if (sct.init(numPairs)) {
    // forbidden pairs (2-clauses)
    const size_t numOverfullClauses2 = sct.build2Clauses();
    LOG_IF(verbose, "  found %'9d overfull-separator 2-clauses", numOverfullClauses2);

    if (numPairs == 3 || params.unsatLevel >= 2) {
      // forbidden triples (3-clauses)
      const size_t numClauses3 = sct.build3Clauses();
      LOG_IF(verbose, "  found %'9d 3-clauses", numClauses3);
    }
  }
}
