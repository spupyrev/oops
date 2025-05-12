#include "graph_algorithms.h"
#include "one_planar.h"

/// Check if a given graph with specified pairs of crossed edges is planar
bool isPlanarWithCrossings(const InputGraph& graph, const std::vector<std::pair<int, int>>& crossings) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  std::vector<bool> isCrossed(m, false);
  for (const auto& [e1, e2] : crossings) {
    CHECK(0 <= e1 && e1 < m);
    CHECK(0 <= e2 && e2 < m);
    CHECK(e1 != e2);
    CHECK(!isCrossed[e1] && !isCrossed[e2]);
    isCrossed[e1] = true;
    isCrossed[e2] = true;
  }

  std::vector<EdgeTy> planarEdges;
  planarEdges.reserve(graph.edges.size() + 4 * crossings.size());
  for (size_t i = 0; i < edges.size(); i++) {
    if (!isCrossed[i]) {
      planarEdges.push_back(edges[i]);
    }
  }

  int nPlanar = n;
  for (const auto& [e1, e2] : crossings) {
    planarEdges.push_back({edges[e1].first, nPlanar});
    planarEdges.push_back({edges[e1].second, nPlanar});
    planarEdges.push_back({edges[e2].first, nPlanar});
    planarEdges.push_back({edges[e2].second, nPlanar});
    nPlanar++;
  }

  return isPlanar(nPlanar, planarEdges, 1);
}

/// Adjust the result by removing unnecessary crossings
void minimizeCrossings(const InputGraph& graph, Result& result, int verbose) {
  const auto& edges = graph.edges;
  const int orgNumCrossings = (int)result.crossings.size();
  std::vector<std::pair<int, int>> reqCrossings = result.crossings;
  for (size_t i = 0; i < result.crossings.size(); i++) {
    auto cr = result.crossings[i];
    CHECK(std::find(reqCrossings.begin(), reqCrossings.end(), cr) != reqCrossings.end());
    auto copyCrossings = reqCrossings;
    remove_value(copyCrossings, cr);
    CHECK(copyCrossings.size() + 1 == reqCrossings.size());
    if (isPlanarWithCrossings(graph, copyCrossings)) {
      reqCrossings = copyCrossings;
    }
  }
  // Nothing is done
  if (reqCrossings.size() == result.crossings.size())
    return;
  result.crossings = reqCrossings;
  result.isCrossed.assign(edges.size(), false);
  for (const auto& [e1, e2] : result.crossings) {
    result.isCrossed[e1] = true;
    result.isCrossed[e2] = true;
  }
  const int newNumCrossings = (int)result.crossings.size();
  LOG_IF(verbose, "reduced the number of crossings from %d to %d", orgNumCrossings, newNumCrossings);
}

struct BruteForceSolver {
  BruteForceSolver(const InputGraph& graph, const std::vector<std::pair<int, int>>& possibleCrossings, const int timeout) 
    : graph(graph),
      possibleCrossings(possibleCrossings),
      timeout(timeout)
    {}

  Result solve(int numCrossings) {
    isCrossed = std::vector<bool>(graph.edges.size(), false);
    crossings.clear();
    foundSolution = false;
    startTime = std::chrono::steady_clock::now();

    solveRec(0, numCrossings);

    if (foundSolution) {
      Result res(ResultCodeTy::SAT);
      res.crossings = crossings;
      res.isCrossed = isCrossed;
      return res;
    }

    if (timeoutReached) {
      return Result(ResultCodeTy::TIMEOUT);
    }

    return Result(ResultCodeTy::UNSAT);
  }

  void solveRec(size_t curIdx, int remainingCrossings) {
    if (foundSolution) return;
    if (timeoutReached) return;
    
    if (remainingCrossings == 0) {
      if (isPlanarWithCrossings(graph, crossings)) {
        foundSolution = true;
      } else if (isTimeout()) {
        timeoutReached = true;
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

    // try to get current
    const int e1 = possibleCrossings[curIdx].first;
    const int e2 = possibleCrossings[curIdx].second;
    if (!isCrossed[e1] && !isCrossed[e2]) {
      isCrossed[e1] = true;
      isCrossed[e2] = true;
      crossings.push_back(possibleCrossings[curIdx]);
      solveRec(curIdx + 1, remainingCrossings - 1);
      if (foundSolution) return;
      isCrossed[e1] = false;
      isCrossed[e2] = false;
      crossings.pop_back();
    }

    // skip current
    solveRec(curIdx + 1, remainingCrossings);
  }

  bool isTimeout() const {
    if (timeout == 0) return false;
    const auto curTime = std::chrono::steady_clock::now();
    const auto durationSec = std::chrono::duration_cast<std::chrono::seconds>(curTime - startTime).count();
    return durationSec > timeout;
  }

  const InputGraph& graph;
  const std::vector<std::pair<int, int>>& possibleCrossings;
  const int timeout = 0;
  std::chrono::time_point<std::chrono::steady_clock> startTime;
  bool timeoutReached = false;

  std::vector<bool> isCrossed;
  std::vector<std::pair<int, int>> crossings;
  bool foundSolution;
};

int computeSkewness(const InputGraph& graph, const int verbose, const int max_skewnees) {
  CHECK(!graph.isDirected());

  if (max_skewnees == 0)
    return max_skewnees + 1;

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

  for (int c = 1; c <= max_skewnees; c++) {
    LOG_IF(verbose >= 2, "  check skewness = %d; possibleCrossings = %d", c, possibleCrossings.size());
    BruteForceSolver solver(graph, possibleCrossings, 0);
    Result res = solver.solve(c);
    if (res.code == ResultCodeTy::SAT) {
      return c;
    } else {
      CHECK(res.code == ResultCodeTy::UNSAT, "unknown result code");
    }
  }

  return max_skewnees + 1;
}

Result bruteForce(const Params& params, const InputGraph& graph) {  
  CHECK(!graph.isDirected());
  CHECK(!params.useIC && !params.useNIC);

  const int timeout = params.timeout;
  const int n = graph.n;
  const int m = (int)graph.edges.size();
  const int maxCrossings = std::min(m / 2, n - 2);
  LOG_IF(params.verbose, "running brute-force with timeout=%d and maxCrossings=%d", timeout, maxCrossings);

  std::vector<std::pair<int, int>> possibleCrossings;
  for (int e1 = 0; e1 < m; e1++) {
    for (int e2 = e1 + 1; e2 < m; e2++) {
      if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
        continue;
      possibleCrossings.push_back({e1, e2});
    }
  }

  for (int c = 1; c <= maxCrossings; c++) {
    LOG_IF(params.verbose, "  check skewness = %d; possibleCrossings = %d", c, possibleCrossings.size());
    BruteForceSolver solver(graph, possibleCrossings, params.timeout);
    Result res = solver.solve(c);
    if (res.code == ResultCodeTy::SAT) {
      return res;
    } else if (res.code == ResultCodeTy::TIMEOUT) {
      return res;
    } else {
      if (res.code == ResultCodeTy::UNSAT && c == maxCrossings)
        return res;
      CHECK(res.code == ResultCodeTy::UNSAT, "unknown result code");
    }
  }

  return Result();
}