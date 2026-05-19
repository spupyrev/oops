#include "one_planar.h"
#include "graph_algorithms.h"
#include "logging.h"

#include <algorithm>
#include <cstdint>
#include <unordered_map>

using namespace Simp21;

namespace {

class SepCyclesDynamic final : public UserPropagator {
public:
  SepCyclesDynamic(const SATModel& model,
                   const InputGraph& graph,
                   const Params& params)
      : model(model),
        graph(graph),
        n(graph.n),
        m((int)graph.edges.size()),
        verbose(params.verbose) {
    initPossibleCrossings();
    LOG_IF(verbose, "sep-cycle UP watches %d cross2 and %d cross1 variables", 
           crossings.size(), m);
    for (int crossingId = 0; crossingId < (int)this->crossings.size(); crossingId++) {
      crossingIdByLit[toInt(cross2Lit(crossingId))] = crossingId;
    }
    for (int edgeId = 0; edgeId < m; edgeId++) {
      uncrossedEdgeIdByLit[toInt(~cross1Lit(edgeId))] = edgeId;
    }
    initStaticState();
  }

  ~SepCyclesDynamic() override {
    LOG_IF(verbose,
           "sep-cycle UP stats: calls=%'llu; clean-skips=%'llu; pending-calls=%'llu; "
           "small-assignment-skips=%'llu; checked=%'llu; conflicts=%'llu",
           numPropagationCalls, numCleanPropagationCalls, numPendingPropagationCalls,
           numSmallAssignmentSkips, numChecks, numConflicts);
  }

  void onAssignment(Lit p, int decisionLevel) override {
    const int litKey = toInt(p);

    const auto crossingIt = crossingIdByLit.find(litKey);
    if (crossingIt != crossingIdByLit.end()) {
      activeCrossingIds.push_back(crossingIt->second);
      activeCrossingDecisionLevels.push_back(decisionLevel);
      pendingCrossingIds.push_back(crossingIt->second);
      pendingCrossingDecisionLevels.push_back(decisionLevel);
      return;
    }

    const auto uncrossedEdgeIt = uncrossedEdgeIdByLit.find(litKey);
    if (uncrossedEdgeIt != uncrossedEdgeIdByLit.end()) {
      const int edgeId = uncrossedEdgeIt->second;
      setUncrossableEdge(edgeId, true);
      activeUncrossedEdgeIds.push_back(edgeId);
      activeUncrossedEdgeDecisionLevels.push_back(decisionLevel);
    }
  }

  void onBacktrack(int decisionLevel) override {
    removeBacktracked(activeCrossingIds, activeCrossingDecisionLevels, decisionLevel);
    removeBacktracked(pendingCrossingIds, pendingCrossingDecisionLevels, decisionLevel);
    removeBacktrackedUncrossedEdges(decisionLevel);
  }

  bool findConflict(std::vector<Lit>& clause) override {
    clause.clear();

    numPropagationCalls++;
    if (pendingCrossingIds.empty()) {
      numCleanPropagationCalls++;
      return false;
    }
    CHECK(!activeCrossingIds.empty());
    numPendingPropagationCalls++;

    const int numActiveCross2 = (int)activeCrossingIds.size();
    const int numActiveCross1False = (int)activeUncrossedEdgeIds.size();
    const bool smallAssignment = numActiveCross2 + numActiveCross1False <= 2;
    if (smallAssignment) {
      numSmallAssignmentSkips++;
      clearPendingCrossings();
      return false;
    }

    numChecks++;

    std::vector<Lit> reason;
    if (findCurrentAssignmentClause(reason)) {
      sort_unique(reason);
      CHECK(reason.size() >= 2,
            "sep-cycle UP found a unit/empty reason; expected at least a binary conflict clause");
      CHECK(reason.size() >= 3 || crossings.size() > MAX_CROSSINGS_FOR_2GT,
            "sep-cycle UP found a binary reason even though static sep-cycle 2-clauses should be complete");

      clause = std::move(reason);

      LOG_IF(verbose >= 2, "sep-cycle UP conflict with %d crossings", clause.size());
      numConflicts++;
      return true;
    }

    clearPendingCrossings();
    return false;
  }

private:
  static constexpr int MAX_DYNAMIC_CYCLE_LEN = 8;

  void initPossibleCrossings() {
    for (int e1 = 0; e1 < m; e1++) {
      for (int e2 = e1 + 1; e2 < m; e2++) {
        if (!canBeMerged(e1 + n, e2 + n, n, graph.edges))
          continue;
        crossings.push_back({e1, e2});
      }
    }
  }

  Lit cross2Lit(int crossingId) const {
    const auto [e1, e2] = crossings[crossingId];
    return model.getSolverLit(model.getCross2Var(e1 + n, e2 + n, true));
  }

  Lit cross1Lit(int edgeId) const {
    return model.getSolverLit(model.getCross1Var(edgeId + n, true));
  }

  void clearPendingCrossings() {
    pendingCrossingIds.clear();
    pendingCrossingDecisionLevels.clear();
  }

  void removeBacktracked(std::vector<int>& ids, std::vector<int>& decisionLevels, int decisionLevel) {
    CHECK(ids.size() == decisionLevels.size());
    int write = 0;
    for (int read = 0; read < (int)ids.size(); read++) {
      if (decisionLevels[read] <= decisionLevel) {
        ids[write] = ids[read];
        decisionLevels[write] = decisionLevels[read];
        write++;
      }
    }
    ids.resize(write);
    decisionLevels.resize(write);
  }

  void removeBacktrackedUncrossedEdges(int decisionLevel) {
    CHECK(activeUncrossedEdgeIds.size() == activeUncrossedEdgeDecisionLevels.size());
    int write = 0;
    for (int read = 0; read < (int)activeUncrossedEdgeIds.size(); read++) {
      const int edgeId = activeUncrossedEdgeIds[read];
      const int edgeDecisionLevel = activeUncrossedEdgeDecisionLevels[read];
      if (edgeDecisionLevel > decisionLevel) {
        setUncrossableEdge(edgeId, false);
      } else {
        activeUncrossedEdgeIds[write] = edgeId;
        activeUncrossedEdgeDecisionLevels[write] = edgeDecisionLevel;
        write++;
      }
    }
    activeUncrossedEdgeIds.resize(write);
    activeUncrossedEdgeDecisionLevels.resize(write);
  }

  void initStaticState() {
    isEdge.assign(n, std::vector<bool>(n, false));
    adjList.assign(n, std::vector<int>());
    isCrossEdge.assign(n, std::vector<bool>(n, false));
    isKiteEdge.assign(n, std::vector<int>(n, 0));
    isUncrossableEdge.assign(n, std::vector<bool>(n, false));
    uncrossableEdgeCount.assign(m, 0);
    crossEdgeProvider.assign(n * n, -1);
    kiteEdgeProviders.assign(n * n, std::vector<int>());
    edgeIdByKey.assign(n * n, -1);
    cycleBlocked.assign(n, 0);
    cycleUsed.assign(n, 0);
    cycleBuffer.reserve(MAX_DYNAMIC_CYCLE_LEN);

    for (int edgeId = 0; edgeId < m; edgeId++) {
      const auto& [u, v] = graph.edges[edgeId];
      CHECK(u < v, "graph edge (%d, %d) is not normalized", u, v);
      isEdge[u][v] = true;
      isEdge[v][u] = true;
      adjList[u].push_back(v);
      adjList[v].push_back(u);
      edgeIdByKey[edgeKey(u, v)] = edgeId;
    }
  }

  int edgeKey(int u, int v) const {
    if (u > v) {
      return v * n + u;
    }
    return u * n + v;
  }

  bool canMarkCrossing(int u, int v, int x, int y) const {
    if (isCrossEdge[u][v] || isCrossEdge[x][y])
      return false;
    if (isKiteEdge[u][v] || isKiteEdge[x][y])
      return false;
    if (isCrossEdge[u][x] || isCrossEdge[x][v] || isCrossEdge[v][y] || isCrossEdge[y][u])
      return false;
    return true;
  }

  void setUncrossableEdge(int edgeId, bool value) {
    const auto [u, v] = graph.edges[edgeId];
    CHECK(isEdge[u][v]);
    if (value) {
      uncrossableEdgeCount[edgeId]++;
    } else {
      CHECK(uncrossableEdgeCount[edgeId] > 0);
      uncrossableEdgeCount[edgeId]--;
    }
    const bool isUncrossable = uncrossableEdgeCount[edgeId] > 0;
    isUncrossableEdge[u][v] = isUncrossable;
    isUncrossableEdge[v][u] = isUncrossable;
  }

  void markKiteEdge(int crossingId, int u, int v) {
    isKiteEdge[u][v] += 1;
    isKiteEdge[v][u] += 1;
    kiteEdgeProviders[edgeKey(u, v)].push_back(crossingId);

    if (isKiteEdge[u][v] == 1 && !isEdge[u][v]) {
      adjList[u].push_back(v);
      adjList[v].push_back(u);
    }
  }

  void unmarkKiteEdge(int crossingId, int u, int v) {
    CHECK(isKiteEdge[u][v] > 0);
    CHECK(isKiteEdge[u][v] == isKiteEdge[v][u]);
    CHECK(kiteEdgeProviders[edgeKey(u, v)].back() == crossingId);

    isKiteEdge[u][v] -= 1;
    isKiteEdge[v][u] -= 1;
    kiteEdgeProviders[edgeKey(u, v)].pop_back();

    if (isKiteEdge[u][v] == 0 && !isEdge[u][v]) {
      CHECK(adjList[u].back() == v);
      adjList[u].pop_back();
      CHECK(adjList[v].back() == u);
      adjList[v].pop_back();
    }
  }

  void markCrossing(int crossingId, int u, int v, int x, int y) {
    isCrossEdge[u][v] = true;
    isCrossEdge[v][u] = true;
    isCrossEdge[x][y] = true;
    isCrossEdge[y][x] = true;
    crossEdgeProvider[edgeKey(u, v)] = crossingId;
    crossEdgeProvider[edgeKey(x, y)] = crossingId;

    // Same orientation convention as SepCyclesStatic: u--x--v--y.
    markKiteEdge(crossingId, u, x);
    markKiteEdge(crossingId, x, v);
    markKiteEdge(crossingId, v, y);
    markKiteEdge(crossingId, y, u);
  }

  void unmarkCrossing(int crossingId) {
    const auto [e1, e2] = crossings[crossingId];
    const auto [u, v] = graph.edges[e1];
    const auto [x, y] = graph.edges[e2];

    CHECK(isCrossEdge[u][v] && isCrossEdge[x][y]);
    CHECK(crossEdgeProvider[edgeKey(u, v)] == crossingId);
    CHECK(crossEdgeProvider[edgeKey(x, y)] == crossingId);

    isCrossEdge[u][v] = false;
    isCrossEdge[v][u] = false;
    isCrossEdge[x][y] = false;
    isCrossEdge[y][x] = false;
    crossEdgeProvider[edgeKey(u, v)] = -1;
    crossEdgeProvider[edgeKey(x, y)] = -1;

    unmarkKiteEdge(crossingId, y, u);
    unmarkKiteEdge(crossingId, v, y);
    unmarkKiteEdge(crossingId, x, v);
    unmarkKiteEdge(crossingId, u, x);
  }

  bool isFreeEdge(int u, int v) const {
    return !isCrossEdge[u][v] && !isKiteEdge[u][v] && !isUncrossableEdge[u][v];
  }

  bool isCurrentEdge(int u, int v) const {
    return isEdge[u][v] || isKiteEdge[u][v] > 0;
  }

  bool isEdgeOnCycle(const std::vector<int>& cycle, const EdgeTy& edge) const {
    for (size_t i = 0; i < cycle.size(); i++) {
      if (edge == make_edge(cycle[i], cycle[(i + 1) % cycle.size()])) {
        return true;
      }
    }
    return false;
  }

  void addCrossingReason(std::vector<Lit>& reason, int crossingId) const {
    if (crossingId >= 0) {
      reason.push_back(~cross2Lit(crossingId));
    }
  }

  void addUncrossableEdgeReason(std::vector<Lit>& reason, int u, int v) const {
    if (!isUncrossableEdge[u][v] || isCrossEdge[u][v] || isKiteEdge[u][v]) {
      return;
    }
    const int edgeId = edgeIdByKey[edgeKey(u, v)];
    CHECK(edgeId >= 0);
    reason.push_back(cross1Lit(edgeId));
  }

  std::vector<EdgeTy> getCrossCycleEdges(const std::vector<int>& cycle, std::vector<Lit>& reason) const {
    std::vector<EdgeTy> crossCycleEdges;
    for (int crossingId : activeCrossingIds) {
      const auto [e1, e2] = crossings[crossingId];
      const EdgeTy edge1 = graph.edges[e1];
      const EdgeTy edge2 = graph.edges[e2];
      const bool edge1OnCycle = isEdgeOnCycle(cycle, edge1);
      const bool edge2OnCycle = isEdgeOnCycle(cycle, edge2);
      if (edge1OnCycle) {
        crossCycleEdges.push_back(edge2);
        addCrossingReason(reason, crossingId);
      }
      if (edge2OnCycle) {
        crossCycleEdges.push_back(edge1);
        addCrossingReason(reason, crossingId);
      }
    }
    sort_unique(crossCycleEdges);
    return crossCycleEdges;
  }

  void addCycleProviders(const std::vector<int>& cycle, std::vector<Lit>& reason) const {
    for (size_t i = 0; i < cycle.size(); i++) {
      const int u = cycle[i];
      const int v = cycle[(i + 1) % cycle.size()];
      const int key = edgeKey(u, v);

      addCrossingReason(reason, crossEdgeProvider[key]);
      for (int reasonCrossingId : kiteEdgeProviders[key]) {
        addCrossingReason(reason, reasonCrossingId);
      }
      addUncrossableEdgeReason(reason, u, v);
    }
  }

  bool processCycle(int mainCrossingId, int u, int v, size_t numFree,
                    std::vector<Lit>& reason) const {
    CHECK(cycleBuffer.size() >= 3);

    if (numFree >= (size_t)std::min(graph.degree(u) - 1, graph.degree(v) - 1)) {
      return false;
    }

    std::vector<Lit> candidateReason;
    addCrossingReason(candidateReason, mainCrossingId);
    addCycleProviders(cycleBuffer, candidateReason);
    const auto crossCycleEdges = getCrossCycleEdges(cycleBuffer, candidateReason);
    const size_t numPaths = countEdgeDisjointPaths({u}, {v}, graph.adj, cycleBuffer, crossCycleEdges, numFree + 1);

    if (numPaths > numFree) {
      CHECK((int)activeCrossingIds.size() + (int)activeUncrossedEdgeIds.size() > 2,
            "sep-cycle found with at most two active cross1=false/cross2 assignments");
      reason = std::move(candidateReason);
      return true;
    }

    return false;
  }

  bool findSepCycleForCrossing(int mainCrossingId, int x, int y, int u, int v,
                               int maxCycleLen, std::vector<Lit>& reason) {
    const size_t degreeBound = size_t(std::min(graph.degree(u) - 1, graph.degree(v) - 1));
    const int maxDepth = maxCycleLen - 2;

    cycleBlocked[x] = cycleBlocked[y] = cycleBlocked[u] = cycleBlocked[v] = 1;
    cycleUsed[x] = 1;
    cycleBuffer.clear();
    cycleBuffer.push_back(x);

    auto dfs = [&](auto&& self, int current, int depth, size_t numFree) -> bool {
      if (depth >= 1 && isCurrentEdge(current, y)) {
        const size_t finalFree = numFree + size_t(isFreeEdge(current, y));
        if (finalFree < degreeBound) {
          cycleBuffer.push_back(y);
          const bool stop = processCycle(mainCrossingId, u, v, finalFree, reason);
          cycleBuffer.pop_back();
          if (stop)
            return true;
        }
      }

      if (depth == maxDepth)
        return false;

      for (int next : adjList[current]) {
        if (cycleBlocked[next] || cycleUsed[next])
          continue;

        const size_t nextFree = numFree + size_t(isFreeEdge(current, next));
        if (nextFree >= degreeBound)
          continue;

        cycleUsed[next] = 1;
        cycleBuffer.push_back(next);
        const bool stop = self(self, next, depth + 1, nextFree);
        cycleBuffer.pop_back();
        cycleUsed[next] = 0;
        if (stop)
          return true;
      }
      return false;
    };

    const bool found = dfs(dfs, x, 0, size_t(isFreeEdge(y, x)));

    cycleBuffer.clear();
    cycleUsed[x] = 0;
    cycleBlocked[x] = cycleBlocked[y] = cycleBlocked[u] = cycleBlocked[v] = 0;

    return found;
  }

  bool findSepCycleClause(int maxCycleLen, std::vector<Lit>& reason) {
    for (int crossingId : pendingCrossingIds) {
      const auto [e1, e2] = crossings[crossingId];
      const auto [u, v] = graph.edges[e1];
      const auto [x, y] = graph.edges[e2];

      if (findSepCycleForCrossing(crossingId, x, y, u, v, maxCycleLen, reason)) {
        return true;
      }
      if (findSepCycleForCrossing(crossingId, u, v, x, y, maxCycleLen, reason)) {
        return true;
      }
    }
    return !reason.empty();
  }

  bool findCurrentAssignmentClause(std::vector<Lit>& reason) {
    std::vector<int> markedCrossings;
    bool incompatibleCrossings = false;
    for (int crossingId : activeCrossingIds) {
      const auto [e1, e2] = crossings[crossingId];
      const auto [u, v] = graph.edges[e1];
      const auto [x, y] = graph.edges[e2];
      if (!canMarkCrossing(u, v, x, y)) {
        incompatibleCrossings = true;
        for (int reasonCrossingId : activeCrossingIds) {
          addCrossingReason(reason, reasonCrossingId);
        }
        break;
      }
      markCrossing(crossingId, u, v, x, y);
      markedCrossings.push_back(crossingId);
    }

    if (reason.empty() && !incompatibleCrossings) {
      findSepCycleClause(MAX_DYNAMIC_CYCLE_LEN, reason);
    }

    for (auto it = markedCrossings.rbegin(); it != markedCrossings.rend(); ++it) {
      unmarkCrossing(*it);
    }

    return !reason.empty();
  }

private:
  const SATModel& model;
  const InputGraph& graph;
  const int n;
  const int m;
  const int verbose;
  std::vector<std::pair<int, int>> crossings;
  std::unordered_map<int, int> crossingIdByLit;
  std::unordered_map<int, int> uncrossedEdgeIdByLit;
  uint64_t numPropagationCalls = 0;
  uint64_t numCleanPropagationCalls = 0;
  uint64_t numPendingPropagationCalls = 0;
  uint64_t numSmallAssignmentSkips = 0;
  uint64_t numChecks = 0;
  uint64_t numConflicts = 0;

  std::vector<int> activeCrossingIds;
  std::vector<int> activeCrossingDecisionLevels;
  std::vector<int> pendingCrossingIds;
  std::vector<int> pendingCrossingDecisionLevels;
  std::vector<int> activeUncrossedEdgeIds;
  std::vector<int> activeUncrossedEdgeDecisionLevels;
  std::vector<std::vector<bool>> isEdge;
  std::vector<std::vector<int>> adjList;
  std::vector<std::vector<bool>> isCrossEdge;
  std::vector<std::vector<int>> isKiteEdge;
  std::vector<std::vector<bool>> isUncrossableEdge;
  std::vector<int> uncrossableEdgeCount;
  std::vector<int> crossEdgeProvider;
  std::vector<std::vector<int>> kiteEdgeProviders;
  std::vector<int> edgeIdByKey;
  std::vector<char> cycleBlocked;
  std::vector<char> cycleUsed;
  std::vector<int> cycleBuffer;
};

} // namespace

std::unique_ptr<UserPropagator> createSepCyclesDynamic(
    const SATModel& model,
    const InputGraph& graph,
    const Params& params) {
  CHECK(params.solverType == SolverType::STACK, "sep-cycle user propagation supports stack-planarity only");
  CHECK(params.useSATConstraints, "sep-cycle user propagation requires cross variables");
  CHECK(params.useUNSATConstraints, "sep-cycle user propagation requires -unsat=1");

  return std::make_unique<SepCyclesDynamic>(model, graph, params);
}
