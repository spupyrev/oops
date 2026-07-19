#include "cmd_options.h"
#include "graph_algorithms.h"
#include "io.h"
#include "logging.h"
#include "one_planar.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <functional>
#include <memory>
#include <vector>

std::unique_ptr<GraphList> genGraphs(CMDOptions& options);

namespace {

/// Result of contracting one 6-cycle for the order-28 expansion proof.
struct SixCycleContraction {
  AdjListTy adjList;
  std::vector<EdgeTy> uncrossedAttachments;
};

enum class FiveCycleReduction { NO_FIVE_CYCLE, SUCCESS, FAILURE };

/// Test whether some 5-cycle can be replaced by a three-vertex path while
/// leaving a simple connected cubic graph of girth at least five.
FiveCycleReduction findFiveCycleReduction(const AdjListTy& adjList) {
  constexpr int cycleLength = 5;
  const int n = static_cast<int>(adjList.size());
  CHECK(n <= 63, "5-cycle path reduction uses a 64-bit vertex mask");
  for (const auto& neighbors : adjList)
    CHECK(neighbors.size() == 3, "5-cycle path reduction requires a cubic graph");
  std::vector<int> cycle;
  bool found = false;
  bool sawFiveCycle = false;

  auto tryReduction = [&](const int middlePosition) {
    uint64_t cycleMask = 0;
    for (int vertex : cycle)
      cycleMask |= uint64_t(1) << vertex;
    std::vector<int> outsideNeighbors;
    uint64_t outsideNeighborMask = 0;
    for (int position = 0; position < cycleLength; position++) {
      const int previous = cycle[(position + cycleLength - 1) % cycleLength];
      const int next = cycle[(position + 1) % cycleLength];
      int outsideNeighbor = -1;
      for (int neighbor : adjList[cycle[position]]) {
        if (neighbor == previous || neighbor == next)
          continue;
        if ((cycleMask & (uint64_t(1) << neighbor)) != 0 || outsideNeighbor != -1)
          return false;
        outsideNeighbor = neighbor;
      }
      if (outsideNeighbor == -1 ||
          (outsideNeighborMask & (uint64_t(1) << outsideNeighbor)) != 0)
        return false;
      outsideNeighborMask |= uint64_t(1) << outsideNeighbor;
      outsideNeighbors.push_back(outsideNeighbor);
    }

    std::vector<int> remap(n, -1);
    int reducedN = 0;
    for (int vertex = 0; vertex < n; vertex++)
      if ((cycleMask & (uint64_t(1) << vertex)) == 0)
        remap[vertex] = reducedN++;
    const int a = reducedN++;
    const int b = reducedN++;
    const int c = reducedN++;
    CHECK(reducedN == n - 2);

    AdjListTy reduced(reducedN);
    auto addEdge = [&](const int first, const int second) {
      if (first == second || contains(reduced[first], second))
        return false;
      reduced[first].push_back(second);
      reduced[second].push_back(first);
      return true;
    };
    for (int first = 0; first < n; first++) {
      if (remap[first] == -1)
        continue;
      for (int second : adjList[first])
        if (first < second && remap[second] != -1)
          CHECK(addEdge(remap[first], remap[second]));
    }
    CHECK(addEdge(a, b) && addEdge(b, c));
    const std::array<int, 5> attachmentOwners = {a, a, b, c, c};
    for (int offset = -2; offset <= 2; offset++) {
      const int position = (middlePosition + offset + cycleLength) % cycleLength;
      if (!addEdge(attachmentOwners[offset + 2], remap[outsideNeighbors[position]]))
        return false;
    }
    for (const auto& neighbors : reduced)
      if (neighbors.size() != 3)
        return false;
    return computeGirth(reducedN, adj_to_edges(reduced)) >= 5;
  };

  std::function<void(int, int)> enumerate = [&](const int start, const int current) {
    if (found)
      return;
    if (static_cast<int>(cycle.size()) == cycleLength) {
      if (cycle[1] > cycle.back() || !contains(adjList[current], start))
        return;
      sawFiveCycle = true;
      for (int middle = 0; middle < cycleLength && !found; middle++)
        found = tryReduction(middle);
      return;
    }
    for (int next : adjList[current]) {
      if (next <= start || contains(cycle, next))
        continue;
      cycle.push_back(next);
      enumerate(start, next);
      cycle.pop_back();
    }
  };
  for (int start = 0; start < n && !found; start++) {
    cycle = {start};
    enumerate(start, start);
  }
  if (found)
    return FiveCycleReduction::SUCCESS;
  return sawFiveCycle ? FiveCycleReduction::FAILURE : FiveCycleReduction::NO_FIVE_CYCLE;
}

/// Contract one chordless 6-cycle and return the four attachment edges that
/// the 6-cycle expansion lemma requires to remain uncrossed.
SixCycleContraction contractSixCycle(const AdjListTy& adjList) {
  constexpr int cycleLength = 6;
  const int n = static_cast<int>(adjList.size());
  CHECK(n <= 63, "6-cycle contraction uses a 64-bit vertex mask");
  for (const auto& neighbors : adjList)
    CHECK(neighbors.size() == 3, "6-cycle contraction is defined here only for cubic graphs");

  std::vector<int> outsideNeighbors;
  uint64_t selectedMask = 0;
  bool found = false;
  std::vector<int> cycle;
  std::function<void(int, int)> findCycle = [&](const int start, const int current) {
    if (found)
      return;
    if (static_cast<int>(cycle.size()) == cycleLength) {
      if (cycle[1] > cycle.back() || !contains(adjList[current], start))
        return;
      uint64_t cycleMask = 0;
      for (int vertex : cycle)
        cycleMask |= uint64_t(1) << vertex;
      std::vector<int> candidateOutsideNeighbors;
      uint64_t outsideNeighborMask = 0;
      for (int i = 0; i < cycleLength; i++) {
        const int previous = cycle[(i + cycleLength - 1) % cycleLength];
        const int next = cycle[(i + 1) % cycleLength];
        int outsideNeighbor = -1;
        for (int neighbor : adjList[cycle[i]]) {
          if (neighbor == previous || neighbor == next)
            continue;
          if ((cycleMask & (uint64_t(1) << neighbor)) != 0 || outsideNeighbor != -1)
            return;
          outsideNeighbor = neighbor;
        }
        if (outsideNeighbor == -1)
          return;
        const uint64_t outsideNeighborBit = uint64_t(1) << outsideNeighbor;
        if ((outsideNeighborMask & outsideNeighborBit) != 0)
          return;
        outsideNeighborMask |= outsideNeighborBit;
        candidateOutsideNeighbors.push_back(outsideNeighbor);
      }
      outsideNeighbors = std::move(candidateOutsideNeighbors);
      selectedMask = cycleMask;
      found = true;
      return;
    }
    for (int next : adjList[current]) {
      if (next <= start || contains(cycle, next))
        continue;
      cycle.push_back(next);
      findCycle(start, next);
      cycle.pop_back();
    }
  };
  for (int start = 0; start < n && !found; start++) {
    cycle = {start};
    findCycle(start, start);
  }
  if (!found)
    return {};

  // Rotate the cycle labeling before choosing s1, s2, s4, and s5 as the
  // attachment edges required to be uncrossed.  The expansion lemma is
  // invariant under this relabeling, while the fixed edge indices materially
  // affect SAT branching.
  std::rotate(outsideNeighbors.begin(), outsideNeighbors.end() - 1, outsideNeighbors.end());

  std::vector<int> remap(n, -1);
  int reducedN = 0;
  for (int vertex = 0; vertex < n; vertex++)
    if ((selectedMask & (uint64_t(1) << vertex)) == 0)
      remap[vertex] = reducedN++;
  const int contracted = reducedN++;
  auto mappedVertex = [&](const int vertex) {
    return (selectedMask & (uint64_t(1) << vertex)) != 0 ? contracted : remap[vertex];
  };

  std::vector<EdgeTy> edges;
  for (int u = 0; u < n; u++) {
    for (int v : adjList[u]) {
      if (u >= v)
        continue;
      const int mappedU = mappedVertex(u);
      const int mappedV = mappedVertex(v);
      if (mappedU != mappedV)
        edges.push_back(make_edge(mappedU, mappedV));
    }
  }
  const size_t edgesBeforeDeduplication = edges.size();
  sort_unique(edges);
  CHECK(edges.size() == edgesBeforeDeduplication, "6-cycle contraction created parallel edges");
  AdjListTy reduced = edges_to_adj(reducedN, edges);
  CHECK(reduced[contracted].size() == cycleLength, "contracted 6-cycle vertex has incorrect degree");

  std::vector<EdgeTy> uncrossedAttachments;
  for (int position : {1, 2, 4, 5}) {
    const int outsideNeighbor = outsideNeighbors[position];
    uncrossedAttachments.push_back(make_edge(contracted, remap[outsideNeighbor]));
  }
  return {std::move(reduced), std::move(uncrossedAttachments)};
}

/// Reconstruct the drawing represented by a satisfying assignment and check
/// independently that its proposed crossings produce a planarization.
Result verifiedDrawing(const Params& params, const InputGraph& graph, SATModel& model, Solver& solver) {
  Result drawing(ResultCodeTy::SAT);
  fillResultStack(model, solver, graph, params, drawing);
  CHECK(isPlanarWithCrossings(graph, drawing.crossings), "SAT produced an invalid 1-planar drawing");
  return drawing;
}

/// Find and validate one drawing in which all listed edges are uncrossed.
void verifyDrawing(
    const Params& params, const InputGraph& graph, const std::vector<int>& uncrossedEdges = {}) {
  initCrossablePairs(params, graph);

  SATModel model;
  encodeStackPlanar(model, graph, params);
  Solver solver;
  initSATSolver(params, graph, model, solver);
  vec<Lit> assumptions;
  for (int edge : uncrossedEdges) {
    CHECK(0 <= edge && edge < static_cast<int>(graph.edges.size()), "uncrossed-edge requirement is out of range");
    assumptions.push(model.getSolverLit(model.getCross1Var(graph.n + edge, false)));
  }
  CHECK(solver.solveLimited(assumptions) == l_True, "failed to find a required 1-planar drawing");

  const Result drawing = verifiedDrawing(params, graph, model, solver);
  for (int edge : uncrossedEdges)
    CHECK(!drawing.isCrossed[edge], "required edge %d is crossed", edge);
}

/// Return the sorted set of original edges crossed in a verified drawing.
std::vector<int> crossedEdges(const Result& drawing) {
  std::vector<int> crossed;
  for (const auto& [first, second] : drawing.crossings) {
    crossed.push_back(first);
    crossed.push_back(second);
  }
  sort_unique(crossed);
  return crossed;
}

/// Test whether a prescribed edge set is disjoint from a drawing's crossed-edge set.
bool areDisjoint(const std::vector<int>& prescribed, const std::vector<int>& crossed) {
  return std::none_of(prescribed.begin(), prescribed.end(), [&](const int edge) {
    return std::binary_search(crossed.begin(), crossed.end(), edge);
  });
}

/// Accumulate verified drawings until they prove k-flexibility or cover all
/// prescribed sets of edges that must be simultaneously uncrossed.
void verifyFlexibility(
    const Params& params, const InputGraph& graph, const int k,
    const std::vector<std::vector<int>>& prescribedSets = {}) {
  CHECK((k > 0) != !prescribedSets.empty(), "specify either k-flexibility or prescribed edge sets");

  initCrossablePairs(params, graph);
  SATModel model;
  encodeStackPlanar(model, graph, params);
  Solver solver;
  initSATSolver(params, graph, model, solver);
  CHECK(solveSATModel(params, model, solver) == l_True, "failed to find an initial 1-planar drawing");

  Result drawing = verifiedDrawing(params, graph, model, solver);
  std::vector<std::vector<int>> crossedEdgeSets = {crossedEdges(drawing)};

  auto findUncoveredEdgeSet = [&]() {
    if (!prescribedSets.empty()) {
      for (const auto& prescribed : prescribedSets) {
        const bool isCovered =
            std::any_of(crossedEdgeSets.begin(), crossedEdgeSets.end(),
                        [&](const std::vector<int>& crossed) { return areDisjoint(prescribed, crossed); });
        if (!isCovered)
          return prescribed;
      }
      return std::vector<int>();
    }

    std::vector<int> hittingSet;
    std::function<bool(int)> findHittingSet = [&](const int limit) {
      const std::vector<int>* unhitCrossedSet = nullptr;
      for (const auto& crossed : crossedEdgeSets) {
        const bool hit = std::any_of(hittingSet.begin(), hittingSet.end(), [&](const int edge) {
          return std::binary_search(crossed.begin(), crossed.end(), edge);
        });
        if (!hit) {
          unhitCrossedSet = &crossed;
          break;
        }
      }
      if (unhitCrossedSet == nullptr)
        return true;
      if (static_cast<int>(hittingSet.size()) == limit || unhitCrossedSet->empty())
        return false;
      for (int edge : *unhitCrossedSet) {
        hittingSet.push_back(edge);
        if (findHittingSet(limit))
          return true;
        hittingSet.pop_back();
      }
      return false;
    };

    for (int limit = 1; limit <= k; limit++) {
      hittingSet.clear();
      if (findHittingSet(limit))
        return hittingSet;
    }
    return std::vector<int>();
  };

  while (true) {
    const std::vector<int> prescribed = findUncoveredEdgeSet();
    if (prescribed.empty())
      return;

    vec<Lit> assumptions;
    for (int edge : prescribed)
      assumptions.push(model.getSolverLit(model.getCross1Var(graph.n + edge, false)));
    CHECK(solver.solveLimited(assumptions) == l_True, "failed uncrossed-edge requirement");

    drawing = verifiedDrawing(params, graph, model, solver);
    for (int edge : prescribed)
      CHECK(!drawing.isCrossed[edge], "required edge %d is crossed", edge);
    crossedEdgeSets.push_back(crossedEdges(drawing));
  }
}

/// Convert endpoint pairs into sorted edge indices of the reduced graph.
std::vector<int> edgeIndices(const InputGraph& graph, const std::vector<EdgeTy>& edges) {
  std::vector<int> indices;
  for (const auto& [u, v] : edges) {
    const int edge = graph.findEdgeIndex(u, v);
    CHECK(edge != -1, "required edge (%d,%d) is absent", u, v);
    indices.push_back(edge);
  }
  sort_unique(indices);
  return indices;
}

/// Check that an input record belongs to the required graph family for its order.
void validateInput(const AdjListTy& adjList, const std::vector<EdgeTy>& edges) {
  const int n = static_cast<int>(adjList.size());
  CHECK(4 <= n && n <= 28 && n % 2 == 0, "verification expects an even order from 4 through 28");
  CHECK(std::all_of(adjList.begin(), adjList.end(), [](const auto& neighbors) { return neighbors.size() == 3; }),
        "verification expects a cubic graph");
  CHECK(isConnected(n, edges), "verification inputs contain connected graphs");
  if (n >= 24) {
    CHECK(computeGirth(n, edges) >= 5, "orders 24, 26, and 28 require girth at least five");
  }
}

/// Verify Claim 1 for a nonplanar input through order 22.
void verifyCubic22(const int n, const std::vector<EdgeTy>& edges, const Params& params) {
  CHECK(n <= 22);
  InputGraph graph(n, edges, {});
  verifyFlexibility(params, graph, 3);
}

/// Verify Claim 2 for a nonplanar order-24 input.
void verifyCubic24(const int n, const std::vector<EdgeTy>& edges, const Params& params) {
  CHECK(n == 24);
  InputGraph graph(n, edges, {});
  verifyFlexibility(params, graph, 2);
}

/// Verify Claim 3: every vertex has a drawing with its three incident edges uncrossed.
void verifyCubic26(const int n, const std::vector<EdgeTy>& edges, const Params& params) {
  CHECK(n == 26);
  InputGraph graph(n, edges, {});
  std::vector<std::vector<int>> incidentEdges(n);
  for (int edge = 0; edge < static_cast<int>(edges.size()); edge++) {
    incidentEdges[edges[edge].first].push_back(edge);
    incidentEdges[edges[edge].second].push_back(edge);
  }
  verifyFlexibility(params, graph, 0, incidentEdges);
}

/// Method used to verify a nonplanar order-28 input.
enum class Cubic28Method {
  FIVE_CYCLE_EXPANSION,
  DIRECT_GIRTH_FIVE,
  SIX_CYCLE_EXPANSION,
  DIRECT_NO_SIX_CYCLE
};

/// Verify Claims 4 and 5 by a 5-cycle expansion, a 6-cycle expansion, or a direct drawing.
Cubic28Method verifyCubic28(
    const int n, const std::vector<EdgeTy>& edges, const AdjListTy& adj, const Params& params) {
  CHECK(n == 28);
  const FiveCycleReduction fiveCycleReduction = findFiveCycleReduction(adj);
  if (fiveCycleReduction == FiveCycleReduction::SUCCESS)
    return Cubic28Method::FIVE_CYCLE_EXPANSION;
  if (fiveCycleReduction == FiveCycleReduction::FAILURE) {
    InputGraph graph(n, edges, {});
    verifyDrawing(params, graph);
    return Cubic28Method::DIRECT_GIRTH_FIVE;
  }

  const SixCycleContraction sixCycle = contractSixCycle(adj);
  if (!sixCycle.adjList.empty()) {
    InputGraph graph(static_cast<int>(sixCycle.adjList.size()), adj_to_edges(sixCycle.adjList), {});
    CHECK(sixCycle.uncrossedAttachments.size() == 4, "invalid 6-cycle contraction");
    if (!isPlanar(graph.n, graph.edges, 0)) {
      verifyDrawing(params, graph, edgeIndices(graph, sixCycle.uncrossedAttachments));
    }
    return Cubic28Method::SIX_CYCLE_EXPANSION;
  }

  InputGraph graph(n, edges, {});
  verifyDrawing(params, graph);
  return Cubic28Method::DIRECT_NO_SIX_CYCLE;
}

}  // namespace

/// Verify the claim or claims prescribed by the order of each input graph.
void verifyCubic(CMDOptions& options) {
  Params params;
  params.cubicVerification = true;
  params.crossPriority = true;

  auto graphs = genGraphs(options);
  const int numGraphs = graphs->size();
  const auto start = std::chrono::steady_clock::now();

  int numPlanar = 0;
  int num1Planar = 0;
  int num3Flexible = 0;
  int num2Flexible = 0;
  int numClaim3Verified = 0;
  int numFiveCycleExpanded = 0;
  int numFiveCycleDirect = 0;
  int numSixCycleExpanded = 0;
  int numNoSixCycleDirect = 0;

  for (int i = 0; i < numGraphs; i++) {
    const auto& adj = graphs->next().second;
    const int n = static_cast<int>(adj.size());
    const std::vector<EdgeTy> edges = adj_to_edges(adj);
    validateInput(adj, edges);

    if (isPlanar(n, edges, 0)) {
      numPlanar++;
    } else if (n <= 22) {
      verifyCubic22(n, edges, params);
      num3Flexible++;
      num1Planar++;
    } else if (n == 24) {
      verifyCubic24(n, edges, params);
      num2Flexible++;
      num1Planar++;
    } else if (n == 26) {
      verifyCubic26(n, edges, params);
      numClaim3Verified++;
      num1Planar++;
    } else {
      const Cubic28Method method = verifyCubic28(n, edges, adj, params);
      if (method == Cubic28Method::FIVE_CYCLE_EXPANSION)
        numFiveCycleExpanded++;
      else if (method == Cubic28Method::DIRECT_GIRTH_FIVE)
        numFiveCycleDirect++;
      else if (method == Cubic28Method::SIX_CYCLE_EXPANSION)
        numSixCycleExpanded++;
      else
        numNoSixCycleDirect++;
      num1Planar++;
    }

    LOG_EVERY_MS(30000, "verified %'d of %'d graphs", i + 1, numGraphs);
  }

  LOG("processed %'d graphs in %s", numGraphs, ms_to_str(start, std::chrono::steady_clock::now()).c_str());
  LOG("#planar = %'d; #1-planar = %'d", numPlanar, num1Planar);
  LOG("#3-flexible = %'d; #2-flexible = %'d", num3Flexible, num2Flexible);
  LOG("#claim3-verified = %'d", numClaim3Verified);
  LOG("#5-cycle-expanded = %'d; #direct-girth5 = %'d; #6-cycle-expanded = %'d; "
      "#direct-no-6-cycle = %'d",
      numFiveCycleExpanded, numFiveCycleDirect, numSixCycleExpanded, numNoSixCycleDirect);
}
