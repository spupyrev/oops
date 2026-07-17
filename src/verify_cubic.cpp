#include "cmd_options.h"
#include "graph_algorithms.h"
#include "io.h"
#include "logging.h"
#include "one_planar.h"

#include <algorithm>
#include <chrono>
#include <functional>
#include <memory>
#include <vector>

std::unique_ptr<GraphList> genGraphs(CMDOptions& options);

namespace {

/// Result of contracting a compatible packing of pentagons.
struct PackedPentagonContraction {
  AdjListTy adjList;
  std::vector<EdgeTy> protectedAttachments;
};

/// Result of contracting one pentagon for the order-26 flexibility proof.
struct PentagonFlexibilityContraction {
  AdjListTy adjList;
  std::vector<std::vector<EdgeTy>> requirements;
};

/// Result of contracting one hexagon for the order-28 expansion proof.
struct HexagonContraction {
  AdjListTy adjList;
  std::vector<EdgeTy> protectedAttachments;
};

/// Contract a deterministic compatible packing of chordless pentagons and
/// select two consecutive protected attachment edges for each pentagon.
PackedPentagonContraction contractPackedPentagons(const AdjListTy& adjList) {
  constexpr int cycleLength = 5;
  constexpr int protectedSpokes = 2;
  const int n = static_cast<int>(adjList.size());
  CHECK(n <= 63, "pentagon packing uses a 64-bit vertex mask");
  for (const auto& neighbors : adjList)
    CHECK(neighbors.size() == 3, "pentagon contraction is defined here only for cubic graphs");

  struct Candidate {
    std::vector<int> vertices;
    std::vector<int> outside;
    uint64_t mask = 0;
  };
  std::vector<Candidate> candidates;
  std::vector<int> cycle;
  std::function<void(int, int)> enumerate = [&](const int start, const int current) {
    if (static_cast<int>(cycle.size()) == cycleLength) {
      if (cycle[1] > cycle.back() ||
          std::find(adjList[current].begin(), adjList[current].end(), start) == adjList[current].end())
        return;

      Candidate candidate;
      candidate.vertices = cycle;
      for (int vertex : cycle)
        candidate.mask |= uint64_t(1) << vertex;
      uint64_t outsideMask = 0;
      for (int i = 0; i < cycleLength; i++) {
        const int previous = cycle[(i + cycleLength - 1) % cycleLength];
        const int next = cycle[(i + 1) % cycleLength];
        int external = -1;
        for (int neighbor : adjList[cycle[i]]) {
          if (neighbor == previous || neighbor == next)
            continue;
          if ((candidate.mask & (uint64_t(1) << neighbor)) != 0 || external != -1)
            return;
          external = neighbor;
        }
        if (external == -1)
          return;
        const uint64_t externalBit = uint64_t(1) << external;
        if ((outsideMask & externalBit) != 0)
          return;
        outsideMask |= externalBit;
        candidate.outside.push_back(external);
      }
      candidates.push_back(std::move(candidate));
      return;
    }

    for (int next : adjList[current]) {
      if (next <= start || std::find(cycle.begin(), cycle.end(), next) != cycle.end())
        continue;
      cycle.push_back(next);
      enumerate(start, next);
      cycle.pop_back();
    }
  };
  for (int start = 0; start < n; start++) {
    cycle = {start};
    enumerate(start, start);
  }
  if (candidates.empty())
    return {};

  auto compatible = [&](const int first, const int second) {
    if ((candidates[first].mask & candidates[second].mask) != 0)
      return false;
    int joiningEdges = 0;
    for (int outside : candidates[first].outside)
      joiningEdges += (candidates[second].mask & (uint64_t(1) << outside)) != 0;
    return joiningEdges <= 1;
  };

  auto hasProtectedPairs = [&](const std::vector<int>& packing) {
    uint64_t packedVertices = 0;
    for (int candidate : packing)
      packedVertices |= candidates[candidate].mask;
    for (int candidate : packing) {
      bool found = false;
      const uint64_t otherPentagons = packedVertices & ~candidates[candidate].mask;
      for (int offset = 0; offset < cycleLength; offset++) {
        bool clean = true;
        for (int delta = 0; delta < protectedSpokes; delta++)
          clean &= (otherPentagons &
                    (uint64_t(1) << candidates[candidate].outside[(offset + delta) % cycleLength])) == 0;
        found |= clean;
      }
      if (!found)
        return false;
    }
    return true;
  };

  // Choose a deterministic, low-conflict maximal packing.  Maximality is not
  // needed for correctness; contracting more pentagons only improves speed.
  std::vector<int> available = identity(candidates.size());
  std::vector<int> selected;
  while (!available.empty()) {
    std::vector<std::pair<int, int>> ranked;
    for (int candidate : available) {
      bool fits = true;
      for (int chosen : selected)
        fits &= compatible(candidate, chosen);
      std::vector<int> tentative = selected;
      tentative.push_back(candidate);
      if (!fits || !hasProtectedPairs(tentative))
        continue;
      int conflicts = 0;
      for (int other : available)
        conflicts += candidate != other && !compatible(candidate, other);
      ranked.push_back({conflicts, candidate});
    }
    if (ranked.empty())
      break;
    std::sort(ranked.begin(), ranked.end());
    const int chosen = ranked.front().second;
    selected.push_back(chosen);
    available.erase(std::remove_if(available.begin(), available.end(), [&](const int candidate) {
      return candidate == chosen || !compatible(chosen, candidate);
    }), available.end());
  }
  CHECK(!selected.empty());

  std::vector<int> owner(n, -1);
  for (int packed = 0; packed < static_cast<int>(selected.size()); packed++)
    for (int vertex : candidates[selected[packed]].vertices)
      owner[vertex] = packed;

  std::vector<int> remap(n, -1);
  int reducedN = 0;
  for (int vertex = 0; vertex < n; vertex++)
    if (owner[vertex] == -1)
      remap[vertex] = reducedN++;
  std::vector<int> contracted(selected.size());
  for (int& vertex : contracted)
    vertex = reducedN++;
  auto mappedVertex = [&](const int vertex) {
    return owner[vertex] == -1 ? remap[vertex] : contracted[owner[vertex]];
  };

  std::vector<EdgeTy> edges;
  for (int u = 0; u < n; u++) {
    for (int v : adjList[u]) {
      if (u >= v)
        continue;
      const int mappedU = mappedVertex(u);
      const int mappedV = mappedVertex(v);
      if (mappedU != mappedV)
        edges.push_back({std::min(mappedU, mappedV), std::max(mappedU, mappedV)});
    }
  }
  const size_t edgesBeforeDeduplication = edges.size();
  sort_unique(edges);
  CHECK(edges.size() == edgesBeforeDeduplication, "pentagon contraction created parallel edges");
  AdjListTy reduced = edges_to_adj(reducedN, edges);
  for (int vertex : contracted)
    CHECK(reduced[vertex].size() == cycleLength, "contracted pentagon vertex has incorrect degree");

  std::vector<EdgeTy> protectedAttachments;
  for (int packed = 0; packed < static_cast<int>(selected.size()); packed++) {
    const Candidate& candidate = candidates[selected[packed]];
    int safeOffset = -1;
    for (int offset = 0; offset < cycleLength && safeOffset < 0; offset++) {
      bool clean = true;
      for (int delta = 0; delta < protectedSpokes; delta++)
        clean &= owner[candidate.outside[(offset + delta) % cycleLength]] == -1;
      if (clean)
        safeOffset = offset;
    }
    CHECK(safeOffset >= 0, "packed pentagon has no clean consecutive attachment pair");
    for (int delta = 0; delta < protectedSpokes; delta++) {
      const int outside = candidate.outside[(safeOffset + delta) % cycleLength];
      const EdgeTy edge(std::min(contracted[packed], remap[outside]),
                        std::max(contracted[packed], remap[outside]));
      CHECK(std::find(edges.begin(), edges.end(), edge) != edges.end());
      protectedAttachments.push_back(edge);
    }
  }

  return {std::move(reduced), std::move(protectedAttachments)};
}

/// Contract one chordless pentagon and construct the reduced-graph edge
/// requirements that preserve every possible parent edge.
PentagonFlexibilityContraction contractPentagonForFlexibility(const AdjListTy& adjList) {
  constexpr int cycleLength = 5;
  const int n = static_cast<int>(adjList.size());
  CHECK(n <= 63, "pentagon flexibility contraction uses a 64-bit vertex mask");
  for (const auto& neighbors : adjList)
    CHECK(neighbors.size() == 3, "pentagon flexibility contraction requires a cubic graph");

  std::vector<int> selected;
  std::vector<int> selectedOutside;
  uint64_t selectedMask = 0;
  std::vector<int> cycle;
  std::function<void(int, int)> findCycle = [&](const int start, const int current) {
    if (!selected.empty())
      return;
    if (static_cast<int>(cycle.size()) == cycleLength) {
      if (cycle[1] > cycle.back() ||
          std::find(adjList[current].begin(), adjList[current].end(), start) == adjList[current].end())
        return;
      uint64_t cycleMask = 0;
      for (int vertex : cycle)
        cycleMask |= uint64_t(1) << vertex;
      std::vector<int> outside;
      uint64_t outsideMask = 0;
      for (int i = 0; i < cycleLength; i++) {
        const int previous = cycle[(i + cycleLength - 1) % cycleLength];
        const int next = cycle[(i + 1) % cycleLength];
        int external = -1;
        for (int neighbor : adjList[cycle[i]]) {
          if (neighbor == previous || neighbor == next)
            continue;
          if ((cycleMask & (uint64_t(1) << neighbor)) != 0 || external != -1)
            return;
          external = neighbor;
        }
        if (external == -1)
          return;
        const uint64_t externalBit = uint64_t(1) << external;
        if ((outsideMask & externalBit) != 0)
          return;
        outsideMask |= externalBit;
        outside.push_back(external);
      }
      selected = cycle;
      selectedOutside = std::move(outside);
      selectedMask = cycleMask;
      return;
    }
    for (int next : adjList[current]) {
      if (next <= start || std::find(cycle.begin(), cycle.end(), next) != cycle.end())
        continue;
      cycle.push_back(next);
      findCycle(start, next);
      cycle.pop_back();
    }
  };
  for (int start = 0; start < n && selected.empty(); start++) {
    cycle = {start};
    findCycle(start, start);
  }
  if (selected.empty())
    return {};

  std::vector<int> remap(n, -1);
  int reducedN = 0;
  for (int vertex = 0; vertex < n; vertex++)
    if ((selectedMask & (uint64_t(1) << vertex)) == 0)
      remap[vertex] = reducedN++;
  const int contracted = reducedN++;
  auto mappedVertex = [&](const int vertex) {
    return (selectedMask & (uint64_t(1) << vertex)) != 0 ? contracted : remap[vertex];
  };
  std::vector<EdgeTy> parentEdges;
  std::vector<EdgeTy> reducedEdges;
  for (int first = 0; first < n; first++) {
    for (int second : adjList[first]) {
      if (first >= second)
        continue;
      parentEdges.push_back({first, second});
      const int mappedFirst = mappedVertex(first);
      const int mappedSecond = mappedVertex(second);
      if (mappedFirst != mappedSecond)
        reducedEdges.push_back(make_edge(mappedFirst, mappedSecond));
    }
  }
  const size_t edgesBeforeDeduplication = reducedEdges.size();
  sort_unique(reducedEdges);
  CHECK(reducedEdges.size() == edgesBeforeDeduplication,
        "pentagon flexibility contraction created parallel edges");
  AdjListTy reduced = edges_to_adj(reducedN, reducedEdges);
  CHECK(reduced[contracted].size() == cycleLength, "contracted pentagon vertex has incorrect degree");

  std::vector<EdgeTy> spokes;
  for (const int outside : selectedOutside)
    spokes.push_back(make_edge(contracted, remap[outside]));

  std::vector<std::vector<EdgeTy>> requirements;
  requirements.reserve(parentEdges.size());

  // A cycle target c_i is preserved by the checked local table when the two
  // endpoint spokes s_i,s_{i+1} are uncrossed.
  for (int position = 0; position < cycleLength; position++)
    requirements.push_back({spokes[position], spokes[(position + 1) % cycleLength]});

  // For an attachment target s_i, rotate the ordinary pentagon patch so its
  // protected consecutive spokes are s_{i+1},s_{i+2}.  The target spoke s_i
  // is then also uncrossed by that same patch.  Thus no separate
  // target-specific pentagon table is needed.
  for (int position = 0; position < cycleLength; position++)
    requirements.push_back({spokes[position], spokes[(position + 1) % cycleLength],
                            spokes[(position + 2) % cycleLength]});

  // Every other target edge is unchanged by the local patch.  Protect it
  // together with one fixed consecutive pair at the contracted pentagon.
  for (const EdgeTy& edge : parentEdges) {
    if ((selectedMask & (uint64_t(1) << edge.first)) != 0 ||
        (selectedMask & (uint64_t(1) << edge.second)) != 0)
      continue;
    requirements.push_back({spokes[0], spokes[1], make_edge(remap[edge.first], remap[edge.second])});
  }
  CHECK(requirements.size() == parentEdges.size(), "one parent-edge requirement was omitted");
  for (auto& requirement : requirements) {
    sort_unique(requirement);
    CHECK(2 <= requirement.size() && requirement.size() <= 3);
    for (const EdgeTy& edge : requirement)
      CHECK(std::binary_search(reducedEdges.begin(), reducedEdges.end(), edge),
            "pentagon flexibility requirement is not a reduced edge");
  }
  std::sort(requirements.begin(), requirements.end());

  return {std::move(reduced), std::move(requirements)};
}

/// Contract one chordless hexagon and return the four attachment edges that
/// the expansion certificate requires to remain uncrossed.
HexagonContraction contractHexagon(const AdjListTy& adjList) {
  constexpr int cycleLength = 6;
  const int n = static_cast<int>(adjList.size());
  CHECK(n <= 63, "hexagon contraction uses a 64-bit vertex mask");
  for (const auto& neighbors : adjList)
    CHECK(neighbors.size() == 3, "hexagon contraction is defined here only for cubic graphs");

  std::vector<int> selected;
  std::vector<int> selectedOutside;
  uint64_t selectedMask = 0;
  std::vector<int> cycle;
  std::function<void(int, int)> findCycle = [&](const int start, const int current) {
    if (!selected.empty())
      return;
    if (static_cast<int>(cycle.size()) == cycleLength) {
      if (cycle[1] > cycle.back() ||
          std::find(adjList[current].begin(), adjList[current].end(), start) == adjList[current].end())
        return;
      uint64_t cycleMask = 0;
      for (int vertex : cycle)
        cycleMask |= uint64_t(1) << vertex;
      std::vector<int> outside;
      uint64_t outsideMask = 0;
      for (int i = 0; i < cycleLength; i++) {
        const int previous = cycle[(i + cycleLength - 1) % cycleLength];
        const int next = cycle[(i + 1) % cycleLength];
        int external = -1;
        for (int neighbor : adjList[cycle[i]]) {
          if (neighbor == previous || neighbor == next)
            continue;
          if ((cycleMask & (uint64_t(1) << neighbor)) != 0 || external != -1)
            return;
          external = neighbor;
        }
        if (external == -1)
          return;
        const uint64_t externalBit = uint64_t(1) << external;
        if ((outsideMask & externalBit) != 0)
          return;
        outsideMask |= externalBit;
        outside.push_back(external);
      }
      selected = cycle;
      selectedOutside = std::move(outside);
      selectedMask = cycleMask;
      return;
    }
    for (int next : adjList[current]) {
      if (next <= start || std::find(cycle.begin(), cycle.end(), next) != cycle.end())
        continue;
      cycle.push_back(next);
      findCycle(start, next);
      cycle.pop_back();
    }
  };
  for (int start = 0; start < n && selected.empty(); start++) {
    cycle = {start};
    findCycle(start, start);
  }
  if (selected.empty())
    return {};

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
        edges.push_back({std::min(mappedU, mappedV), std::max(mappedU, mappedV)});
    }
  }
  const size_t edgesBeforeDeduplication = edges.size();
  sort_unique(edges);
  CHECK(edges.size() == edgesBeforeDeduplication, "hexagon contraction created parallel edges");
  AdjListTy reduced = edges_to_adj(reducedN, edges);
  CHECK(reduced[contracted].size() == cycleLength, "contracted hexagon vertex has incorrect degree");

  std::vector<EdgeTy> protectedAttachments;
  for (int position : {1, 2, 4, 5}) {
    const int outside = selectedOutside[position];
    protectedAttachments.push_back({std::min(contracted, remap[outside]), std::max(contracted, remap[outside])});
  }
  return {std::move(reduced), std::move(protectedAttachments)};
}
/// Reconstruct the drawing represented by a satisfying assignment and check
/// independently that its proposed crossings produce a planarization.
Result checkedDrawing(const Params& params, const InputGraph& graph, SATModel& model, Solver& solver) {
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

  const Result drawing = checkedDrawing(params, graph, model, solver);
  for (int edge : uncrossedEdges)
    CHECK(!drawing.isCrossed[edge], "required edge %d is crossed", edge);
}

/// Return the sorted set of original edges crossed in a checked drawing.
std::vector<int> crossedEdges(const Result& drawing) {
  std::vector<int> crossed;
  for (const auto& [first, second] : drawing.crossings) {
    crossed.push_back(first);
    crossed.push_back(second);
  }
  sort_unique(crossed);
  return crossed;
}

/// Test whether a drawing leaves every edge in a requirement uncrossed.
bool covers(const std::vector<int>& requirement, const std::vector<int>& crossed) {
  return std::none_of(requirement.begin(), requirement.end(), [&](const int edge) {
    return std::binary_search(crossed.begin(), crossed.end(), edge);
  });
}

/// Accumulate checked drawings until they cover every required edge set.
void verifyFlexibility(
    const Params& params, const InputGraph& graph, const int flexibility,
    const std::vector<std::vector<int>>& explicitRequirements = {}) {
  CHECK((flexibility > 0) != !explicitRequirements.empty(),
        "specify either a flexibility value or explicit requirements");
  for (const auto& requirement : explicitRequirements) {
    CHECK(!requirement.empty() && all_unique(requirement), "invalid uncrossed-edge requirement");
    for (int edge : requirement) {
      CHECK(0 <= edge && edge < static_cast<int>(graph.edges.size()), "uncrossed-edge requirement is out of range");
    }
  }

  initCrossablePairs(params, graph);
  SATModel model;
  encodeStackPlanar(model, graph, params);
  Solver solver;
  initSATSolver(params, graph, model, solver);
  CHECK(solveSATModel(params, model, solver) == l_True, "failed to find an initial 1-planar drawing");

  Result drawing = checkedDrawing(params, graph, model, solver);
  std::vector<std::vector<int>> witnesses = {crossedEdges(drawing)};

  auto findUncoveredRequirement = [&]() {
    if (!explicitRequirements.empty()) {
      for (const auto& requirement : explicitRequirements) {
        const bool isCovered = std::any_of(witnesses.begin(), witnesses.end(), [&](const std::vector<int>& crossed) {
          return covers(requirement, crossed);
        });
        if (!isCovered)
          return requirement;
      }
      return std::vector<int>();
    }

    std::vector<int> candidate;
    std::function<bool(int)> findHittingSet = [&](const int limit) {
      const std::vector<int>* unhitWitness = nullptr;
      for (const auto& crossed : witnesses) {
        const bool hit = std::any_of(candidate.begin(), candidate.end(), [&](const int edge) {
          return std::binary_search(crossed.begin(), crossed.end(), edge);
        });
        if (!hit) {
          unhitWitness = &crossed;
          break;
        }
      }
      if (unhitWitness == nullptr)
        return true;
      if (static_cast<int>(candidate.size()) == limit || unhitWitness->empty())
        return false;
      for (int edge : *unhitWitness) {
        candidate.push_back(edge);
        if (findHittingSet(limit))
          return true;
        candidate.pop_back();
      }
      return false;
    };

    for (int limit = 1; limit <= flexibility; limit++) {
      candidate.clear();
      if (findHittingSet(limit))
        return candidate;
    }
    return std::vector<int>();
  };

  while (true) {
    const std::vector<int> requirement = findUncoveredRequirement();
    if (requirement.empty())
      return;

    vec<Lit> assumptions;
    for (int edge : requirement)
      assumptions.push(model.getSolverLit(model.getCross1Var(graph.n + edge, false)));
    CHECK(solver.solveLimited(assumptions) == l_True, "failed uncrossed-edge requirement");

    drawing = checkedDrawing(params, graph, model, solver);
    for (int edge : requirement)
      CHECK(!drawing.isCrossed[edge], "required edge %d is crossed", edge);
    witnesses.push_back(crossedEdges(drawing));
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

/// Convert endpoint-based requirements into edge-index requirements.
std::vector<std::vector<int>> edgeRequirements(
    const InputGraph& graph, const std::vector<std::vector<EdgeTy>>& requirements) {
  std::vector<std::vector<int>> result;
  for (const auto& requirement : requirements)
    result.push_back(edgeIndices(graph, requirement));
  return result;
}

/// Check that an input record belongs to the census class for its order.
void validateInput(const int n, const std::vector<EdgeTy>& edges) {
  CHECK(4 <= n && n <= 28 && n % 2 == 0, "verification expects an even order from 4 through 28");
  CHECK(static_cast<int>(edges.size()) * 2 == 3 * n &&
            minDegree(n, edges) == 3 && maxDegree(n, edges) == 3,
        "verification expects a cubic graph");
  CHECK(isConnected(n, edges), "verification censuses contain connected graphs");
  if (n >= 24) {
    CHECK(computeGirth(n, edges) >= 5, "orders 24, 26, and 28 require girth at least five");
  }
}

/// Verify three-edge flexibility for a nonplanar input through order 22.
void verifyCubic22(const int n, const std::vector<EdgeTy>& edges, const Params& params) {
  CHECK(n <= 22);
  InputGraph graph(n, edges, {});
  verifyFlexibility(params, graph, 3);
}

/// Verify two-edge flexibility for a nonplanar order-24 input.
void verifyCubic24(const int n, const std::vector<EdgeTy>& edges, const Params& params) {
  CHECK(n == 24);
  InputGraph graph(n, edges, {});
  verifyFlexibility(params, graph, 2);
}

/// Verify one-edge flexibility for a nonplanar order-26 input; return true
/// when the proof uses a contracted pentagon and false for a direct check.
bool verifyCubic26(
    const int n, const std::vector<EdgeTy>& edges, const AdjListTy& adj, const Params& params) {
  CHECK(n == 26);
  const PentagonFlexibilityContraction contraction = contractPentagonForFlexibility(adj);
  if (contraction.adjList.empty()) {
    InputGraph graph(n, edges, {});
    verifyFlexibility(params, graph, 1);
    return false;
  }

  InputGraph graph(static_cast<int>(contraction.adjList.size()), adj_to_edges(contraction.adjList), {});
  CHECK(contraction.requirements.size() == edges.size(), "pentagon contraction omitted a parent edge");
  verifyFlexibility(params, graph, 0, edgeRequirements(graph, contraction.requirements));
  return true;
}

/// Certificate selected for a nonplanar order-28 input.
enum class Cubic28Certificate { PENTAGON, HEXAGON, DIRECT };

/// Verify a nonplanar order-28 input by pentagon expansion, hexagon expansion,
/// or a direct drawing, and return the certificate that was used.
Cubic28Certificate verifyCubic28(
    const int n, const std::vector<EdgeTy>& edges, const AdjListTy& adj, const Params& params) {
  CHECK(n == 28);
  const PackedPentagonContraction pentagons = contractPackedPentagons(adj);
  if (!pentagons.adjList.empty()) {
    InputGraph graph(static_cast<int>(pentagons.adjList.size()), adj_to_edges(pentagons.adjList), {});
    const size_t removed = adj.size() - pentagons.adjList.size();
    CHECK(removed % 4 == 0 &&
              pentagons.protectedAttachments.size() == removed / 2,
          "invalid packed-pentagon contraction");
    if (!isPlanar(graph.n, graph.edges, 0)) {
      verifyDrawing(params, graph, edgeIndices(graph, pentagons.protectedAttachments));
    }
    return Cubic28Certificate::PENTAGON;
  }

  const HexagonContraction hexagon = contractHexagon(adj);
  if (!hexagon.adjList.empty()) {
    InputGraph graph(static_cast<int>(hexagon.adjList.size()), adj_to_edges(hexagon.adjList), {});
    CHECK(hexagon.protectedAttachments.size() == 4, "invalid hexagon contraction");
    if (!isPlanar(graph.n, graph.edges, 0)) {
      verifyDrawing(params, graph, edgeIndices(graph, hexagon.protectedAttachments));
    }
    return Cubic28Certificate::HEXAGON;
  }

  InputGraph graph(n, edges, {});
  verifyDrawing(params, graph);
  return Cubic28Certificate::DIRECT;
}

}  // namespace

/// Verify the unique cubic-28 proof obligation selected by each graph's order.
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
  int numPentagonFlexible = 0;
  int numDirect1Flexible = 0;
  int numPentagonExpanded = 0;
  int numHexagonExpanded = 0;
  int numDirect28 = 0;

  for (int i = 0; i < numGraphs; i++) {
    const auto& adj = graphs->next().second;
    const int n = static_cast<int>(adj.size());
    const std::vector<EdgeTy> edges = adj_to_edges(adj);
    validateInput(n, edges);

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
      if (verifyCubic26(n, edges, adj, params))
        numPentagonFlexible++;
      else
        numDirect1Flexible++;
      num1Planar++;
    } else {
      const Cubic28Certificate certificate = verifyCubic28(n, edges, adj, params);
      if (certificate == Cubic28Certificate::PENTAGON)
        numPentagonExpanded++;
      else if (certificate == Cubic28Certificate::HEXAGON)
        numHexagonExpanded++;
      else
        numDirect28++;
      num1Planar++;
    }

    LOG_EVERY_MS(30000, "verified %'d of %'d graphs", i + 1, numGraphs);
  }

  LOG("processed %'d graphs in %s", numGraphs, ms_to_str(start, std::chrono::steady_clock::now()).c_str());
  LOG("#planar = %'d; #1-planar = %'d", numPlanar, num1Planar);
  LOG("#3-flexible = %'d; #2-flexible = %'d", num3Flexible, num2Flexible);
  LOG("#pentagon-flexible = %'d; #direct-1-flexible = %'d", numPentagonFlexible, numDirect1Flexible);
  LOG("#pentagon-expanded = %'d; #hexagon-expanded = %'d; #direct-28 = %'d",
      numPentagonExpanded, numHexagonExpanded, numDirect28);
}
