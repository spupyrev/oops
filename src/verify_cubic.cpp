#include "cmd_options.h"
#include "graph_algorithms.h"
#include "io.h"
#include "logging.h"
#include "one_planar.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <climits>
#include <cstdlib>
#include <functional>
#include <memory>
#include <vector>

std::unique_ptr<GraphList> genGraphs(CMDOptions& options);

namespace {

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
    const std::vector<EdgeTy> reducedEdges = adj_to_edges(reduced);
    if (computeGirth(reducedN, reducedEdges) < 5)
      return false;

    // A planar reduced graph has the required three-edge star uncrossed.  A
    // nonplanar reduced graph is covered by Claim 3 only when biconnected.
    return is2Connected(reducedN, reducedEdges) || isPlanar(reducedN, reducedEdges, 0);
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

/// An uncrossed-edge condition consisting of fixed edges and a lower bound
/// on how many candidate edges must also be uncrossed.
struct CoverageRequirement {
  std::vector<int> fixed;
  std::vector<int> candidates;
  int minimumUncrossed = 0;
  std::vector<std::pair<int, int>> forbiddenCrossedPairs;
};

/// Test whether a drawing satisfies one coverage requirement.
bool isCovered(
    const CoverageRequirement& requirement, const std::vector<int>& crossed) {
  if (!areDisjoint(requirement.fixed, crossed))
    return false;
  const int crossedCandidates = static_cast<int>(std::count_if(
      requirement.candidates.begin(), requirement.candidates.end(),
      [&](const int edge) {
        return std::binary_search(crossed.begin(), crossed.end(), edge);
      }));
  if (crossedCandidates >
      static_cast<int>(requirement.candidates.size()) -
          requirement.minimumUncrossed)
    return false;
  return std::none_of(
      requirement.forbiddenCrossedPairs.begin(),
      requirement.forbiddenCrossedPairs.end(),
      [&](const auto& pair) {
        return std::binary_search(crossed.begin(), crossed.end(), pair.first) &&
               std::binary_search(crossed.begin(), crossed.end(), pair.second);
      });
}

/// Accumulate verified drawings until they prove k-flexibility or cover all
/// prescribed sets of edges that must be simultaneously uncrossed.
void verifyFlexibility(
    const Params& params, const InputGraph& graph, const int k,
    const std::vector<std::vector<int>>& prescribedSets = {},
    const bool minimizeWitnesses = false,
    const std::vector<CoverageRequirement>& coverageRequirements = {}) {
  CHECK(
      static_cast<int>(k > 0) + static_cast<int>(!prescribedSets.empty()) +
              static_cast<int>(!coverageRequirements.empty()) ==
          1,
      "specify exactly one flexibility obligation");

  initCrossablePairs(params, graph);
  SATModel model;
  encodeStackPlanar(model, graph, params);
  Solver solver;
  initSATSolver(params, graph, model, solver);

  std::vector<Lit> activations;
  for (const CoverageRequirement& requirement : coverageRequirements) {
    CHECK(
        0 <= requirement.minimumUncrossed &&
            requirement.minimumUncrossed <=
                static_cast<int>(requirement.candidates.size()),
        "invalid uncrossed-edge lower bound");
    const Lit activation = mkLit(solver.newVar(true, false));
    activations.push_back(activation);
    const int forbiddenCrossed =
        static_cast<int>(requirement.candidates.size()) -
        requirement.minimumUncrossed + 1;
    std::vector<int> selected;
    std::function<void(int)> addClauses = [&](const int next) {
      if (static_cast<int>(selected.size()) == forbiddenCrossed) {
        vec<Lit> clause;
        clause.push(~activation);
        for (int candidate : selected)
          clause.push(~model.getSolverLit(
              model.getCross1Var(graph.n + candidate, true)));
        CHECK(solver.addClause(clause), "failed to add a guarded coverage clause");
        return;
      }
      for (int candidate = next;
           candidate < static_cast<int>(requirement.candidates.size());
           candidate++) {
        selected.push_back(requirement.candidates[candidate]);
        addClauses(candidate + 1);
        selected.pop_back();
      }
    };
    addClauses(0);
    for (const auto& [first, second] : requirement.forbiddenCrossedPairs) {
      vec<Lit> clause;
      clause.push(~activation);
      clause.push(~model.getSolverLit(
          model.getCross1Var(graph.n + first, true)));
      clause.push(~model.getSolverLit(
          model.getCross1Var(graph.n + second, true)));
      CHECK(solver.addClause(clause), "failed to add a guarded crossing-pair clause");
    }
  }

  if (coverageRequirements.size() == 1) {
    vec<Lit> assumptions;
    assumptions.push(activations.front());
    for (int edge : coverageRequirements.front().fixed)
      assumptions.push(
          model.getSolverLit(model.getCross1Var(graph.n + edge, false)));
    CHECK(solver.solveLimited(assumptions) == l_True,
          "failed uncrossed-edge coverage requirement");
    const Result drawing = verifiedDrawing(params, graph, model, solver);
    CHECK(
        isCovered(coverageRequirements.front(), crossedEdges(drawing)),
        "drawing does not satisfy its coverage requirement");
    return;
  }

  if (activations.empty()) {
    CHECK(solveSATModel(params, model, solver) == l_True,
          "failed to find an initial 1-planar drawing");
  } else {
    vec<Lit> assumptions;
    for (Lit activation : activations)
      assumptions.push(~activation);
    CHECK(solver.solveLimited(assumptions) == l_True,
          "failed to find an initial 1-planar drawing");
  }

  Result drawing = verifiedDrawing(params, graph, model, solver);
  if (minimizeWitnesses)
    minimizeCrossings(graph, drawing, 0);
  std::vector<std::vector<int>> crossedEdgeSets = {crossedEdges(drawing)};

  int selectedCoverage = -1;
  auto findUncoveredEdgeSet = [&]() {
    selectedCoverage = -1;
    if (!coverageRequirements.empty()) {
      int bestScore = INT_MAX;
      for (int index = 0;
           index < static_cast<int>(coverageRequirements.size()); index++) {
        const CoverageRequirement& requirement = coverageRequirements[index];
        const bool isCovered = std::any_of(
            crossedEdgeSets.begin(), crossedEdgeSets.end(),
            [&](const std::vector<int>& crossed) {
              return ::isCovered(requirement, crossed);
            });
        if (isCovered)
          continue;
        int score = 0;
        for (const auto& crossed : crossedEdgeSets) {
          for (int edge : requirement.fixed)
            score += std::binary_search(crossed.begin(), crossed.end(), edge);
          const int crossedCandidates = static_cast<int>(std::count_if(
              requirement.candidates.begin(), requirement.candidates.end(),
              [&](const int edge) {
                return std::binary_search(crossed.begin(), crossed.end(), edge);
              }));
          score += std::max(
              0, crossedCandidates -
                     (static_cast<int>(requirement.candidates.size()) -
                      requirement.minimumUncrossed));
        }
        if (score < bestScore) {
          selectedCoverage = index;
          bestScore = score;
        }
      }
      return selectedCoverage == -1
                 ? std::vector<int>()
                 : coverageRequirements[selectedCoverage].fixed;
    }

    if (!prescribedSets.empty()) {
      // Prefer the requirement that conflicts least with the drawings already
      // found; this avoids difficult assumption queries without changing the
      // set of requirements that must be covered.
      const std::vector<int>* best = nullptr;
      int bestScore = INT_MAX;
      for (const auto& prescribed : prescribedSets) {
        const bool isCovered =
            std::any_of(crossedEdgeSets.begin(), crossedEdgeSets.end(),
                        [&](const std::vector<int>& crossed) { return areDisjoint(prescribed, crossed); });
        if (isCovered)
          continue;
        int score = 0;
        for (const auto& crossed : crossedEdgeSets)
          for (int edge : prescribed)
            score += std::binary_search(crossed.begin(), crossed.end(), edge);
        if (score < bestScore) {
          best = &prescribed;
          bestScore = score;
        }
      }
      return best == nullptr ? std::vector<int>() : *best;
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
    for (int index = 0; index < static_cast<int>(activations.size()); index++)
      assumptions.push(
          index == selectedCoverage ? activations[index] : ~activations[index]);
    for (int edge : prescribed)
      assumptions.push(model.getSolverLit(model.getCross1Var(graph.n + edge, false)));
    CHECK(solver.solveLimited(assumptions) == l_True, "failed uncrossed-edge requirement");

    drawing = verifiedDrawing(params, graph, model, solver);
    if (minimizeWitnesses)
      minimizeCrossings(graph, drawing, 0);
    for (int edge : prescribed)
      CHECK(!drawing.isCrossed[edge], "required edge %d is crossed", edge);
    const std::vector<int> crossed = crossedEdges(drawing);
    if (selectedCoverage != -1)
      CHECK(isCovered(coverageRequirements[selectedCoverage], crossed),
            "drawing does not satisfy its coverage requirement");
    crossedEdgeSets.push_back(crossed);
  }
}

enum class VerificationInput {
  CUBIC,
  FIVE_CYCLE_CORE,
  SIX_CYCLE_CORE,
  STAR_MARKER
};

/// Check that an input record is cubic or is one of the two contracted graphs
/// used for Claim 3.
VerificationInput validateInput(const AdjListTy& adjList, const std::vector<EdgeTy>& edges) {
  const int n = static_cast<int>(adjList.size());
  CHECK(isConnected(n, edges), "verification inputs contain connected graphs");

  if (n == 28) {
    int degreeOne = -1;
    int degreeTwo = -1;
    int degreeFour = -1;
    bool markerDegrees = true;
    for (int vertex = 0; vertex < n; vertex++) {
      const int degree = static_cast<int>(adjList[vertex].size());
      int* location = nullptr;
      if (degree == 1)
        location = &degreeOne;
      else if (degree == 2)
        location = &degreeTwo;
      else if (degree == 4)
        location = &degreeFour;
      else if (degree != 3)
        markerDegrees = false;
      if (location != nullptr) {
        if (*location != -1)
          markerDegrees = false;
        *location = vertex;
      }
    }
    if (markerDegrees && degreeOne != -1 && degreeTwo != -1 &&
        degreeFour != -1 && contains(adjList[degreeOne], degreeTwo) &&
        contains(adjList[degreeTwo], degreeFour))
      return VerificationInput::STAR_MARKER;
  }

  int exceptionalVertex = -1;
  int exceptionalDegree = -1;
  for (int vertex = 0; vertex < n; vertex++) {
    const int degree = static_cast<int>(adjList[vertex].size());
    if (degree == 3)
      continue;
    CHECK(exceptionalVertex == -1, "verification permits at most one noncubic vertex");
    exceptionalVertex = vertex;
    exceptionalDegree = degree;
  }
  if (exceptionalVertex != -1) {
    if (n == 22 && exceptionalDegree == 5)
      return VerificationInput::FIVE_CYCLE_CORE;
    CHECK(n == 21 && exceptionalDegree == 6, "invalid contracted graph for Claim 3");
    return VerificationInput::SIX_CYCLE_CORE;
  }

  CHECK(4 <= n && n <= 28 && n % 2 == 0, "verification expects an even order from 4 through 28");
  if (n >= 24)
    CHECK(computeGirth(n, edges) >= 5, "orders 24, 26, and 28 require girth at least five");
  return VerificationInput::CUBIC;
}

/// Return the edge indices incident with one vertex.
std::vector<int> incidentEdges(const InputGraph& graph, const int vertex) {
  std::vector<int> result;
  for (int edge = 0; edge < static_cast<int>(graph.edges.size()); edge++)
    if (graph.edges[edge].first == vertex || graph.edges[edge].second == vertex)
      result.push_back(edge);
  return result;
}

/// Remove a pendant two-edge marker and verify that the marked vertex has an
/// uncrossed three-edge star in the underlying order-26 graph.
void verifyMarkedStar(const AdjListTy& adjList, const Params& params) {
  const int n = static_cast<int>(adjList.size());
  int target = -1;
  std::vector<int> remap(n, -1);
  int baseN = 0;
  for (int vertex = 0; vertex < n; vertex++) {
    const int degree = static_cast<int>(adjList[vertex].size());
    if (degree == 4)
      target = vertex;
    if (degree >= 3)
      remap[vertex] = baseN++;
  }
  CHECK(baseN == 26 && target != -1, "invalid marked-star graph");

  std::vector<EdgeTy> baseEdges;
  for (int first = 0; first < n; first++)
    for (int second : adjList[first])
      if (first < second && remap[first] != -1 && remap[second] != -1)
        baseEdges.push_back(make_edge(remap[first], remap[second]));
  InputGraph graph(baseN, baseEdges, {});
  CHECK(graph.edges.size() == 39, "marked-star base is not cubic");
  const std::vector<int> star = incidentEdges(graph, remap[target]);
  CHECK(star.size() == 3, "marked vertex does not have a three-edge star");
  verifyDrawing(params, graph, star);
}

/// Verify the cyclic-order-independent conditions for expanding a contracted
/// 5-cycle while preserving every three-edge star of the original graph.
void verifyFiveCycleCore(
    const std::vector<EdgeTy>& edges, const AdjListTy& adjList, const Params& params) {
  const int hub = static_cast<int>(
      std::find_if(adjList.begin(), adjList.end(),
                   [](const auto& neighbors) { return neighbors.size() == 5; }) -
      adjList.begin());
  CHECK(hub < static_cast<int>(adjList.size()), "5-cycle core has no degree-5 vertex");

  InputGraph graph(static_cast<int>(adjList.size()), edges, {});
  const std::vector<int> attachments = incidentEdges(graph, hub);
  CHECK(attachments.size() == 5, "5-cycle core has an invalid attachment set");

  std::vector<CoverageRequirement> requirements = {
      CoverageRequirement{attachments, {}, 0, {}}};
  for (int vertex = 0; vertex < graph.n; vertex++) {
    if (vertex == hub)
      continue;
    const std::vector<int> star = incidentEdges(graph, vertex);
    CHECK(star.size() == 3, "non-hub vertex of a 5-cycle core is not cubic");
    std::vector<int> eligible;
    for (int attachment : attachments) {
      if (!contains(star, attachment))
        eligible.push_back(attachment);
    }
    CHECK(eligible.size() == 4 || eligible.size() == 5,
          "invalid attachment incidence in a 5-cycle core");
    requirements.push_back(
        CoverageRequirement{star, std::move(eligible), 3, {}});
  }
  CHECK(requirements.size() == 22,
        "unexpected number of 5-cycle-core requirements");
  verifyFlexibility(params, graph, 0, {}, true, requirements);
}

/// Verify the conditions for expanding one member of a pair of
/// vertex-disjoint 6-cycles in the Claim 3 proof.  Vertices 0,...,5 are the
/// attachment endpoints in their cyclic order.
void verifySixCycleCore(
    const std::vector<EdgeTy>& edges, const AdjListTy& adjList, const Params& params) {
  const int hub = static_cast<int>(
      std::find_if(adjList.begin(), adjList.end(),
                   [](const auto& neighbors) { return neighbors.size() == 6; }) -
      adjList.begin());
  CHECK(hub < static_cast<int>(adjList.size()), "6-cycle core has no degree-6 vertex");

  InputGraph graph(static_cast<int>(adjList.size()), edges, {});
  const std::vector<int> attachments = incidentEdges(graph, hub);
  CHECK(attachments.size() == 6, "6-cycle core has an invalid attachment set");
  for (int position = 0; position < 6; position++)
    CHECK(graph.findEdgeIndex(position, hub) == attachments[position],
          "6-cycle attachments are not labeled in cyclic order");

  std::vector<CoverageRequirement> requirements;
  const char* diagnosticVertex = std::getenv("OOPS_CUBIC_VERTEX");
  for (int vertex = 0; vertex < graph.n; vertex++) {
    if (vertex == hub)
      continue;
    if (vertex < 6)
      continue;
    if (diagnosticVertex != nullptr && vertex != std::atoi(diagnosticVertex))
      continue;
    const std::vector<int> star = incidentEdges(graph, vertex);
    CHECK(star.size() == 3, "non-hub vertex of a 6-cycle core is not cubic");
    std::vector<std::pair<int, int>> forbiddenPairs;
    for (int first = 0; first < 6; first++)
      for (int second = first + 1; second < 6; second++)
        if (second != first + 3)
          forbiddenPairs.emplace_back(
              attachments[first], attachments[second]);
    requirements.push_back(CoverageRequirement{
        star, {}, 0, std::move(forbiddenPairs)});
  }
  CHECK(diagnosticVertex != nullptr || requirements.size() == 14,
        "unexpected number of 6-cycle-core requirements");
  for (const CoverageRequirement& requirement : requirements)
    verifyFlexibility(params, graph, 0, {}, false, {requirement});
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
  verifyFlexibility(params, graph, 0, incidentEdges, true);
}

/// Method used to verify a nonplanar order-28 input.
enum class Cubic28Method {
  FIVE_CYCLE_EXPANSION,
  DIRECT_GIRTH_FIVE
};

/// Verify Claim 4 by a 5-cycle expansion or a direct drawing.
Cubic28Method verifyCubic28(
    const int n, const std::vector<EdgeTy>& edges, const AdjListTy& adj, const Params& params) {
  CHECK(n == 28);
  const FiveCycleReduction fiveCycleReduction = findFiveCycleReduction(adj);
  if (fiveCycleReduction == FiveCycleReduction::SUCCESS)
    return Cubic28Method::FIVE_CYCLE_EXPANSION;
  CHECK(fiveCycleReduction == FiveCycleReduction::FAILURE,
        "Claim 4 verification requires girth exactly five");
  InputGraph graph(n, edges, {});
  verifyDrawing(params, graph);
  return Cubic28Method::DIRECT_GIRTH_FIVE;
}

}  // namespace

/// Verify the claim or claims prescribed by the order of each input graph.
void verifyCubic(CMDOptions& options) {
  Params params;
  params.cubicVerification = true;
  params.useSATConstraints = true;
  params.useUNSATConstraints = true;
  params.crossPriority = false;

  auto graphs = genGraphs(options);
  const int numGraphs = graphs->size();
  const auto start = std::chrono::steady_clock::now();

  int numPlanar = 0;
  int num1Planar = 0;
  int num3Flexible = 0;
  int num2Flexible = 0;
  int numFiveCycleCores = 0;
  int numSixCycleCores = 0;
  int numMarkedStars = 0;
  int numDirectStarVerified = 0;
  int numFiveCycleExpanded = 0;
  int numFiveCycleDirect = 0;

  for (int i = 0; i < numGraphs; i++) {
    if (std::getenv("OOPS_CUBIC_TRACE") != nullptr)
      LOG("starting graph %d", i);
    const auto& adj = graphs->next().second;
    const int n = static_cast<int>(adj.size());
    const std::vector<EdgeTy> edges = adj_to_edges(adj);
    const VerificationInput input = validateInput(adj, edges);
    if (input == VerificationInput::FIVE_CYCLE_CORE)
      numFiveCycleCores++;
    else if (input == VerificationInput::SIX_CYCLE_CORE)
      numSixCycleCores++;
    else if (input == VerificationInput::STAR_MARKER)
      numMarkedStars++;

    if (isPlanar(n, edges, 0)) {
      numPlanar++;
    } else if (input == VerificationInput::FIVE_CYCLE_CORE) {
      verifyFiveCycleCore(edges, adj, params);
      num1Planar++;
    } else if (input == VerificationInput::SIX_CYCLE_CORE) {
      verifySixCycleCore(edges, adj, params);
      num1Planar++;
    } else if (input == VerificationInput::STAR_MARKER) {
      verifyMarkedStar(adj, params);
      num1Planar++;
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
      numDirectStarVerified++;
      num1Planar++;
    } else {
      const Cubic28Method method = verifyCubic28(n, edges, adj, params);
      if (method == Cubic28Method::FIVE_CYCLE_EXPANSION)
        numFiveCycleExpanded++;
      else
        numFiveCycleDirect++;
      num1Planar++;
    }

    LOG_EVERY_MS(30000, "verified %'d of %'d graphs", i + 1, numGraphs);
  }

  LOG("processed %'d graphs in %s", numGraphs, ms_to_str(start, std::chrono::steady_clock::now()).c_str());
  LOG("#planar = %'d; #1-planar = %'d", numPlanar, num1Planar);
  LOG("#3-flexible = %'d; #2-flexible = %'d", num3Flexible, num2Flexible);
  LOG("#5-cycle-cores = %'d; #6-cycle-cores = %'d; #marked-stars = %'d; "
      "#direct-star-verified = %'d",
      numFiveCycleCores, numSixCycleCores, numMarkedStars,
      numDirectStarVerified);
  LOG("#5-cycle-expanded = %'d; #direct-girth5 = %'d",
      numFiveCycleExpanded, numFiveCycleDirect);
}
