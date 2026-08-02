#include "cmd_options.h"
#include "graph_algorithms.h"
#include "io.h"
#include "logging.h"
#include "one_planar.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <fstream>
#include <functional>
#include <memory>
#include <set>
#include <string>
#include <vector>

std::unique_ptr<GraphList> genGraphs(CMDOptions& options);

namespace {

struct VerificationTimedOut {
  std::string stage;
  uint64_t conflicts;
  uint64_t decisions;
};

void requireSatisfiable(
    Solver& solver, const lbool result, const std::string& failure) {
  if (result == l_Undef)
    throw VerificationTimedOut{
        failure, solver.conflicts, solver.decisions};
  CHECK(result == l_True, "%s", failure.c_str());
}

/// Encode an order-at-most-62 adjacency list as one graph6 record.
std::string graph6(const AdjListTy& adjList) {
  const int n = static_cast<int>(adjList.size());
  CHECK(1 <= n && n <= 62, "cubic residue requires small graph6 order");
  std::string result(1, static_cast<char>(n + 63));
  int value = 0;
  int bits = 0;
  for (int second = 1; second < n; second++) {
    for (int first = 0; first < second; first++) {
      value = (value << 1) | contains(adjList[first], second);
      if (++bits == 6) {
        result.push_back(static_cast<char>(value + 63));
        value = 0;
        bits = 0;
      }
    }
  }
  if (bits != 0) {
    value <<= 6 - bits;
    result.push_back(static_cast<char>(value + 63));
  }
  return result;
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
  requireSatisfiable(
      solver,
      solver.solveLimited(assumptions),
      "failed to find a required 1-planar drawing");

  const Result drawing = verifiedDrawing(params, graph, model, solver);
  for (int edge : uncrossedEdges)
    CHECK(!drawing.isCrossed[edge], "required edge %d is crossed", edge);
}

/// Remove crossing candidates incident with edges known to be uncrossed
/// before encoding, then verify the same condition in the final drawing.
void verifyDrawingWithFixedEdges(
    Params& params, const InputGraph& graph,
    const std::vector<int>& uncrossedEdges) {
  CHECK(params.fixCross1.empty(),
        "marked verification expects no fixed crossing");
  struct ClearFixedEdges {
    Params& params;
    ~ClearFixedEdges() { params.fixCross1.clear(); }
  } clearFixedEdges{params};
  for (int edge : uncrossedEdges) {
    if (!params.fixCross1.empty())
      params.fixCross1 += ';';
    const auto& [first, second] = graph.edges[edge];
    params.fixCross1 += std::to_string(first) + ':' +
        std::to_string(second) + '-';
  }
  verifyDrawing(params, graph, uncrossedEdges);
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
  return true;
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
  }

  if (coverageRequirements.size() == 1) {
    vec<Lit> assumptions;
    assumptions.push(activations.front());
    for (int edge : coverageRequirements.front().fixed)
      assumptions.push(
          model.getSolverLit(model.getCross1Var(graph.n + edge, false)));
    requireSatisfiable(
        solver,
        solver.solveLimited(assumptions),
        "failed uncrossed-edge coverage requirement");
    const Result drawing = verifiedDrawing(params, graph, model, solver);
    CHECK(
        isCovered(coverageRequirements.front(), crossedEdges(drawing)),
        "drawing does not satisfy its coverage requirement");
    return;
  }

  if (activations.empty()) {
    requireSatisfiable(
        solver,
        solveSATModel(params, model, solver),
        "failed to find an initial 1-planar drawing");
  } else {
    vec<Lit> assumptions;
    for (Lit activation : activations)
      assumptions.push(~activation);
    requireSatisfiable(
        solver,
        solver.solveLimited(assumptions),
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
    const std::string stage =
        selectedCoverage == -1
            ? "flexibility requirement"
            : "coverage requirement " + std::to_string(selectedCoverage);
    requireSatisfiable(
        solver,
        solver.solveLimited(assumptions),
        stage);

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
  FIVE_CYCLE_PATH_STAR,
  SIX_HUB,
  SEVEN_HUB
};

/// Check that an input record is cubic or is a marked reduced graph.
VerificationInput validateInput(const AdjListTy& adjList, const std::vector<EdgeTy>& edges) {
  const int n = static_cast<int>(adjList.size());
  CHECK(isConnected(n, edges), "verification inputs contain connected graphs");

  if (n == 27) {
    int marker = -1;
    for (int vertex = 0; vertex < n; vertex++) {
      const int degree = static_cast<int>(adjList[vertex].size());
      if (degree == 1) {
        CHECK(marker == -1, "marked path reduction has two markers");
        marker = vertex;
      }
    }
    CHECK(marker != -1, "marked path reduction has no marker");
    const int target = adjList[marker].front();
    CHECK(adjList[target].size() == 4,
          "marked path reduction has an invalid target");
    for (int vertex = 0; vertex < n; vertex++)
      CHECK(vertex == marker || vertex == target || adjList[vertex].size() == 3,
            "marked path reduction has an invalid degree");
    return VerificationInput::FIVE_CYCLE_PATH_STAR;
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
    if (n == 21 && exceptionalDegree == 6)
      return VerificationInput::SIX_HUB;
    if (n == 20 && exceptionalDegree == 7)
      return VerificationInput::SEVEN_HUB;
    CHECK(n == 22 && exceptionalDegree == 5,
          "invalid contracted graph for Claim 3");
    return VerificationInput::FIVE_CYCLE_CORE;
  }

  CHECK(4 <= n && n <= 26 && n % 2 == 0,
        "verification expects an even order from 4 through 26");
  if (n >= 24)
    CHECK(computeGirth(n, edges) >= 5,
          "orders 24 and 26 require girth at least five");
  return VerificationInput::CUBIC;
}

std::string normalizedRotation(std::vector<int> rotation) {
  const auto zero = std::find(rotation.begin(), rotation.end(), 0);
  std::rotate(rotation.begin(), zero, rotation.end());
  std::vector<int> reverse = {0};
  reverse.insert(reverse.end(), rotation.rbegin(), rotation.rend() - 1);
  const std::vector<int>& selected = std::min(rotation, reverse);
  std::string result;
  for (int label : selected)
    result.push_back(static_cast<char>('0' + label));
  return result;
}

std::string hubRotation(
    const std::vector<int>& positions, const int hub,
    const std::vector<int>& divisions) {
  // In the two-page embedding, the incident arcs on each page are nested.
  // Reading around the hub therefore gives the left division vertices from
  // farthest to nearest, followed by the right division vertices from
  // farthest to nearest.  Reversing the whole list represents reflection.
  std::vector<std::pair<int, int>> left;
  std::vector<std::pair<int, int>> right;
  for (int label = 0; label < static_cast<int>(divisions.size()); label++) {
    auto& side = positions[divisions[label]] < positions[hub] ? left : right;
    side.emplace_back(positions[divisions[label]], label);
  }
  std::sort(left.begin(), left.end());
  std::sort(right.rbegin(), right.rend());
  std::vector<int> rotation;
  for (const auto& [_, label] : left)
    rotation.push_back(label);
  for (const auto& [_, label] : right)
    rotation.push_back(label);
  return normalizedRotation(std::move(rotation));
}

/// Return the edge indices incident with one vertex.
std::vector<int> incidentEdges(const InputGraph& graph, const int vertex) {
  std::vector<int> result;
  for (int edge = 0; edge < static_cast<int>(graph.edges.size()); edge++)
    if (graph.edges[edge].first == vertex || graph.edges[edge].second == vertex)
      result.push_back(edge);
  return result;
}

void verifySixHub(
    const std::vector<EdgeTy>& edges, const AdjListTy& adjList,
    Params& params) {
  const int hub = static_cast<int>(
      std::find_if(adjList.begin(), adjList.end(),
                   [](const auto& neighbors) { return neighbors.size() == 6; }) -
      adjList.begin());
  CHECK(hub < static_cast<int>(adjList.size()), "degree-6 quotient has no hub");
  InputGraph graph(static_cast<int>(adjList.size()), edges, {});
  const std::vector<int> attachments = incidentEdges(graph, hub);
  CHECK(attachments.size() == 6, "degree-6 quotient has invalid attachments");
  verifyDrawingWithFixedEdges(params, graph, attachments);
}

std::string findSevenHubDrawing(
    const std::vector<EdgeTy>& edges, const AdjListTy& adjList,
    const Params& params, const std::set<std::string>& forbidden) {
  const int hub = static_cast<int>(
      std::find_if(adjList.begin(), adjList.end(),
                   [](const auto& neighbors) { return neighbors.size() == 7; }) -
      adjList.begin());
  CHECK(hub < static_cast<int>(adjList.size()), "degree-7 quotient has no hub");
  InputGraph graph(static_cast<int>(adjList.size()), edges, {});
  const std::vector<int> attachments = incidentEdges(graph, hub);
  CHECK(attachments.size() == 7, "degree-7 quotient has invalid attachments");
  std::vector<int> divisions;
  for (int edge : attachments)
    divisions.push_back(graph.n + edge);

  initCrossablePairs(params, graph);
  SATModel model;
  encodeStackPlanar(model, graph, params);
  std::vector<int> order = {hub};
  order.insert(order.end(), divisions.begin(), divisions.end());
  std::sort(order.begin(), order.end());
  int forbiddenLinearOrders = 0;
  do {
    std::vector<int> positions(graph.n + static_cast<int>(graph.edges.size()));
    for (int position = 0; position < static_cast<int>(order.size()); position++)
      positions[order[position]] = position;
    if (!forbidden.count(hubRotation(positions, hub, divisions)))
      continue;
    forbiddenLinearOrders++;
    MClause clause;
    for (int index = 0; index + 1 < static_cast<int>(order.size()); index++)
      clause.vars.push_back(
          model.getRelVar(order[index], order[index + 1], false));
    model.addClause(std::move(clause));
  } while (std::next_permutation(order.begin(), order.end()));
  CHECK(
      forbiddenLinearOrders == 112 * static_cast<int>(forbidden.size()),
      "degree-7 rotation encoding has an incomplete forbidden-order orbit");

  Solver solver;
  initSATSolver(params, graph, model, solver);
  vec<Lit> assumptions;
  for (int edge : attachments)
    assumptions.push(
        model.getSolverLit(model.getCross1Var(graph.n + edge, false)));
  requireSatisfiable(
      solver,
      solver.solveLimited(assumptions),
      "failed degree-7 quotient requirement");
  const Result drawing = verifiedDrawing(params, graph, model, solver);
  std::vector<int> positions(graph.n + static_cast<int>(graph.edges.size()), -1);
  for (int position = 0; position < static_cast<int>(drawing.order.size()); position++)
    for (int vertex : drawing.order[position])
      positions[vertex] = position;
  const std::string rotation = hubRotation(positions, hub, divisions);
  CHECK(!forbidden.count(rotation),
        "degree-7 quotient has forbidden rotation %s", rotation.c_str());
  return rotation;
}

void verifySevenHub(
    const std::vector<EdgeTy>& edges, const AdjListTy& adjList,
    const Params& params) {
  const std::set<std::string> patchForbidden = {
      "0245136", "0246135", "0264153", "0315426", "0316425"};
  std::vector<std::vector<int>> assignments;
  std::vector<int> assignment = {0, 1, 2, 3, 4, 5, 6};
  do {
    assignments.push_back(assignment);
  } while (std::next_permutation(assignment.begin(), assignment.end()));
  std::vector<bool> covered(assignments.size(), false);

  auto patchRotation = [](const std::string& quotientRotation,
                          const std::vector<int>& assignment) {
    std::vector<int> inverse(7);
    for (int patch = 0; patch < 7; patch++)
      inverse[assignment[patch]] = patch;
    std::vector<int> rotation;
    for (char label : quotientRotation)
      rotation.push_back(inverse[label - '0']);
    return normalizedRotation(std::move(rotation));
  };

  std::set<std::string> forbidden;
  while (true) {
    const std::string rotation =
        findSevenHubDrawing(edges, adjList, params, forbidden);
    for (int index = 0; index < static_cast<int>(assignments.size()); index++)
      if (!patchForbidden.count(patchRotation(rotation, assignments[index])))
        covered[index] = true;
    const auto uncovered = std::find(covered.begin(), covered.end(), false);
    if (uncovered == covered.end())
      return;
    const std::vector<int>& selected =
        assignments[static_cast<int>(uncovered - covered.begin())];
    forbidden.clear();
    for (const std::string& patchOrder : patchForbidden) {
      std::vector<int> mapped;
      for (char label : patchOrder)
        mapped.push_back(selected[label - '0']);
      forbidden.insert(normalizedRotation(std::move(mapped)));
    }
  }
}

/// Verify the cyclic-order-independent conditions for expanding a contracted
/// 5-cycle while preserving any prescribed edge of the original graph.
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
      CoverageRequirement{attachments, {}, 0}};
  for (int edge = 0; edge < static_cast<int>(graph.edges.size()); edge++) {
    std::vector<int> candidates = attachments;
    if (contains(attachments, edge))
      candidates.erase(std::find(candidates.begin(), candidates.end(), edge));
    requirements.push_back(
        CoverageRequirement{{edge}, std::move(candidates), 3});
  }
  CHECK(requirements.size() == 35,
        "unexpected number of 5-cycle-core requirements");
  verifyFlexibility(params, graph, 0, {}, true, requirements);
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

/// Verify Claim 3 directly when the input has no 5-cycle.
void verifyCubic26(const int n, const std::vector<EdgeTy>& edges, const Params& params) {
  CHECK(n == 26);
  InputGraph graph(n, edges, {});
  verifyFlexibility(params, graph, 1);
}

/// Remove the pendant marker and verify the uncrossed three-edge star needed
/// by the order-28 5-cycle-to-path expansion.
void verifyFiveCyclePathStar(const AdjListTy& adjList, Params& params) {
  const int n = static_cast<int>(adjList.size());
  CHECK(n == 27, "marked path reduction has the wrong order");
  int marker = -1;
  for (int vertex = 0; vertex < n; vertex++)
    if (adjList[vertex].size() == 1)
      marker = vertex;
  CHECK(marker != -1, "marked path reduction has no marker");
  const int target = adjList[marker].front();

  std::vector<int> remap(n, -1);
  int baseN = 0;
  for (int vertex = 0; vertex < n; vertex++)
    if (vertex != marker)
      remap[vertex] = baseN++;
  std::vector<EdgeTy> baseEdges;
  for (int first = 0; first < n; first++)
    for (int second : adjList[first])
      if (first < second && first != marker && second != marker)
        baseEdges.emplace_back(
            std::min(remap[first], remap[second]),
            std::max(remap[first], remap[second]));
  InputGraph graph(baseN, baseEdges, {});
  CHECK(baseN == 26 && graph.edges.size() == 39,
        "marked path reduction has an invalid base");
  const std::vector<int> star = incidentEdges(graph, remap[target]);
  CHECK(star.size() == 3, "marked path reduction has an invalid star");
  verifyDrawingWithFixedEdges(params, graph, star);
}

/// Expand a degree-five Claim 3 core in every cyclic order, retaining the
/// biconnected, nonplanar, girth-at-least-five parents covered by Claim 3.
std::vector<AdjListTy> claim3Parents(
    const std::vector<EdgeTy>& coreEdges, const AdjListTy& coreAdj) {
  const int coreN = static_cast<int>(coreAdj.size());
  CHECK(coreN == 22, "Claim 3b expects order-22 cores");
  const int hub = static_cast<int>(
      std::find_if(coreAdj.begin(), coreAdj.end(),
                   [](const auto& neighbors) { return neighbors.size() == 5; }) -
      coreAdj.begin());
  CHECK(hub < coreN, "Claim 3b core has no degree-five hub");

  std::vector<int> remap(coreN, -1);
  int reducedN = 0;
  for (int vertex = 0; vertex < coreN; vertex++)
    if (vertex != hub)
      remap[vertex] = reducedN++;
  CHECK(reducedN == 21, "Claim 3b core remapping failed");

  std::vector<EdgeTy> baseEdges;
  for (const auto& [first, second] : coreEdges)
    if (first != hub && second != hub)
      baseEdges.emplace_back(remap[first], remap[second]);

  std::vector<int> attachments;
  for (int neighbor : coreAdj[hub])
    attachments.push_back(remap[neighbor]);
  std::sort(attachments.begin(), attachments.end());
  CHECK(attachments.size() == 5, "Claim 3b has invalid attachments");

  const int parentN = 26;
  const int cycleStart = 21;
  std::vector<AdjListTy> parents;
  std::vector<int> order = attachments;
  int numCyclicOrders = 0;
  do {
    // Fix rotation and reflection.  The remaining 4!/2 orders are the 12
    // cyclic orders of the five incident hub edges.
    if (order.front() != attachments.front() || order[1] > order.back())
      continue;
    numCyclicOrders++;
    std::vector<EdgeTy> edges = baseEdges;
    for (int position = 0; position < 5; position++) {
      const int cycleVertex = cycleStart + position;
      edges.emplace_back(
          std::min(order[position], cycleVertex),
          std::max(order[position], cycleVertex));
      const int nextCycleVertex = cycleStart + (position + 1) % 5;
      edges.emplace_back(
          std::min(cycleVertex, nextCycleVertex),
          std::max(cycleVertex, nextCycleVertex));
    }
    std::sort(edges.begin(), edges.end());
    const AdjListTy parent = edges_to_adj(parentN, edges);
    CHECK(
        std::all_of(
            parent.begin(), parent.end(),
            [](const auto& neighbors) { return neighbors.size() == 3; }),
        "Claim 3b expansion is not cubic");
    if (!is2Connected(parentN, edges) ||
        computeGirth(parentN, edges) < 5 ||
        isPlanar(parentN, edges, 0))
      continue;
    parents.push_back(parent);
  } while (std::next_permutation(order.begin(), order.end()));
  CHECK(numCyclicOrders == 12, "Claim 3b missed a cyclic order");
  CHECK(!parents.empty(), "Claim 3b core has no relevant parent");
  CHECK(parents.size() <= 12, "Claim 3b generated too many parents");
  return parents;
}

void verifyClaim3b(
    GraphList& graphs, const Params& params, std::ofstream& residue,
    const std::string& residueFilename) {
  const int numRecords = static_cast<int>(graphs.size());
  const auto start = std::chrono::steady_clock::now();
  int numCompletedRecords = 0;
  int numParents = 0;
  int num1Flexible = 0;
  int numUnknown = 0;

  for (int recordIndex = 0; recordIndex < numRecords; recordIndex++) {
    const auto [graphName, recordAdj] = graphs.next();
    const std::vector<EdgeTy> recordEdges = adj_to_edges(recordAdj);
    const VerificationInput input =
        validateInput(recordAdj, recordEdges);
    std::vector<AdjListTy> parents;
    if (input == VerificationInput::FIVE_CYCLE_CORE) {
      parents = claim3Parents(recordEdges, recordAdj);
    } else {
      CHECK(
          input == VerificationInput::CUBIC &&
              recordAdj.size() == 26,
          "Claim 3b accepts only Claim 3 cores or order-26 graphs");
      parents.push_back(recordAdj);
    }
    bool completed = true;
    for (int parentIndex = 0;
         parentIndex < static_cast<int>(parents.size()); parentIndex++) {
      const AdjListTy& parentAdj = parents[parentIndex];
      const std::vector<EdgeTy> parentEdges = adj_to_edges(parentAdj);
      numParents++;
      try {
        verifyCubic26(26, parentEdges, params);
        num1Flexible++;
      } catch (const VerificationTimedOut& timeout) {
        CHECK(residue.good(), "Claim 3b timeout has no residue file");
        residue << graph6(parentAdj) << '\n';
        residue.flush();
        CHECK(residue.good(), "failed to write Claim 3b residue");
        numUnknown++;
        completed = false;
        LOG(
            "Claim 3b timed out for %s parent %d/%d after %d seconds "
            "during %s (%'llu conflicts, %'llu decisions); wrote it to %s",
            graphName.c_str(), parentIndex + 1,
            static_cast<int>(parents.size()), params.timeout,
            timeout.stage.c_str(), timeout.conflicts, timeout.decisions,
            residueFilename.c_str());
      }
    }
    if (completed)
      numCompletedRecords++;
    LOG_EVERY_MS(
        30000, "Claim 3b processed %'d of %'d records",
        recordIndex + 1, numRecords);
  }

  CHECK(
      num1Flexible + numUnknown == numParents,
      "Claim 3b counters do not account for every parent");
  LOG(
      "Claim 3b processed %'d records and %'d parents in %s",
      numRecords, numParents,
      ms_to_str(start, std::chrono::steady_clock::now()).c_str());
  LOG(
      "#claim3b-records = %'d; #claim3b-completed = %'d; "
      "#claim3b-parents = %'d; #1-flexible = %'d; #unknown = %'d",
      numRecords, numCompletedRecords, numParents,
      num1Flexible, numUnknown);
}

}  // namespace

/// Verify the claim or claims prescribed by the order of each input graph.
void verifyCubic(CMDOptions& options) {
  Params params;
  params.cubicVerification = true;
  params.crossPriority = true;
  params.timeout = options.getInt("-timeout");
  CHECK(params.timeout >= 0, "verification timeout must be nonnegative");

  const std::string residueFilename =
      options.getCustomValue("verify-cubic-residue");
  CHECK(
      params.timeout == 0 || !residueFilename.empty(),
      "-verify-cubic with -timeout requires "
      "-Cverify-cubic-residue=FILE");
  std::ofstream residue;
  if (!residueFilename.empty()) {
    residue.open(residueFilename, std::ios::out | std::ios::trunc);
    CHECK(residue.good(), "cannot open cubic-verification residue file");
  }

  auto graphs = genGraphs(options);
  if (options.hasCustomValue("claim3b")) {
    verifyClaim3b(
        *graphs, params, residue, residueFilename);
    return;
  }
  const int numGraphs = graphs->size();
  const auto start = std::chrono::steady_clock::now();

  int numPlanar = 0;
  int num1Planar = 0;
  int num3Flexible = 0;
  int num2Flexible = 0;
  int num1Flexible = 0;
  int numFiveCycleCores = 0;
  int numFiveCyclePathStars = 0;
  int numSixHubs = 0;
  int numSevenHubs = 0;
  int numUnknown = 0;

  for (int i = 0; i < numGraphs; i++) {
    const auto [graphName, adj] = graphs->next();
    const int n = static_cast<int>(adj.size());
    const std::vector<EdgeTy> edges = adj_to_edges(adj);
    const VerificationInput input = validateInput(adj, edges);
    if (input == VerificationInput::FIVE_CYCLE_CORE)
      numFiveCycleCores++;
    else if (input == VerificationInput::FIVE_CYCLE_PATH_STAR)
      numFiveCyclePathStars++;
    else if (input == VerificationInput::SIX_HUB)
      numSixHubs++;
    else if (input == VerificationInput::SEVEN_HUB)
      numSevenHubs++;

    try {
      if (input == VerificationInput::SEVEN_HUB) {
        verifySevenHub(edges, adj, params);
        num1Planar++;
      } else if (isPlanar(n, edges, 0)) {
        numPlanar++;
      } else if (input == VerificationInput::FIVE_CYCLE_CORE) {
        verifyFiveCycleCore(edges, adj, params);
        num1Planar++;
      } else if (input == VerificationInput::FIVE_CYCLE_PATH_STAR) {
        verifyFiveCyclePathStar(adj, params);
        num1Planar++;
      } else if (input == VerificationInput::SIX_HUB) {
        verifySixHub(edges, adj, params);
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
        num1Flexible++;
        num1Planar++;
      } else {
        CHECK(false, "unsupported cubic verification input");
      }
    } catch (const VerificationTimedOut& timeout) {
      CHECK(residue.good(), "cubic-verification timeout has no residue file");
      residue << graph6(adj) << '\n';
      residue.flush();
      CHECK(residue.good(), "failed to write cubic-verification residue");
      numUnknown++;
      LOG(
          "verification timed out for %s after %d seconds during %s "
          "(%'llu conflicts, %'llu decisions); wrote it to %s",
          graphName.c_str(), params.timeout, timeout.stage.c_str(),
          timeout.conflicts, timeout.decisions, residueFilename.c_str());
    }

    LOG_EVERY_MS(30000, "verified %'d of %'d graphs", i + 1, numGraphs);
  }

  CHECK(
      numPlanar + num1Planar + numUnknown == numGraphs,
      "cubic-verification counters do not account for every input");
  LOG("processed %'d graphs in %s", numGraphs, ms_to_str(start, std::chrono::steady_clock::now()).c_str());
  LOG("#planar = %'d; #1-planar = %'d", numPlanar, num1Planar);
  LOG("#unknown = %'d", numUnknown);
  LOG("#3-flexible = %'d; #2-flexible = %'d", num3Flexible, num2Flexible);
  LOG("#1-flexible = %'d; #5-cycle-cores = %'d; "
      "#5-cycle-path-stars = %'d; #degree-6-hubs = %'d; "
      "#degree-7-hubs = %'d",
      num1Flexible, numFiveCycleCores, numFiveCyclePathStars,
      numSixHubs, numSevenHubs);
}
