#include "cmd_options.h"
#include "graph_algorithms.h"
#include "io.h"
#include "logging.h"
#include "one_planar.h"

#include <algorithm>
#include <chrono>
#include <climits>
#include <functional>
#include <memory>
#include <vector>

std::unique_ptr<GraphList> genGraphs(CMDOptions& options);

namespace {

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
  FIVE_CYCLE_STAR_CORE
};

/// Check that an input record is cubic or is a marked reduced graph.
VerificationInput validateInput(const AdjListTy& adjList, const std::vector<EdgeTy>& edges) {
  const int n = static_cast<int>(adjList.size());
  CHECK(isConnected(n, edges), "verification inputs contain connected graphs");

  if (n == 23) {
    int leaf = -1;
    int hub = -1;
    bool valid = true;
    for (int vertex = 0; vertex < n; vertex++) {
      const int degree = static_cast<int>(adjList[vertex].size());
      int* location = nullptr;
      if (degree == 1)
        location = &leaf;
      else if (degree == 6)
        location = &hub;
      else if (degree != 3)
        valid = false;
      if (location != nullptr) {
        if (*location != -1)
          valid = false;
        *location = vertex;
      }
    }
    if (valid && leaf != -1 && hub != -1 &&
        contains(adjList[leaf], hub))
      return VerificationInput::FIVE_CYCLE_STAR_CORE;
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
    CHECK(n == 22 && exceptionalDegree == 5,
          "invalid contracted graph for Claim 3");
    return VerificationInput::FIVE_CYCLE_CORE;
  }

  CHECK(4 <= n && n <= 28 && n % 2 == 0,
        "verification expects an even order from 4 through 28");
  if (n >= 24)
    CHECK(computeGirth(n, edges) >= 5,
          "orders 24 through 28 require girth at least five");
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

/// Remove the pendant marker and verify the all-attachments requirement used
/// jointly by the order-26 contraction and the order-28 path reduction.
void verifyFiveCycleStarCore(
    const AdjListTy& adjList, const Params& params) {
  const int n = static_cast<int>(adjList.size());
  int markedHub = -1;
  std::vector<int> remap(n, -1);
  int baseN = 0;
  for (int vertex = 0; vertex < n; vertex++) {
    const int degree = static_cast<int>(adjList[vertex].size());
    if (degree == 6)
      markedHub = vertex;
    if (degree >= 3)
      remap[vertex] = baseN++;
  }
  CHECK(baseN == 22 && markedHub != -1,
        "invalid marked 5-cycle core");

  std::vector<EdgeTy> baseEdges;
  for (int first = 0; first < n; first++)
    for (int second : adjList[first])
      if (first < second && remap[first] != -1 && remap[second] != -1)
        baseEdges.push_back(make_edge(remap[first], remap[second]));
  InputGraph graph(baseN, baseEdges, {});
  CHECK(graph.edges.size() == 34,
        "marked 5-cycle core has an invalid base");
  const std::vector<int> attachments =
      incidentEdges(graph, remap[markedHub]);
  CHECK(attachments.size() == 5,
        "marked 5-cycle core has an invalid attachment set");
  verifyDrawing(params, graph, attachments);
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

/// Directly verify a residual order-28 graph from the joint reduction.
void verifyCubic28(
    const int n, const std::vector<EdgeTy>& edges, const Params& params) {
  CHECK(n == 28);
  InputGraph graph(n, edges, {});
  verifyDrawing(params, graph, {});
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
  int num1Flexible = 0;
  int numFiveCycleCores = 0;
  int numFiveCycleStarCores = 0;

  for (int i = 0; i < numGraphs; i++) {
    const auto& adj = graphs->next().second;
    const int n = static_cast<int>(adj.size());
    const std::vector<EdgeTy> edges = adj_to_edges(adj);
    const VerificationInput input = validateInput(adj, edges);
    if (input == VerificationInput::FIVE_CYCLE_CORE)
      numFiveCycleCores++;
    else if (input == VerificationInput::FIVE_CYCLE_STAR_CORE)
      numFiveCycleStarCores++;

    if (isPlanar(n, edges, 0)) {
      numPlanar++;
    } else if (input == VerificationInput::FIVE_CYCLE_CORE) {
      verifyFiveCycleCore(edges, adj, params);
      num1Planar++;
    } else if (input == VerificationInput::FIVE_CYCLE_STAR_CORE) {
      verifyFiveCycleStarCore(adj, params);
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
    } else if (n == 28) {
      verifyCubic28(n, edges, params);
      num1Planar++;
    } else {
      CHECK(false, "unsupported cubic verification input");
    }

    LOG_EVERY_MS(30000, "verified %'d of %'d graphs", i + 1, numGraphs);
  }

  LOG("processed %'d graphs in %s", numGraphs, ms_to_str(start, std::chrono::steady_clock::now()).c_str());
  LOG("#planar = %'d; #1-planar = %'d", numPlanar, num1Planar);
  LOG("#3-flexible = %'d; #2-flexible = %'d", num3Flexible, num2Flexible);
  LOG("#1-flexible = %'d; #5-cycle-cores = %'d; "
      "#5-cycle-star-cores = %'d",
      num1Flexible, numFiveCycleCores, numFiveCycleStarCores);
}
