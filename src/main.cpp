#include "logging.h"
#include "cmd_options.h"
#include "io.h"
#include "graph_algorithms.h"
#include "graph_generator.h"
#include "one_planar.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <queue>
#include <chrono>
#include <ctime>
#include <optional>
#include <regex>

using namespace std;

/// Prepare command-line arguments
void prepareOptions(CMDOptions& options) {
  options.setUsageMessage("Usage: oops -i=[input_graph] [options]");

  // IO
  options.addAllowedOption("-i", "File name with input graph(s); supported formats are dot/gml/graphml/s6/g6");
  options.addAllowedOption("-o", "", "File name to output solution; supported formats are txt/dot/gml/svg");
  options.addAllowedOption("-part", "", "The part of the input to process in the form of part_idx/num_parts");
  options.addAllowedOption("-max-n", "-1", "The maximum number of vertices in the processed graph");
  options.addAllowedOption("-max-degree", "-1", "Max vertex degree in the processed graph");
  options.addAllowedOption("-directed", "false", "Whether the input graph is directed");

  // Debug
  options.addAllowedOption("-verbose", "1", "Verbosity level");
  options.addAllowedOption("-seed", "0", "Random seed");

  // Solver options
  options.addAllowedOption("-solver", "stack", "Type of the solver to use");
  options.addAllowedValue("-solver", "stack");
  options.addAllowedValue("-solver", "move");
  options.addAllowedValue("-solver", "brute-force");
  options.addAllowedOption("-skewness", "1", "Maximum value of skewness");
  options.addAllowedOption("-satsuma", "false", "Whether to apply Satsuma symmetry detection");
  options.addAllowedOption("-breakID", "false", "Whether to apply BreakID symmetry detection");
  options.addAllowedOption("-timeout", "0", "Maximum time (in seconds) to solve SAT");
  options.addAllowedOption("-sparsify", "false", "Try to eliminate crossings in the resulting solution");
  options.addAllowedOption("-sat", "true", "Add auxiliary constraints to speedup finding 1-planar instances");
  options.addAllowedOption("-unsat", "false", "Add auxiliary constraints to speedup finding non-1-planar instances");
  options.addAllowedOption("-ic", "false", "Enforce IC constraints");
  options.addAllowedOption("-nic", "false", "Enforce NIC constraints");

  // External SAT solver
  options.addAllowedOption("-dimacs", "", "Output dimacs file");
  options.addAllowedOption("-dimacs-result", "", "Input dimacs file with a solution");

  // Graph generation
  options.addAllowedOption("-graphs", "1", "How many graphs to generate");
  options.addAllowedOption("-n", "20", "The number of vertices in the generated graph");

  // Experimental
  options.addAllowedOption("-forbid-crossings", "false", "[Experimental] Forbid all crossings");
  options.addAllowedOption("-skip-reducible-triangles", "false", "[Experimental] Skip reducible triangles");
  options.addAllowedOption("-swap-constraints", "", "[Experimental] Add swap constraints: num_pairs/num_reorder");
}

/// Check if 1-planarity test can be skipped because of graphs size or density
bool skipTestingBySize(const InputGraph& graph, const int verbose, ResultCodeTy& resultCode) {
  const int n = graph.n;
  const int m = (int)graph.edges.size();
  // check planarity
  if (n < 5 || m < 9) {
    LOG_IF(verbose, "the graph is planar due to size");
    resultCode = ResultCodeTy::SAT;
    return true;
  }
  // check 1-planarity
  if (n < 7 || m < 18) {
    LOG_IF(verbose, "the graph is 1-planar due to size");
    resultCode = ResultCodeTy::SAT;
    return true;
  }
  CHECK(n >= 7);

  // density of general
  if (m > 4 * n - 8) {
    LOG_IF(verbose, "the graph is not 1-planar due to edge density (%d > %d)", m, 4 * n - 8);
    resultCode = ResultCodeTy::UNSAT;
    return true;
  }
  // density of bipartite
  if (m > 3 * n - 8 - n % 2 && isBipartite(n, graph.edges)) {
    LOG_IF(verbose, "the graph is not 1-planar due to edge density (%d > %d)", m, 3 * n - 8 - n % 2);
    resultCode = ResultCodeTy::UNSAT;
    return true;
  }
  return false;
}

ResultCodeTy isOnePlanar(
    CMDOptions& options, 
    const Params& params,
    const InputGraph& graph,
    const bool splitIntoBiconnected, 
    const std::string& graphName) {
  const int verbose = options.getInt("-verbose");
  const bool EnableEarlyExit = !graph.isDirected() && !params.useIC && !params.alwaysCreateSolution;

  const int n = graph.n;
  const auto& edges = graph.edges;

  LOG_IF(verbose, "testing 1-planarity of a graph with |V| = %d and |E| = %d", n, edges.size());
  if (verbose) {
    printGraphStats(n, edges, std::vector<bool>());
  }

  // print input
  printInput(options.getStr("-o"), graphName, graph, verbose);

  // split into biconnected components
  if (EnableEarlyExit && splitIntoBiconnected) {
    auto biComponents = biconnectedComponents(n, edges);
    LOG_IF(verbose, "the number of bi-componets: %d", biComponents.size());
    if (biComponents.size() > 1) {
      std::sort(
          biComponents.begin(), 
          biComponents.end(), 
          [](const std::vector<EdgeTy>& l, const std::vector<EdgeTy>& r) { 
              return l.size() < r.size(); 
      });
      for (auto& compEdges : biComponents) {
        // create a graph for the component
        InputGraph compGraph(compEdges);
        CHECK(compGraph.n == 2 || compGraph.minDegree() >= 2, 
          "the bi-connected component contains low-degree vertices");
        // start test for the component
        ResultCodeTy res = isOnePlanar(options, params, compGraph, false, graphName);
        if (res != ResultCodeTy::SAT) {
          return res;
        }
      }
      return ResultCodeTy::SAT;
    }
  }

  // check if planarity test can be skipped
  ResultCodeTy resultCode = ResultCodeTy::ERROR;
  if (EnableEarlyExit && skipTestingBySize(graph, verbose, resultCode)) {
    CHECK(resultCode == ResultCodeTy::SAT || resultCode == ResultCodeTy::UNSAT);
    return resultCode;
  }

  // init pairs for edge crossings
  initCrossablePairs(params, graph);

  if (EnableEarlyExit && params.solverType == SolverType::BRUTE_FORCE) {
    CHECK(EnableEarlyExit, "cannot combine specified parameters with brute-force");
    Result result = bruteForce(params, graph);
    return result.code;
  }

  // check skewness (max number of crossings)
  if (EnableEarlyExit) {
    const int maxSkewness = options.getInt("-skewness");
    const int skewness = computeSkewness(graph, verbose, maxSkewness);
    if (skewness <= maxSkewness) {
      LOG_IF(verbose, "the graph has skewness %d (<= %d)", skewness, maxSkewness);
      return ResultCodeTy::SAT;
    }
    LOG_IF(verbose, "the graph has skewness >= %d", skewness);
  }

  auto startTimeSAT = chrono::steady_clock::now();
  Result result = runSolver(params, graph);
  auto endTimeSAT = chrono::steady_clock::now();
  LOG_IF(verbose >= 1, "SAT solving took %s for %s [timeout = %d sec]", 
      ms_to_str(startTimeSAT, endTimeSAT).c_str(), 
      graphName.c_str(), 
      params.timeout
  );

  if (result.code != ResultCodeTy::SAT) {
    return result.code;
  }

  // Verify that the result is planar (after subdivision of crossings)
  if (!isPlanarWithCrossings(graph, result.crossings)) {
    LOG(TextColor::red, "the constructed embedding is not 1-planar");
    ERROR(false);
  }

  // Try to remove crossings, if possible
  const bool sparsify = options.getBool("-sparsify");
  if (sparsify) {
    minimizeCrossings(graph, result, verbose);
  }

  // print output
  printOutput(options.getStr("-o"), graphName, graph, result, verbose);

  return result.code;
}

std::unique_ptr<GraphList> genGraphs(CMDOptions& options) {  
  const int verbose = options.getInt("-verbose");
  const std::string input = options.getStr("-i");

  // gen random
  if (input.substr(0, 4) == "gen-") {
    const string graphClass = input.substr(4);
    const int numGraphs = options.getInt("-graphs");
    const int startSeed = options.getInt("-seed");
    const int n0 = options.getInt("-n");

    std::vector<std::pair<std::string, AdjListTy>> graphs;
    for (int t = 0; t < numGraphs; t++) {
      const int seed = t + startSeed;
      Rand::setSeed(seed);

      int n = n0;
      std::vector<EdgeTy> edges;
      genByClass(graphClass, n, edges);
      LOG_IF(verbose, "generated %s graph with |V| = %d  |E| = %d", graphClass.c_str(), n, edges.size());

      graphs.push_back({
        "seed_" + to_string(t), 
        edges_to_adj(n, edges)
      });
    }
    return std::make_unique<GraphListRaw>(graphs);
  }

  const std::string filename = options.getStr("-i");
  const std::string extension = filename.substr(filename.find_last_of(".") + 1);

  const int maxN = options.getInt("-max-n");
  const int maxD = options.getInt("-max-degree");
  auto graphFilter = [maxN,maxD](const int n, const std::vector<EdgeTy>& edges) -> bool {
    if (maxN != -1 && n > maxN)
      return false;
    if (maxD != -1 && maxDegree(n, edges) > maxD)
      return false;
    return true;
  };

  // read a list from a file
  if (extension == "g6") {
    const std::string part = options.getStr("-part");
    return std::make_unique<GraphListG6>(filename, part, graphFilter);
  }
  if (extension == "s6") {
    // TODO: replace with lazy list
    const std::string part = options.getStr("-part");
    return std::make_unique<GraphListRaw>(readS6(filename, part, graphFilter));
  }
  if (extension == "cfg") {
    const std::string part = options.getStr("-part");
    return std::make_unique<GraphListRaw>(readCfg(filename, part, graphFilter));
  }

  // read a single instance
  IOGraph ioGraph;
  if (extension == "dot") {
    LOG("parsing dot graph from '%s'", filename.c_str());
    GraphParser parser;
    CHECK(parser.readDotGraph(filename, ioGraph), "dot file parsing failed");
  } else if (extension == "gml") {
    LOG("parsing gml graph from '%s'", filename.c_str());
    GraphParser parser;
    CHECK(parser.readGmlGraph(filename, ioGraph), "gml file parsing failed");
  } else if (extension == "graphml") {
    LOG("parsing graphml graph from '%s'", filename.c_str());
    GraphParser parser;
    CHECK(parser.readGraphmlGraph(filename, ioGraph), "graphml file parsing failed");
  } else {
    ERROR("unknown file extension: " + extension);
  }

  int n;
  std::vector<EdgeTy> edges;
  ioGraph.extractEdges(n, edges);
  const int numGraphs = options.getInt("-graphs");
  std::vector<std::pair<std::string, AdjListTy>> graphs;
  for (int i = 0; i < numGraphs; i++) {
    std::string graphName = (numGraphs == 1 ? filename : filename + "_" + to_string(i));
    graphs.push_back({
      filename,
      edges_to_adj(n, edges)
    });
    // TODO: all but the first graph have random directions
  }

  return std::make_unique<GraphListRaw>(graphs);
}

// TODO: merge with the above method
void genDirections(CMDOptions& options, const int n, std::vector<EdgeTy>& edges, std::vector<bool>& directions) {
  // read from file
  // const std::string filename = options.getStr("-i");
  // const std::string extension = filename.substr(filename.find_last_of(".") + 1);
  // CHECK(extension == "gml");

  // // read a single instance
  // PlanarGraph graph;
  // LOG("parsing gml graph from '%s'", filename.c_str());
  // GraphParser parser;
  // CHECK(parser.readGmlGraph(filename, graph.ioGraph), "gml file parsing failed");

  // CHECK(n == (int)graph.ioGraph.nodes.size());
  // CHECK(edges.size() == graph.ioGraph.edges.size());

  // edges.clear();
  // directions.clear();
  // std::vector<EdgeTy> dirEdges;
  // for (auto& e : graph.ioGraph.edges) {
  //   const auto u = graph.ioGraph.getNode(e.source);
  //   const auto v = graph.ioGraph.getNode(e.target);
  //   CHECK(u->index != v->index);
  //   // LOG("parsed edge (%d -> %d)", u->index, v->index);

  //   dirEdges.push_back(EdgeTy(u->index, v->index));
  //   if (u->index < v->index) {
  //     edges.push_back(EdgeTy(u->index, v->index));
  //     directions.push_back(true);
  //   } else {
  //     edges.push_back(EdgeTy(v->index, u->index));
  //     directions.push_back(false);
  //   }
  // }
  // CHECK(isAcyclic(n, dirEdges, 2), "the GML graph is not acyclic");

  // gen random
  directions.reserve(edges.size());
  auto perm = Rand::permutation(n);
  // LOG("vertex perm: %s", to_string(perm).c_str());
  for (const auto& [u, v] : edges) {
    CHECK(u != v);
    if (perm[u] < perm[v])
      directions.push_back(true);
    else
      directions.push_back(false);
  }
}

/// Create parameters for SAT solving
void initSATParams(CMDOptions& options, Params& params) {
  params.verbose = options.getInt("-verbose");;
  params.timeout = options.getInt("-timeout");;
  //params.directed = options.getBool("-directed");
  params.applyBreakID = options.getBool("-breakID");
  params.applySatsuma = options.getBool("-satsuma");
  params.modelFile = options.getStr("-dimacs");
  params.resultFile = options.getStr("-dimacs-result");

  params.useSATConstraints = options.getBool("-sat");
  if (options.getBool("-unsat")) {
    params.useUNSATConstraints = true;
    params.useSATConstraints = true;
    params.applySatsuma = true;
  }
  if (options.getBool("-ic")) {
    params.useIC = true;
    params.useSATConstraints = true;
  }
  if (options.getBool("-nic")) {
    params.useNIC = true;
    params.useSATConstraints = true;
  }

  if (options.getStr("-solver") == "move") {
    params.solverType = SolverType::MOVE;
    params.useSATConstraints = true;
  } else if (options.getStr("-solver") == "brute-force") {
    params.solverType = SolverType::BRUTE_FORCE;
  } else {
    CHECK(options.getStr("-solver") == "stack");
    params.solverType = SolverType::STACK;
  }

  params.forbidCrossings = options.getBool("-forbid-crossings");
  params.swapConstraints = options.getStr("-swap-constraints");

  const std::string outFile = options.getStr("-o");
  if (outFile != "") {
    const std::string extension = outFile.substr(outFile.find_last_of(".") + 1);
    CHECK(extension == "txt" || extension == "gml" || extension == "dot" || extension == "svg", 
          "unsupported output format: '%s'", extension.c_str());
    params.alwaysCreateSolution = true;
  }

  //CHECK(!params.directed || params.useMovePlanarity, "directed edges should be used with move-planarity");
  CHECK(!params.useIC || !params.useNIC, "cannot simultanosly use IC and NIC modes");
  CHECK(!params.useUNSATConstraints || params.useSATConstraints, "`-unsat` constraints should be used with `-sat`");
  CHECK(!params.useIC || params.useSATConstraints, "`-useIC` constraints should be used with `-sat`");
  CHECK(!params.useNIC || params.useSATConstraints, "`-useNIC` constraints should be used with `-sat`");
}

/// Gen a random 1-planar graph and verify the 1-planar SAT model
void testOnePlanar(CMDOptions& options) {
  const int startSeed = options.getInt("-seed");
  const int verbose = options.getInt("-verbose");
  const bool directed = options.getBool("-directed");

  Params params;
  initSATParams(options, params);

  auto graphs = genGraphs(options);
  const int numGraphs = graphs->size();
  // start time for every graph
  std::vector<chrono::time_point<chrono::steady_clock>> times;
  times.reserve(numGraphs + 1);
  times.push_back(chrono::steady_clock::now());
  // processing times of each instance
  std::vector<uint64_t> processingTimes;
  processingTimes.reserve(numGraphs);

  int numPlanar = 0;
  int num1Planar = 0;
  int numNon1Planar = 0;
  int numUnknown = 0;
  int numSkipped = 0;
  for (int t = 0; t < numGraphs; t++) {
    const int seed = t + startSeed;
    Rand::setSeed(seed);
    const auto& [graphName, graphAdj] = graphs->next();
    LOG_IF(verbose, "processing graph %d (%s)", t, graphName.c_str());

    const int n = (int)graphAdj.size();
    std::vector<EdgeTy> edges = adj_to_edges(graphAdj);
    std::vector<bool> directions;
    if (directed) {
      genDirections(options, n, edges, directions);
    }

    // only for cubic
    if (options.getBool("-skip-reducible-triangles") && hasReducibleTriangle(graphAdj)) {
      if (verbose)
        LOG(TextColor::green, "the graph contains reducible triangle");
      numSkipped++;
    } else {
      // test planarity
      if (directions.empty() && isPlanar(n, edges, 0)) {
        if (verbose)
          LOG(TextColor::green, "the graph is planar");
        numPlanar++;
      } else {
        // test 1-planarity
        CHECK(n >= 5, "the graph is too small");
        InputGraph graph(n, edges, directions);
        ResultCodeTy res = isOnePlanar(options, params, graph, true, graphName);
        if (res == ResultCodeTy::SAT) {
          if (verbose)
            LOG(TextColor::green, "graph '%s' (index %d) with |V| = %d and |E| = %d is 1-planar", graphName.c_str(), t, n, edges.size());
          num1Planar++;
        } else if (res == ResultCodeTy::UNSAT) {
          if (verbose)
            LOG(TextColor::red, "graph '%s' (index %d) with |V| = %d and |E| = %d is not 1-planar", graphName.c_str(), t, n, edges.size());
          numNon1Planar++;
        } else if (res == ResultCodeTy::TIMEOUT) {
          LOG(TextColor::red, "graph '%s' (index %d) with |V| = %d and |E| = %d timed out", graphName.c_str(), t, n, edges.size());
          numUnknown++;
        } else {
          ERROR("unreachable");
        }
      }
    }

    times.push_back(chrono::steady_clock::now());
    processingTimes.push_back(chrono::duration_cast<chrono::milliseconds>(times[t + 1] - times[t]).count());
    LOG_IF(verbose, "processed graph %d (%s) in %s\n", t, graphName.c_str(), ms_to_str(processingTimes.back()).c_str());
    if (t + 1 >= numGraphs) continue;

    const int prevIndex = t >= 1000 ? t - 1000 : 0;
    const auto durationMs = chrono::duration_cast<chrono::milliseconds>(times[t + 1] - times[prevIndex]).count();
    const double msPerGraph = double(durationMs) / (t + 1 - prevIndex);
    const uint64_t remainingMs = uint64_t(double(numGraphs - t - 1) * msPerGraph);
    LOG_EVERY_MS(
        30000, 
        "processed %'3d (%5.2lf%%) graphs with n = %d; expected completion in %s", 
        t + 1, 100.0 * (t + 1) / numGraphs, n, ms_to_str(remainingMs).c_str()
    );
    LOG_EVERY_MS(
        180000, 
        "#planar = %'d; #1-planar = %'d; #non-1-planar = %'d; #unknown = %'d; #skipped = %'d", 
        numPlanar, num1Planar, numNon1Planar, numUnknown, numSkipped
    );
  }
  CHECK((int)times.size() == numGraphs + 1);

  LOG("processed all %'d graphs; total runtime is %s; mean processing time is %zu ± %zu ms", 
      numGraphs, ms_to_str(times.front(), times.back()).c_str(),
      uint64_t(average(processingTimes)), uint64_t(confidence_interval(processingTimes)));
  LOG("#planar = %'d; #1-planar = %'d; #non-1-planar = %'d; #unknown = %'d; #skipped = %'d", 
      numPlanar, num1Planar, numNon1Planar, numUnknown, numSkipped);
}

int main(int argc, char* argv[]) {
  auto options = CMDOptions::create();

  try {
    prepareOptions(*options);
    options->parse(argc, argv);
    testOnePlanar(*options);
  } catch (int code) {
    return code;
  }

  return 0;
}
