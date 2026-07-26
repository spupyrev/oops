#include "graph_algorithms.h"

#include <nauty.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kMaximumOrder = 32;

struct Graph {
  int n = 0;
  std::vector<uint64_t> adj;

  explicit Graph(const int order = 0) : n(order), adj(order) {}

  bool hasEdge(const int first, const int second) const {
    return (adj[first] & (uint64_t(1) << second)) != 0;
  }

  void addEdge(const int first, const int second) {
    if (first == second || hasEdge(first, second))
      throw std::runtime_error("attempted to add a loop or parallel edge");
    adj[first] |= uint64_t(1) << second;
    adj[second] |= uint64_t(1) << first;
  }

  int degree(const int vertex) const {
    return __builtin_popcountll(adj[vertex]);
  }

  std::vector<int> neighbors(const int vertex) const {
    std::vector<int> result;
    for (int neighbor = 0; neighbor < n; neighbor++)
      if (hasEdge(vertex, neighbor))
        result.push_back(neighbor);
    return result;
  }

  std::vector<EdgeTy> edges() const {
    std::vector<EdgeTy> result;
    for (int first = 0; first < n; first++)
      for (int second = first + 1; second < n; second++)
        if (hasEdge(first, second))
          result.emplace_back(first, second);
    return result;
  }
};

Graph parseGraph6(const std::string& line) {
  if (line.empty() || line[0] == '>')
    throw std::runtime_error("expected one small graph6 record per line");
  const int n = static_cast<unsigned char>(line[0]) - 63;
  if (n < 1 || n > kMaximumOrder)
    throw std::runtime_error("only graph6 records of order at most 32 are supported");
  Graph graph(n);
  int byte = 1;
  int bit = 5;
  for (int second = 1; second < n; second++) {
    for (int first = 0; first < second; first++) {
      if (byte >= static_cast<int>(line.size()))
        throw std::runtime_error("truncated graph6 record");
      const int value = static_cast<unsigned char>(line[byte]) - 63;
      if ((value >> bit) & 1)
        graph.addEdge(first, second);
      if (--bit < 0) {
        bit = 5;
        byte++;
      }
    }
  }
  return graph;
}

std::string graph6(const Graph& graph) {
  std::string result(1, static_cast<char>(graph.n + 63));
  int value = 0;
  int bits = 0;
  for (int second = 1; second < graph.n; second++) {
    for (int first = 0; first < second; first++) {
      value = (value << 1) | graph.hasEdge(first, second);
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

std::string canonicalGraph6(const Graph& input) {
  if (input.n > kMaximumOrder)
    throw std::runtime_error("canonical labeling order exceeds fixed workspace");
  const int m = SETWORDSNEEDED(input.n);
  if (m != 1)
    throw std::runtime_error("unexpected nauty word count");

  ::graph graph[kMaximumOrder] = {};
  ::graph canonical[kMaximumOrder] = {};
  int labels[kMaximumOrder] = {};
  int partition[kMaximumOrder] = {};
  int orbits[kMaximumOrder] = {};
  for (const auto& [first, second] : input.edges())
    ADDONEEDGE(graph, first, second, m);

  static DEFAULTOPTIONS_GRAPH(options);
  options.getcanon = TRUE;
  statsblk stats;
  nauty_check(WORDSIZE, m, input.n, NAUTYVERSIONID);
  densenauty(
      graph, labels, partition, orbits, &options, &stats, m, input.n,
      canonical);

  Graph result(input.n);
  for (int first = 0; first < input.n; first++)
    for (int second = first + 1; second < input.n; second++)
      if (ISELEMENT(GRAPHROW(canonical, first, m), second))
        result.addEdge(first, second);
  return graph6(result);
}

bool isCubic(const Graph& graph) {
  for (int vertex = 0; vertex < graph.n; vertex++)
    if (graph.degree(vertex) != 3)
      return false;
  return true;
}

bool isConnected(const Graph& graph) {
  uint64_t reached = 1;
  while (true) {
    uint64_t next = reached;
    for (int vertex = 0; vertex < graph.n; vertex++)
      if ((reached >> vertex) & 1)
        next |= graph.adj[vertex];
    if (next == reached)
      break;
    reached = next;
  }
  return __builtin_popcountll(reached) == graph.n;
}

bool isBiconnected(const Graph& graph) {
  if (!isConnected(graph))
    return false;
  for (int removed = 0; removed < graph.n; removed++) {
    int start = removed == 0 ? 1 : 0;
    uint64_t reached = uint64_t(1) << start;
    while (true) {
      uint64_t next = reached;
      for (int vertex = 0; vertex < graph.n; vertex++)
        if (vertex != removed && ((reached >> vertex) & 1))
          next |= graph.adj[vertex] & ~(uint64_t(1) << removed);
      if (next == reached)
        break;
      reached = next;
    }
    if (__builtin_popcountll(reached) != graph.n - 1)
      return false;
  }
  return true;
}

int girth(const Graph& graph) {
  int result = graph.n + 1;
  for (int source = 0; source < graph.n; source++) {
    std::array<int, kMaximumOrder> distance;
    std::array<int, kMaximumOrder> parent;
    distance.fill(-1);
    parent.fill(-1);
    std::array<int, kMaximumOrder> queue;
    int begin = 0;
    int end = 0;
    distance[source] = 0;
    queue[end++] = source;
    while (begin < end) {
      const int current = queue[begin++];
      for (int next : graph.neighbors(current)) {
        if (distance[next] == -1) {
          distance[next] = distance[current] + 1;
          parent[next] = current;
          queue[end++] = next;
        } else if (parent[current] != next) {
          result = std::min(
              result, distance[current] + distance[next] + 1);
        }
      }
    }
  }
  return result;
}

bool isPlanarGraph(const Graph& graph) {
  return isPlanar(graph.n, graph.edges(), 0);
}

std::vector<std::vector<int>> cycles(const Graph& graph, const int length) {
  std::vector<std::vector<int>> result;
  std::vector<int> path;
  std::function<void(int, int)> search = [&](const int start, const int current) {
    if (static_cast<int>(path.size()) == length) {
      if (graph.hasEdge(current, start) && path[1] < path.back())
        result.push_back(path);
      return;
    }
    for (int next : graph.neighbors(current)) {
      if (next <= start ||
          std::find(path.begin(), path.end(), next) != path.end())
        continue;
      path.push_back(next);
      search(start, next);
      path.pop_back();
    }
  };
  for (int start = 0; start < graph.n; start++) {
    path = {start};
    search(start, start);
  }
  return result;
}

uint64_t vertexMask(const std::vector<int>& vertices) {
  uint64_t result = 0;
  for (int vertex : vertices)
    result |= uint64_t(1) << vertex;
  return result;
}

uint64_t edgeMask(const Graph& graph, const std::vector<int>& cycle) {
  uint64_t result = 0;
  const std::vector<EdgeTy> edges = graph.edges();
  for (int position = 0; position < static_cast<int>(cycle.size()); position++) {
    const EdgeTy target = std::minmax(
        cycle[position], cycle[(position + 1) % cycle.size()]);
    const auto found = std::find(edges.begin(), edges.end(), target);
    if (found == edges.end())
      throw std::runtime_error("cycle edge is absent");
    result |= uint64_t(1) << (found - edges.begin());
  }
  return result;
}

Graph contractCycle(const Graph& graph, const std::vector<int>& cycle) {
  const uint64_t selected = vertexMask(cycle);
  std::vector<int> remap(graph.n, -1);
  int reducedN = 0;
  for (int vertex = 0; vertex < graph.n; vertex++)
    if (((selected >> vertex) & 1) == 0)
      remap[vertex] = reducedN++;
  const int hub = reducedN++;
  Graph reduced(reducedN);
  for (const auto& [first, second] : graph.edges()) {
    const bool firstSelected = (selected >> first) & 1;
    const bool secondSelected = (selected >> second) & 1;
    if (firstSelected && secondSelected)
      continue;
    const int reducedFirst = firstSelected ? hub : remap[first];
    const int reducedSecond = secondSelected ? hub : remap[second];
    if (reducedFirst == reducedSecond || reduced.hasEdge(reducedFirst, reducedSecond))
      throw std::runtime_error("cycle contraction is not simple");
    reduced.addEdge(reducedFirst, reducedSecond);
  }
  return reduced;
}

std::string orderedSixCycleContraction(
    const Graph& graph, const std::vector<int>& cycle) {
  const uint64_t selected = vertexMask(cycle);
  std::vector<int> outside;
  for (int position = 0; position < 6; position++) {
    const int previous = cycle[(position + 5) % 6];
    const int next = cycle[(position + 1) % 6];
    int found = -1;
    for (int neighbor : graph.neighbors(cycle[position])) {
      if (neighbor != previous && neighbor != next) {
        if (found != -1)
          throw std::runtime_error("6-cycle vertex has two outside neighbors");
        found = neighbor;
      }
    }
    if (found == -1 || ((selected >> found) & 1) ||
        std::find(outside.begin(), outside.end(), found) != outside.end())
      throw std::runtime_error("invalid 6-cycle attachments");
    outside.push_back(found);
  }

  std::vector<int> remap(graph.n, -1);
  int reducedN = 0;
  for (int vertex : outside)
    remap[vertex] = reducedN++;
  for (int vertex = 0; vertex < graph.n; vertex++)
    if (((selected >> vertex) & 1) == 0 && remap[vertex] == -1)
      remap[vertex] = reducedN++;
  const int hub = reducedN++;
  Graph reduced(reducedN);
  for (const auto& [first, second] : graph.edges()) {
    const bool firstSelected = (selected >> first) & 1;
    const bool secondSelected = (selected >> second) & 1;
    if (firstSelected && secondSelected)
      continue;
    reduced.addEdge(
        firstSelected ? hub : remap[first],
        secondSelected ? hub : remap[second]);
  }
  if (reduced.n != 21 || reduced.degree(hub) != 6)
    throw std::runtime_error("invalid ordered 6-cycle contraction");
  return graph6(reduced);
}

struct PathReduction {
  Graph graph;
  int left;
  int middle;
  int right;
  std::vector<int> oldToNew;
};

std::optional<PathReduction> reduceFiveCycle(
    const Graph& graph, const std::vector<int>& cycle,
    const int middlePosition) {
  const uint64_t selected = vertexMask(cycle);
  std::vector<int> outside;
  uint64_t outsideMask = 0;
  for (int position = 0; position < 5; position++) {
    const int previous = cycle[(position + 4) % 5];
    const int next = cycle[(position + 1) % 5];
    int found = -1;
    for (int neighbor : graph.neighbors(cycle[position])) {
      if (neighbor == previous || neighbor == next)
        continue;
      if (((selected >> neighbor) & 1) || found != -1)
        return std::nullopt;
      found = neighbor;
    }
    if (found == -1 || ((outsideMask >> found) & 1))
      return std::nullopt;
    outsideMask |= uint64_t(1) << found;
    outside.push_back(found);
  }

  std::vector<int> remap(graph.n, -1);
  int reducedN = 0;
  for (int vertex = 0; vertex < graph.n; vertex++)
    if (((selected >> vertex) & 1) == 0)
      remap[vertex] = reducedN++;
  const int left = reducedN++;
  const int middle = reducedN++;
  const int right = reducedN++;
  Graph reduced(reducedN);
  try {
    for (const auto& [first, second] : graph.edges())
      if (remap[first] != -1 && remap[second] != -1)
        reduced.addEdge(remap[first], remap[second]);
    reduced.addEdge(left, middle);
    reduced.addEdge(middle, right);
    const std::array<int, 5> owners = {
        left, left, middle, right, right};
    for (int offset = -2; offset <= 2; offset++) {
      const int position = (middlePosition + offset + 5) % 5;
      reduced.addEdge(owners[offset + 2], remap[outside[position]]);
    }
  } catch (const std::runtime_error&) {
    return std::nullopt;
  }
  if (!isCubic(reduced) || !isConnected(reduced))
    return std::nullopt;
  return PathReduction{
      std::move(reduced), left, middle, right, std::move(remap)};
}

Graph markStar(const Graph& base, const int center) {
  Graph marked(base.n + 2);
  for (const auto& edge : base.edges())
    marked.addEdge(edge.first, edge.second);
  marked.addEdge(center, base.n);
  marked.addEdge(base.n, base.n + 1);
  return marked;
}

Graph markBroom(
    const Graph& base, const int center,
    const int firstEligibleLeaf, const int secondEligibleLeaf) {
  Graph marked(base.n + 4);
  for (const auto& edge : base.edges())
    marked.addEdge(edge.first, edge.second);
  marked.addEdge(center, base.n);
  marked.addEdge(base.n, base.n + 1);
  marked.addEdge(firstEligibleLeaf, base.n + 2);
  marked.addEdge(secondEligibleLeaf, base.n + 3);
  return marked;
}

std::string minimumCanonicalContraction(
    const Graph& graph, const std::vector<std::vector<int>>& fiveCycles,
    std::vector<std::string>* codes = nullptr) {
  std::string minimum;
  if (codes)
    codes->clear();
  for (const auto& cycle : fiveCycles) {
    const std::string code = canonicalGraph6(contractCycle(graph, cycle));
    if (codes)
      codes->push_back(code);
    if (minimum.empty() || code < minimum)
      minimum = code;
  }
  return minimum;
}

void prepareClaim3Five(const Graph& graph) {
  if (graph.n != 26 || !isCubic(graph) || girth(graph) != 5)
    throw std::runtime_error("claim3-five expects cubic order-26 graphs of girth 5");
  const auto fiveCycles = cycles(graph, 5);
  if (fiveCycles.empty())
    throw std::runtime_error("claim3-five input has no 5-cycle");
  std::cout << minimumCanonicalContraction(graph, fiveCycles) << '\n';
}

void prepareClaim3Six(const Graph& graph, std::ofstream& direct) {
  if (graph.n != 26 || !isCubic(graph) || girth(graph) < 6)
    throw std::runtime_error("claim3-six expects cubic order-26 graphs of girth at least 6");
  const auto sixCycles = cycles(graph, 6);
  std::pair<std::string, std::string> best;
  std::pair<int, int> bestCycles;
  std::vector<int> bestUncovered;
  bool found = false;
  std::vector<std::string> codes;
  std::vector<uint64_t> masks;
  std::vector<uint64_t> closedNeighborhoods;
  for (const auto& cycle : sixCycles) {
    codes.push_back(canonicalGraph6(contractCycle(graph, cycle)));
    masks.push_back(edgeMask(graph, cycle));
    uint64_t closed = vertexMask(cycle);
    for (int vertex : cycle)
      closed |= graph.adj[vertex];
    closedNeighborhoods.push_back(closed);
  }
  for (int first = 0; first < static_cast<int>(sixCycles.size()); first++) {
    for (int second = first + 1; second < static_cast<int>(sixCycles.size()); second++) {
      if ((masks[first] & masks[second]) != 0)
        continue;
      std::pair<std::string, std::string> candidate = {
          codes[first], codes[second]};
      if (candidate.first > candidate.second)
        std::swap(candidate.first, candidate.second);
      const uint64_t uncoveredMask =
          closedNeighborhoods[first] & closedNeighborhoods[second];
      std::vector<int> uncovered;
      for (int vertex = 0; vertex < graph.n; vertex++)
        if ((uncoveredMask >> vertex) & 1)
          uncovered.push_back(vertex);
      if (!found || uncovered.size() < bestUncovered.size() ||
          (uncovered.size() == bestUncovered.size() && candidate < best)) {
        best = std::move(candidate);
        bestCycles = {first, second};
        bestUncovered = std::move(uncovered);
        found = true;
      }
    }
  }
  if (!found) {
    for (int vertex = 0; vertex < graph.n; vertex++)
      direct << canonicalGraph6(markStar(graph, vertex)) << '\n';
    return;
  }
  std::cout << orderedSixCycleContraction(graph, sixCycles[bestCycles.first])
            << '\n'
            << orderedSixCycleContraction(graph, sixCycles[bestCycles.second])
            << '\n';
  for (int vertex : bestUncovered)
    direct << canonicalGraph6(markStar(graph, vertex)) << '\n';
}

void prepareClaim4(const Graph& parent, std::ofstream& direct) {
  if (parent.n != 28 || !isCubic(parent) || girth(parent) != 5)
    throw std::runtime_error("claim4 expects cubic order-28 graphs of girth 5");

  bool categoryTwo = false;
  bool coveredByMinimum = false;
  std::vector<std::string> broomMarkers;
  std::vector<std::string> starMarkers;
  const auto firstCycles = cycles(parent, 5);
  for (const auto& firstCycle : firstCycles) {
    for (int firstMiddle = 0; firstMiddle < 5; firstMiddle++) {
      const auto first = reduceFiveCycle(parent, firstCycle, firstMiddle);
      if (!first || girth(first->graph) < 5)
        continue;
      if (isPlanarGraph(first->graph)) {
        coveredByMinimum = true;
        continue;
      }
      if (!isBiconnected(first->graph))
        continue;

      const auto secondCycles = cycles(first->graph, 5);
      std::vector<std::string> contractionCodes;
      const std::string minimum = minimumCanonicalContraction(
          first->graph, secondCycles, &contractionCodes);
      bool middleOnSecondCycle = false;
      for (int cycleIndex = 0;
           cycleIndex < static_cast<int>(secondCycles.size()); cycleIndex++) {
        const auto& secondCycle = secondCycles[cycleIndex];
        const bool containsMiddle =
            std::find(secondCycle.begin(), secondCycle.end(), first->middle) !=
            secondCycle.end();
        if (!containsMiddle)
          continue;
        middleOnSecondCycle = true;
        if (contractionCodes[cycleIndex] == minimum)
          coveredByMinimum = true;

        const bool containsLeft =
            std::find(secondCycle.begin(), secondCycle.end(), first->left) !=
            secondCycle.end();
        const bool containsRight =
            std::find(secondCycle.begin(), secondCycle.end(), first->right) !=
            secondCycle.end();
        if (containsLeft && containsRight) {
          categoryTwo = true;
          continue;
        }
        if (containsLeft == containsRight)
          continue;

        const int remaining = containsLeft ? first->right : first->left;
        const int middlePosition = static_cast<int>(
            std::find(secondCycle.begin(), secondCycle.end(), first->middle) -
            secondCycle.begin());
        const auto second = reduceFiveCycle(
            first->graph, secondCycle, middlePosition);
        if (!second)
          continue;
        const int mappedRemaining = second->oldToNew[remaining];
        if (mappedRemaining < 0)
          continue;
        const std::vector<int> centerNeighbors =
            second->graph.neighbors(second->middle);
        std::vector<int> expected = {
            second->left, second->right, mappedRemaining};
        std::vector<int> actual = centerNeighbors;
        std::sort(expected.begin(), expected.end());
        std::sort(actual.begin(), actual.end());
        if (expected != actual)
          continue;
        broomMarkers.push_back(canonicalGraph6(markBroom(
            second->graph, second->middle, second->left, second->right)));
      }
      if (!middleOnSecondCycle)
        starMarkers.push_back(
            canonicalGraph6(markStar(first->graph, first->middle)));
    }
  }

  if (categoryTwo || coveredByMinimum)
    return;
  if (!broomMarkers.empty()) {
    std::cout << *std::min_element(
        broomMarkers.begin(), broomMarkers.end()) << '\n';
  } else if (!starMarkers.empty()) {
    std::cout << *std::min_element(
        starMarkers.begin(), starMarkers.end()) << '\n';
  } else {
    direct << graph6(parent) << '\n';
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 3) {
      std::cerr << "usage: prepare_cubic MODE DIRECT_OUTPUT\n"
                << "MODE is claim3-five, claim3-six, or claim4\n";
      return 2;
    }
    const std::string mode = argv[1];
    std::ofstream direct(argv[2], std::ios::app);
    if (!direct)
      throw std::runtime_error("cannot open direct-output file");

    std::string line;
    uint64_t processed = 0;
    while (std::getline(std::cin, line)) {
      if (line.empty())
        continue;
      const Graph graph = parseGraph6(line);
      if (mode == "claim3-five")
        prepareClaim3Five(graph);
      else if (mode == "claim3-six")
        prepareClaim3Six(graph, direct);
      else if (mode == "claim4")
        prepareClaim4(graph, direct);
      else
        throw std::runtime_error("unknown preparation mode");
      processed++;
    }
    std::cerr << "processed " << processed << " graph6 records\n";
  } catch (const std::exception& error) {
    std::cerr << "prepare_cubic: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
