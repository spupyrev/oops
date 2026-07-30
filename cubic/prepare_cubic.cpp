#include "graph_algorithms.h"

#include <nauty.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kMaximumOrder = 62;

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
    throw std::runtime_error("only graph6 records of order at most 62 are supported");
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
    const int start = removed == 0 ? 1 : 0;
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

struct PathReduction {
  Graph graph;
  int left;
  int middle;
  int right;
  int leftFar;
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
  const int leftFarPosition = (middlePosition + 3) % 5;
  return PathReduction{
      std::move(reduced), left, middle, right,
      remap[outside[leftFarPosition]]};
}

std::string minimumCanonicalContraction(
    const Graph& graph, const std::vector<std::vector<int>>& fiveCycles) {
  std::string minimum;
  for (const auto& cycle : fiveCycles) {
    const std::string code = canonicalGraph6(contractCycle(graph, cycle));
    if (minimum.empty() || code < minimum)
      minimum = code;
  }
  return minimum;
}

void prepareClaim3(const Graph& graph) {
  if (graph.n != 26 || !isCubic(graph) || girth(graph) < 5)
    throw std::runtime_error(
        "claim3 expects cubic order-26 graphs of girth at least 5");
  const auto fiveCycles = cycles(graph, 5);
  if (fiveCycles.empty())
    std::cout << graph6(graph) << '\n';
  else
    std::cout << minimumCanonicalContraction(graph, fiveCycles) << '\n';
}

std::optional<std::vector<int>> overlapSequence(
    const std::vector<int>& cycle, const PathReduction& reduction) {
  const auto middleIt = std::find(
      cycle.begin(), cycle.end(), reduction.middle);
  if (middleIt == cycle.end())
    return std::nullopt;
  const int middlePosition =
      static_cast<int>(middleIt - cycle.begin());
  for (int direction : {-1, 1}) {
    std::vector<int> sequence;
    for (int offset = 0; offset < static_cast<int>(cycle.size()); offset++) {
      int position =
          (middlePosition + direction * offset) %
          static_cast<int>(cycle.size());
      if (position < 0)
        position += static_cast<int>(cycle.size());
      sequence.push_back(cycle[position]);
    }
    // Starting at b, the required words are b-a-L-O-M and
    // b-a-L-O-O-M.  The last vertex is forced to be the attachment at b:
    // b's other two neighbors are a and the excluded path vertex c.
    if (sequence[1] == reduction.left &&
        sequence[2] == reduction.leftFar &&
        sequence.back() != reduction.right)
      return sequence;
  }
  return std::nullopt;
}

std::optional<Graph> contractOverlap(
    const Graph& graph, const std::vector<int>& sequence,
    const int right) {
  // Contract the complete overlap, including the path endpoint c, to one
  // hub.  A 5-cycle overlap leaves six boundary edges; a 6-cycle overlap
  // leaves seven.
  const uint64_t selected = vertexMask(sequence);
  std::vector<int> outside;
  for (int position = 0;
       position < static_cast<int>(sequence.size()); position++) {
    const int previous =
        sequence[(position + sequence.size() - 1) % sequence.size()];
    const int next = sequence[(position + 1) % sequence.size()];
    int found = -1;
    for (int neighbor : graph.neighbors(sequence[position])) {
      if (neighbor == previous || neighbor == next)
        continue;
      if (found != -1)
        return std::nullopt;
      found = neighbor;
    }
    if (found == -1)
      return std::nullopt;
    outside.push_back(found);
  }
  if (outside.front() != right)
    return std::nullopt;
  std::vector<int> endpoints(outside.begin() + 1, outside.end());
  for (int neighbor : graph.neighbors(right))
    if (neighbor != sequence.front())
      endpoints.push_back(neighbor);
  std::sort(endpoints.begin(), endpoints.end());
  if (std::unique(endpoints.begin(), endpoints.end()) != endpoints.end())
    return std::nullopt;

  const uint64_t removed = selected | (uint64_t(1) << right);
  for (int endpoint : endpoints)
    if ((removed >> endpoint) & 1)
      return std::nullopt;
  std::vector<int> remap(graph.n, -1);
  int reducedN = 0;
  for (int vertex = 0; vertex < graph.n; vertex++)
    if (((removed >> vertex) & 1) == 0)
      remap[vertex] = reducedN++;
  const int hub = reducedN++;
  Graph quotient(reducedN);
  try {
    for (const auto& [first, second] : graph.edges())
      if (remap[first] != -1 && remap[second] != -1)
        quotient.addEdge(remap[first], remap[second]);
    for (int endpoint : endpoints)
      quotient.addEdge(hub, remap[endpoint]);
  } catch (const std::runtime_error&) {
    return std::nullopt;
  }
  if (!isConnected(quotient) ||
      quotient.degree(hub) != static_cast<int>(sequence.size()) + 1)
    return std::nullopt;
  for (int vertex = 0; vertex < hub; vertex++)
    if (quotient.degree(vertex) != 3)
      return std::nullopt;
  return quotient;
}

char prepareClaim4Joint(const Graph& parent) {
  if (parent.n != 28 || !isCubic(parent) || girth(parent) != 5)
    throw std::runtime_error(
        "claim4-joint expects cubic order-28 graphs of girth 5");

  const auto fiveCycles = cycles(parent, 5);
  if (fiveCycles.empty())
    throw std::runtime_error("claim4-joint input has no 5-cycle");
  std::string bestSixHub;
  std::string bestSevenHub;
  for (const auto& cycle : fiveCycles) {
    for (bool reverse : {false, true}) {
      std::vector<int> oriented = cycle;
      if (reverse)
        std::reverse(oriented.begin() + 1, oriented.end());
      for (int middle = 0; middle < 5; middle++) {
        const auto reduction = reduceFiveCycle(parent, oriented, middle);
        if (!reduction)
          throw std::runtime_error("claim4-joint path reduction failed");

        const auto secondCycles = cycles(reduction->graph, 5);
        for (const auto& secondCycle : secondCycles) {
          const auto middleIt = std::find(
              secondCycle.begin(), secondCycle.end(), reduction->middle);
          if (middleIt != secondCycle.end() &&
              std::find(
                  secondCycle.begin(), secondCycle.end(), reduction->left) !=
                  secondCycle.end() &&
              std::find(
                  secondCycle.begin(), secondCycle.end(), reduction->right) !=
                  secondCycle.end()) {
            const auto order24 = reduceFiveCycle(
                reduction->graph, secondCycle,
                static_cast<int>(middleIt - secondCycle.begin()));
            if (order24 && girth(order24->graph) >= 5 &&
                isBiconnected(order24->graph)) {
              std::cout << "D\n";
              return 'D';
            }
          }
          const auto sequence = overlapSequence(secondCycle, *reduction);
          if (sequence) {
            const auto quotient = contractOverlap(
                reduction->graph, *sequence, reduction->right);
            if (quotient) {
              const std::string code = canonicalGraph6(*quotient);
              if (bestSixHub.empty() || code < bestSixHub)
                bestSixHub = code;
            }
          }
        }
        for (const auto& secondCycle : cycles(reduction->graph, 6)) {
          const auto sequence = overlapSequence(secondCycle, *reduction);
          if (!sequence)
            continue;
          const auto quotient = contractOverlap(
              reduction->graph, *sequence, reduction->right);
          if (quotient) {
            const std::string code = canonicalGraph6(*quotient);
            if (bestSevenHub.empty() || code < bestSevenHub)
              bestSevenHub = code;
          }
        }
      }
    }
  }
  if (!bestSixHub.empty()) {
    std::cout << "B " << bestSixHub << '\n';
    return 'B';
  }
  if (!bestSevenHub.empty()) {
    std::cout << "C " << bestSevenHub << '\n';
    return 'C';
  }
  throw std::runtime_error("claim4-joint found no overlap certificate");
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc != 2) {
      std::cerr << "usage: prepare_cubic MODE\n"
                << "MODE is claim3 or claim4-joint\n";
      return 2;
    }
    const std::string mode = argv[1];
    std::string line;
    uint64_t processed = 0;
    uint64_t jointDirect = 0;
    uint64_t jointExtraCores = 0;
    while (std::getline(std::cin, line)) {
      if (line.empty())
        continue;
      const Graph graph = parseGraph6(line);
      if (mode == "claim3")
        prepareClaim3(graph);
      else if (mode == "claim4-joint") {
        const char classification = prepareClaim4Joint(graph);
        if (classification == 'D')
          jointDirect++;
        else if (classification == 'B' || classification == 'C')
          jointExtraCores++;
        else
          throw std::runtime_error("invalid claim4-joint classification");
      } else
        throw std::runtime_error("unknown preparation mode");
      processed++;
    }
    std::cerr << "processed " << processed << " graph6 records\n";
    if (mode == "claim4-joint")
      std::cerr << "claim4-joint: D=" << jointDirect
                << " hubs=" << jointExtraCores << '\n';
  } catch (const std::exception& error) {
    std::cerr << "prepare_cubic: " << error.what() << '\n';
    return 1;
  }
  return 0;
}
