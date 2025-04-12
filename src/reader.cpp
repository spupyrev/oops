#include "io_utils.h"

#include "logging.h"
#include "io_graph.h"
#include "graph_parser.h"
#include "one_planar.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <queue>
#include <set>
#include <cmath>

using namespace std;

std::vector<std::pair<std::string, AdjListTy>> readMagma(const std::string& filename, const int maxN) {
  LOG("readMagma from '%s' with |V| <= %d", filename.c_str(), maxN);

  ifstream ifs(filename);
  CHECK(ifs.good(), "file '%s' not found", filename.c_str());
  std::vector<std::pair<std::string, AdjListTy>> graphs;

  char c;
  int state = -1;
  vector<int> curNumbers;
  int curNumber = 0;

  while (ifs.get(c)) {
    if (c == '<') {
      state = 0;
      continue;
    }

    if (c == '>') {
      state = -1;
      // add a graph
      CHECK(curNumbers.size() % 2 == 1);
      int n = curNumbers[0];
      auto adj = vector<vector<int>>(n, vector<int>());
      for (int i = 1; i < (int)curNumbers.size(); i += 2) {
        int u = curNumbers[i];
        int v = curNumbers[i + 1];
        CHECK(1 <= u && u <= n);
        CHECK(1 <= v && v <= n);
        adj[u - 1].push_back(v - 1);
        adj[v - 1].push_back(u - 1);
      }
      curNumbers.clear();
      if ((int)adj.size() <= maxN)
        graphs.push_back({std::to_string(graphs.size()), adj});
      continue;
    }

    if (state == 0) {
      if (isdigit(c)) {
        curNumber = curNumber * 10 + int(c - '0');
      } else {
        if (curNumber > 0) {
          curNumbers.push_back(curNumber);
        }
        curNumber = 0;
      }
    }
  }

  CHECK(state == -1);
  ifs.close();
  return graphs;
}

std::vector<std::pair<std::string, AdjListTy>> readLst(const string& filename, const int maxN) {
  ifstream ifs(filename);
  CHECK(ifs.good(), "file '%s' not found", filename.c_str());
  std::vector<std::pair<std::string, AdjListTy>> graphs;

  int curN = 0;
  auto adj = vector<vector<int>>();

  std::string line;
  while (std::getline(ifs, line)) {
    auto tmp = SplitNotNull(line, " ");
    if (tmp.empty()) {
      if (curN > 0) {
        if ((int)adj.size() <= maxN)
          graphs.push_back({std::to_string(graphs.size()), adj});
        curN = 0;
        adj.clear();        
      }
      continue;
    }

    CHECK(tmp[0][tmp[0].length() - 1] == ':');
    int v = to_int(tmp[0].substr(0, tmp[0].length() - 1));
    CHECK(v == curN + 1);
    vector<int> neighbors;
    for (int i = 1; i < (int)tmp.size(); i++) {
      int u = to_int(tmp[i]);
      CHECK(u >= 1);
      neighbors.push_back(u - 1);
    }
    adj.push_back(neighbors);
    curN++;
  }

  if (curN > 0) {
    if ((int)adj.size() <= maxN)
      graphs.push_back({std::to_string(graphs.size()), adj});
    curN = 0;
    adj.clear();        
  }

  ifs.close();
  return graphs;
}

std::pair<int, int> parsePartOption(const std::string& part) {
  if (part.length() == 0)
    return {0, 1};
  auto parsed_part = SplitNotNullInt(part, "/");
  CHECK(parsed_part.size() == 1 || parsed_part.size() == 2);
  if (parsed_part.size() == 1) {
    return {parsed_part[0], 123456789};
  }
  CHECK(parsed_part[1] >= 1);
  CHECK(0 <= parsed_part[0] && parsed_part[0] < parsed_part[1]);
  return {parsed_part[0], parsed_part[1]};
}

std::vector<std::pair<std::string, AdjListTy>> readCfg(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter) {
  std::ifstream ifs(filename);
  CHECK(ifs.good(), "file '%s' not found", filename.c_str());

  const auto [partIdx, numParts] = parsePartOption(part);

  std::vector<std::pair<std::string, AdjListTy>> graphs;
  int graphIdx = 0;
  int graphMaxN = 0;
  std::string header;
  while (std::getline(ifs, header)) {
    if (header.find("graph ") != string::npos) {
      vector<string> headerParams = SplitNotNull(header, " ");
      CHECK(headerParams.at(0) == "graph" || headerParams.at(0) == "digraph");
      std::string name = headerParams.at(1);
      name.erase(std::remove(name.begin(), name.end(), '\"'), name.end());

      // starting a new function
      std::vector<EdgeTy> curEdges;
      std::unordered_map<std::string, size_t> b2i;

      std::string line;
      while (std::getline(ifs, line)) {
        if (line == "}")
          break;
        if (line.length() == 0 || line.find(";") == string::npos)
          continue;

        if (line.find("->") == string::npos && line.find("--") == string::npos) {
          // node
          vector<string> params = SplitNotNull(line, " ;");
          if (params.at(0) == "node" || params.at(0) == "edge")
            continue;
          std::string nodeId = params.at(0);
          CHECK(b2i.count(nodeId) == 0, "duplicate node id");
          b2i[nodeId] = b2i.size();
        } else {
          // edge
          std::vector<std::string> params = SplitNotNull(line, " ;");
          CHECK(b2i.count(params.at(0)) == 1, "missing a source node for edge; line = '%s'", line.c_str());
          CHECK(b2i.count(params.at(2)) == 1, "missing a target node for edge");
          const int srcId = std::min(b2i[params[0]], b2i[params[2]]);
          const int dstId = std::max(b2i[params[0]], b2i[params[2]]);
          if (srcId != dstId) {
            // make it acyclic
            if (srcId < dstId)
              curEdges.push_back({srcId, dstId});
            else
              curEdges.push_back({dstId, srcId});
          }
        }
      }

      // skip empty
      if (curEdges.empty()) {
        graphIdx++;
        continue;
      }

      if (graphIdx % numParts != partIdx) {
        graphIdx++;
        continue;
      }

      sort_unique(curEdges);
      const AdjListTy& adj = edges_to_adj(curEdges);
      const int curN = (int)adj.size();
      if (graphFilter(curN, curEdges)) {
        graphs.push_back({name, adj});
        graphMaxN = std::max(graphMaxN, curN);
      }
      graphIdx++;
    }
  }

  ifs.close();

  LOG(TextColor::blue, "extracted %'d graphs (%.2lf%% out of %'d) with maximum |V| = %d from %s", 
    graphs.size(), 100.0 * graphs.size() / graphIdx, graphIdx,
    graphMaxN, filename.c_str()
  );

  return graphs;
}

// TODO: drop this
std::vector<std::pair<std::string, AdjListTy>> readG6(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter) {
  ifstream ifs(filename);
  CHECK(ifs.good(), "file '%s' not found", filename.c_str());

  const auto [partIdx, numParts] = parsePartOption(part);

  int maxN = 0;
  std::vector<std::pair<std::string, AdjListTy>> graphs;
  std::string line;
  int graphIdx = 0;
  while (std::getline(ifs, line)) {
    //std::cout << line << "\n";
    if (graphIdx % numParts != partIdx) {
      graphIdx++;
      continue;
    }
    if (starts_with(line, ">>graph6<<")) {
      line = line.substr(10);
    }

    std::vector<int> data;
    for (char c : line) {
      CHECK(63 <= c && c <= 126);
      data.push_back(int(c - 63));
    }

    int n;
    if (data[0] <= 62) {
      n = data[0];
      data.erase(data.begin(), data.begin() + 1);
    } else if (data[1] <= 62) {
      n = (data[1]<<12) + (data[2]<<6) + data[3];
      data.erase(data.begin(), data.begin() + 4);
    } else {
      ERROR("not implemented yet");
      // return ((data[2]<<30) + (data[3]<<24) + (data[4]<<18) +
      //       (data[5]<<12) + (data[6]<<6) + data[7], data[8:])
    }
    //LOG("n = %d", n);
    CHECK(n >= 1);

    const int nd = (n * (n - 1) / 2 + 5) / 6;
    CHECK((int)data.size() == nd);
    std::vector<bool> bitAdj;
    for (int i = 0; i < nd; i++) {
      for (int j = 0; j < 6; j++) {
        const int bit = 5-j;
        if ((data[i] & (1 << bit)))
          bitAdj.push_back(1);
        else
          bitAdj.push_back(0);
      }
    }
    CHECK((int)bitAdj.size() == nd * 6);

    // the order is (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),...,(n-2,n-1)
    size_t idx = 0;
    std::vector<EdgeTy> edges;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < i; j++) {
        CHECK(idx < bitAdj.size());
        if (bitAdj[idx])
          edges.push_back({j, i});
        idx++;
      }
    }

    if (graphFilter(n, edges)) {
      const std::string base = filename.substr(filename.find_last_of("/") + 1);
      const std::string graphName = base + "_" + to_string(graphIdx);
      graphs.push_back({
        graphName,
        edges_to_adj(n, edges)
      });
      maxN = std::max(n, maxN);
    }
    graphIdx++;
  }

  ifs.close();

  LOG(TextColor::blue, "extracted %'d graphs (%.2lf%% out of %'d) with maximum |V| = %d from %s", 
    graphs.size(), 100.0 * graphs.size() / graphIdx, graphIdx, maxN, filename.c_str()
  );

  return graphs;
}

std::vector<std::pair<std::string, AdjListTy>> readS6(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter) {
  ifstream ifs(filename);
  CHECK(ifs.good(), "file '%s' not found", filename.c_str());

  const auto [partIdx, numParts] = parsePartOption(part);

  int maxN = 0;
  std::vector<std::pair<std::string, AdjListTy>> graphs;
  std::string line;
  int graphIdx = 0;
  while (std::getline(ifs, line)) {
    // std::cout << line << "\n";
    if (graphIdx % numParts != partIdx) {
      graphIdx++;
      continue;
    }

    if (starts_with(line, ">>sparse6<<")) {
      line = line.substr(11);
    }

    CHECK(line[0] == ':');
    line = line.substr(1);
    std::vector<int> data;
    for (char c : line) {
      CHECK(63 <= c && c <= 126);
      data.push_back(int(c - 63));
    }

    int n;
    if (data[0] <= 62) {
      n = data[0];
      data.erase(data.begin(), data.begin() + 1);
    } else if (data[1] <= 62) {
      n = (data[1]<<12) + (data[2]<<6) + data[3];
      data.erase(data.begin(), data.begin() + 4);
    } else {
      ERROR("not implemented yet");
      // return ((data[2]<<30) + (data[3]<<24) + (data[4]<<18) +
      //       (data[5]<<12) + (data[6]<<6) + data[7], data[8:])
    }
    // std::cout << "n = " << n << "\n";
    CHECK(n >= 1);

    int k = 1;
    while ((1 << k) < n)
      k++;
	
    // extract bits
    // CHECK((data.size() * 8) % 6 == 0, "data.size() = %d", data.size());
    std::vector<bool> bits;
    for (int i = 0; i < (int)data.size(); i++) {
      for (int j = 0; j < 6; j++) {
        const int bit = 5-j;
        if ((data[i] & (1 << bit)))
          bits.push_back(1);
        else
          bits.push_back(0);
      }
    }

    std::vector<EdgeTy> edges;
    int v = 0;
    size_t idx = 0;
    while (idx < bits.size()) {
      int b = bits[idx]; idx++;
      // read k bits
      int x = 0;
      for (int i = 0; i < k; i++) {
        x = x * 2 + bits[idx];
        idx++;
      }

      if (b == 1) {
        v += 1;
      }

      if (x >= n || v >= n) {
        break;
      } else if (x > v) {
        v = x;
      } else {
        CHECK(x < v);
        edges.push_back({x, v});
      }
    }

    if (graphFilter(n, edges)) {
      const std::string base = filename.substr(filename.find_last_of("/") + 1);
      const std::string graphName = base + "_" + to_string(graphIdx);
      graphs.push_back({
        graphName,
        edges_to_adj(n, edges)
      });
      maxN = std::max(n, maxN);
    }
    graphIdx++;
  }

  ifs.close();

  LOG(TextColor::blue, "extracted %'d graphs (%.2lf%% out of %'d) with maximum |V| = %d from %s", 
    graphs.size(), 100.0 * graphs.size() / graphIdx, graphIdx, maxN, filename.c_str()
  );

  return graphs;
}

GraphListG6::GraphListG6(
    const std::string& filename_, 
    const std::string& part_,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter_) {
  filename = filename_;
  graphFilter = graphFilter_;
  const auto [partIdx_, numParts_] = parsePartOption(part_);
  partIdx = partIdx_;
  numParts = numParts_;
  graphIdx = 0;

  ifs.open(filename);
  CHECK(ifs.good(), "file '%s' not found", filename.c_str());

  int maxN = 0;
  std::string line;
  while (std::getline(ifs, line)) {
    if (graphIdx % numParts != partIdx) {
      graphIdx++;
      continue;
    }

    const auto& [n, edges] = parse(line);
    if (!graphFilter(n, edges)) {
      graphIdx++;
      std::cerr << "skipped graph filter ::read\n";
      continue;
    }

    graphIdx++;
    numGraphs++;
    maxN = std::max(n, maxN);
  }

  LOG(TextColor::blue, "extracted %'d graphs (%.2lf%% out of %'d) with maximum |V| = %d from %s", 
    numGraphs, 100.0 * numGraphs / graphIdx, graphIdx, maxN, filename.c_str()
  );

  ifs.close();
}

std::pair<std::string, AdjListTy> GraphListG6::next() {
  if (!ifs.is_open()) {
    ifs.open(filename);
    graphIdx = 0;
  }
  CHECK(ifs.good() && ifs.is_open());

  std::string line;
  while (std::getline(ifs, line)) {
    if (graphIdx % numParts != partIdx) {
      graphIdx++;
      continue;
    }

    const auto& [n, edges] = parse(line);
    if (!graphFilter(n, edges)) {
      graphIdx++;
      continue;
    }

    const std::string base = filename.substr(filename.find_last_of("/") + 1);
    const std::string graphName = base + "_" + to_string(graphIdx);
    graphIdx++;
    return std::make_pair(graphName, edges_to_adj(n, edges));
  }
  ERROR("not reachable");
  return {};
}

std::pair<int, std::vector<EdgeTy>> GraphListG6::parse(const std::string& line_) {
  std::string line = line_;
  if (starts_with(line, ">>graph6<<")) {
    line = line.substr(10);
  }

  std::vector<int> data;
  for (char c : line) {
    CHECK(63 <= c && c <= 126);
    data.push_back(int(c - 63));
  }

  int n;
  if (data[0] <= 62) {
    n = data[0];
    data.erase(data.begin(), data.begin() + 1);
  } else if (data[1] <= 62) {
    n = (data[1]<<12) + (data[2]<<6) + data[3];
    data.erase(data.begin(), data.begin() + 4);
  } else {
    ERROR("not implemented yet");
    // return ((data[2]<<30) + (data[3]<<24) + (data[4]<<18) +
    //       (data[5]<<12) + (data[6]<<6) + data[7], data[8:])
  }
  //LOG("parsed n = %d", n);
  CHECK(n >= 1);

  const int nd = (n * (n - 1) / 2 + 5) / 6;
  CHECK((int)data.size() == nd);
  std::vector<bool> bitAdj;
  for (int i = 0; i < nd; i++) {
    for (int j = 0; j < 6; j++) {
      const int bit = 5-j;
      if ((data[i] & (1 << bit)))
        bitAdj.push_back(1);
      else
        bitAdj.push_back(0);
    }
  }
  CHECK((int)bitAdj.size() == nd * 6);

  // the order is (0,1),(0,2),(1,2),(0,3),(1,3),(2,3),...,(n-2,n-1)
  size_t idx = 0;
  std::vector<EdgeTy> edges;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      CHECK(idx < bitAdj.size());
      if (bitAdj[idx])
        edges.push_back({j, i});
      idx++;
    }
  }
  return std::make_pair(n, edges);
}