#pragma once

#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <functional>
#include <string>
#include <map>


class Attr {
 public:
  std::map<std::string, std::string> attr;

 public:
  std::string getAttr(const std::string& key) const {
    if (attr.count(key) == 0) {
      std::cerr << "attribute '" << key << "' not found" << "\n";
    }

    auto it = attr.find(key);
    assert(it != attr.end());
    return it->second;
  }

  std::string getAttr(const std::string& key, const std::string& defaultValue) const {
    if (attr.count(key) == 0) {
      return defaultValue;
    }

    auto it = attr.find(key);
    assert(it != attr.end());
    return it->second;
  }

  double getDoubleAttr(const std::string& key) const {
    auto value = getAttr(key);
    return std::stod(value);
  }

  void setDoubleAttr(const std::string& key, double value) {
    setAttr(key, std::to_string(value));
  }

  void setAttr(const std::string& key, const std::string& value) {
    attr[key] = value;
  }

  void setAttr(const std::map<std::string, std::string>& attr_) {
    attr = attr_;
  }

  bool hasAttr(const std::string& key) const {
    return (attr.find(key) != attr.end());
  }

  void removeAttr(const std::string& key) {
    if (attr.find(key) != attr.end()) {
      attr.erase(attr.find(key));
    }
  }
};

class IONode : public Attr {
 public:
  IONode(const IONode&) = delete;
  IONode(IONode&&) = default;
  IONode& operator=(const IONode&) = delete;
  IONode& operator=(IONode&&) = default;

  int index;
  std::string id;

  IONode(int index, const std::string& id): index(index), id(id) {}

  void setX(double x) {
    setDoubleAttr("x", x);
  }
  void setY(double y) {
    setDoubleAttr("y", y);
  }
};

class IOEdge : public Attr {
 public:
  IOEdge(const IOEdge&) = delete;
  IOEdge(IOEdge&&) = default;
  IOEdge& operator=(const IOEdge&) = delete;
  IOEdge& operator=(IOEdge&&) = default;

  std::string source;
  std::string target;

  IOEdge(const std::string& source, const std::string& target): source(source), target(target) {}
};

class IOGraph {
 public:
  IOGraph(const IOGraph&) = delete;
  IOGraph& operator = (const IOGraph&) = delete;

  IOGraph() {}
  ~IOGraph() {}

  void clear() {
    id2styleIdx.clear();
    id2nodeIdx.clear();
    id2edgeIdx.clear();
    style.clear();
    nodes.clear();
    edges.clear();
  }

  std::vector<IONode> style;
  std::vector<IONode> nodes;
  std::vector<IOEdge> edges;

  std::map<std::string, size_t> id2styleIdx;
  std::map<std::string, size_t> id2nodeIdx;
  std::map<std::pair<std::string, std::string>, size_t> id2edgeIdx;

  IONode* addNode(const std::string& id) {
    assert(getNode(id) == nullptr);
    nodes.emplace_back((int)nodes.size(), id);
    id2nodeIdx[id] = nodes.size() - 1;
    return &nodes.back();
  }

  IOEdge* addEdge(const std::string& source, const std::string& target) {
    assert(getOrCreateNode(source) != nullptr);
    assert(getOrCreateNode(target) != nullptr);
    edges.emplace_back(source, target);
    id2edgeIdx[std::make_pair(source, target)] = edges.size() - 1;
    return &edges.back();
  }

  IONode* getOrCreateNode(const std::string& id) {
    auto it = id2nodeIdx.find(id);
    return it == id2nodeIdx.end() ? addNode(id) : &nodes[it->second];
  }

  IONode* getNode(const std::string& id) {
    auto it = id2nodeIdx.find(id);
    return it == id2nodeIdx.end() ? nullptr : &nodes[it->second];
  }

  IOEdge* getEdge(const std::string& source, const std::string& target) {
    auto pair = std::make_pair(source, target);
    auto it = id2edgeIdx.find(pair);

    if (it == id2edgeIdx.end()) {
      pair = std::make_pair(target, source);
      it = id2edgeIdx.find(pair);
    }

    return it == id2edgeIdx.end() ? nullptr : &edges[it->second];
  }

  IONode* addStyle(const std::string& id) {
    auto it = id2styleIdx.find(id);
    (void)it;
    assert(it == id2styleIdx.end());
    style.emplace_back((int)style.size(), id);
    id2styleIdx[id] = style.size() - 1;
    return &style.back();
  }

  void checkConsistency() const {
    for (size_t i = 0; i < nodes.size(); i++) {
      const auto& node = nodes[i];
      (void)node;
      assert(node.index == (int)i);
      assert(id2nodeIdx.find(node.id)->second == i);
    }

    for (const auto& edge : edges) {
      (void)edge;
      assert(id2nodeIdx.find(edge.source) != id2nodeIdx.end());
      assert(id2nodeIdx.find(edge.target) != id2nodeIdx.end());
    }
  }

  bool empty() const {
    return nodes.empty();
  }

  void extractEdges(int& n, std::vector<std::pair<int, int>>& outEdges) const {
    n = (int)nodes.size();
    for (const auto& edge : edges) {
      const int s = id2nodeIdx.at(edge.source);
      const int t = id2nodeIdx.at(edge.target);
      assert(s != t);
      if (s < t)
        outEdges.push_back({s, t});
      else if (t < s)
        outEdges.push_back({t, s});
    }
  }
};


class GraphParser {
  GraphParser(const GraphParser&);
  GraphParser& operator = (const GraphParser&);

public:
  GraphParser() {}

  bool readGraph(const std::string& filename, IOGraph& graph) const {
    auto func = static_cast<bool (GraphParser::*)(std::istream&, IOGraph&) const>(&GraphParser::readGraph);
    return wrapRead(func, filename, graph);
  }
  
  bool readGraph(std::istream& in, IOGraph& graph) const {
    std::stringstream input;
    input << in.rdbuf();
    input.clear();
    input.seekg(0, std::ios::beg);
    graph.clear();
    bool res = readDotGraph(input, graph);

    if (res) {
      return true;
    }

    input.clear();
    input.seekg(0, std::ios::beg);
    graph.clear();
    res = readGmlGraph(input, graph);

    if (res) {
      return true;
    }

    return false;
  }

  bool readDotGraph(const std::string& filename, IOGraph& graph) const {
    auto func = static_cast<bool (GraphParser::*)(std::istream&, IOGraph&) const>(&GraphParser::readDotGraph);
    return wrapRead(func, filename, graph);
  }
  bool readDotGraph(std::istream& in, IOGraph& graph) const;

  bool readGmlGraph(const std::string& filename, IOGraph& graph) const {
    auto func = static_cast<bool (GraphParser::*)(std::istream&, IOGraph&) const>(&GraphParser::readGmlGraph);
    return wrapRead(func, filename, graph);
  }
  bool readGmlGraph(std::istream& in, IOGraph& graph) const;

  bool readGraphmlGraph(const std::string& filename, IOGraph& graph) const {
    auto func = static_cast<bool (GraphParser::*)(std::istream&, IOGraph&) const>(&GraphParser::readGraphmlGraph);
    return wrapRead(func, filename, graph);
  }
  bool readGraphmlGraph(std::istream& in, IOGraph& graph) const;

  bool writeDotGraph(const std::string& filename, IOGraph& graph) const {
    auto func = static_cast<bool (GraphParser::*)(std::ostream&, IOGraph&) const>(&GraphParser::writeDotGraph);
    return wrapWrite(func, filename, graph);
  }
  bool writeDotGraph(std::ostream& out, IOGraph& graph) const;

  bool writeGmlGraph(const std::string& filename, IOGraph& graph) const {
    auto func = static_cast<bool (GraphParser::*)(std::ostream&, IOGraph&) const>(&GraphParser::writeGmlGraph);
    return wrapWrite(func, filename, graph);
  }
  bool writeGmlGraph(std::ostream& out, IOGraph& graph) const;

  bool writeSvgGraph(const std::string& filename, IOGraph& graph) const {
    auto func = static_cast<bool (GraphParser::*)(std::ostream&, IOGraph&) const>(&GraphParser::writeSvgGraph);
    return wrapWrite(func, filename, graph);
  }
  bool writeSvgGraph(std::ostream& out, IOGraph& graph) const;

private:
  void checkFile(const std::string& filename) const {
    if (filename != "") {
      std::ifstream fileStream;
      fileStream.open(filename.c_str(), std::ios::in);

      if (!fileStream) {
        std::cerr << "input file '" << filename << "' doesn't exist\n";
        throw 20;
      }

      fileStream.close();
    }
  }

  bool wrapRead(std::function<bool(const GraphParser&, std::istream&, IOGraph&)> func, const std::string& filename, IOGraph& graph) const {
    checkFile(filename);

    if (filename != "") {
      std::ifstream fileStream;
      fileStream.open(filename.c_str(), std::ios::in);
      bool result = func(*this, fileStream, graph);
      fileStream.close();
      return result;
    } else {
      return func(*this, std::cin, graph);
    }
  }

  bool wrapWrite(std::function<bool(const GraphParser&, std::ostream&, IOGraph&)> func, const std::string& filename, IOGraph& graph) const {
    if (filename != "") {
      std::ofstream fileStream;
      fileStream.open(filename.c_str(), std::ios::out);
      bool result = func(*this, fileStream, graph);
      fileStream.close();
      return result;
    } else {
      return func(*this, std::cout, graph);
    }
  }
};


using EdgeTy = std::pair<int, int>;
using AdjListTy = std::vector<std::vector<int>>;

class GraphList {
private:  
  GraphList(const GraphList&) = delete;
  GraphList& operator = (const GraphList&) = delete;

public:
  GraphList() = default;
  virtual ~GraphList() = default;

  virtual size_t size() const = 0;

  virtual std::pair<std::string, AdjListTy> next() = 0;
};

class GraphListRaw : public GraphList {
public:
  GraphListRaw(const std::vector<std::pair<std::string, AdjListTy>>& graphs) : graphs(graphs) {}
  
  size_t size() const override {
    return graphs.size();
  }

  std::pair<std::string, AdjListTy> next() override {
    assert(it < graphs.size());
    return graphs[it++];
  }

private:  
  std::vector<std::pair<std::string, AdjListTy>> graphs;
  size_t it = 0;
};

class GraphListG6 : public GraphList {
public:
  GraphListG6(const std::string& filename, 
              const std::string& part,
              const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter);

  ~GraphListG6() {
    ifs.close();
  }

  size_t size() const override {
    return numGraphs;
  }

  std::pair<std::string, AdjListTy> next() override;

private:
  std::pair<int, std::vector<EdgeTy>> parse(const std::string& line);

  size_t numGraphs = 0;
  std::ifstream ifs;
  std::string filename;
  int partIdx;
  int numParts;
  std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter;
  int graphIdx;
};

std::vector<std::pair<std::string, AdjListTy>> readMagma(const std::string& filename, const int maxN);
std::vector<std::pair<std::string, AdjListTy>> readLst(const std::string& filename, const int maxN);

std::vector<std::pair<std::string, AdjListTy>> readCfg(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter);

std::vector<std::pair<std::string, AdjListTy>> readG6(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter);

std::vector<std::pair<std::string, AdjListTy>> readS6(
    const std::string& filename, 
    const std::string& part,
    const std::function<bool(const int, const std::vector<EdgeTy>&)> graphFilter);


struct InputGraph;
struct Result;

void printInput(
    const std::string& filename, 
    const std::string& graphName, 
    const InputGraph& graph, 
    const int verbose);
    
void printInputDot(
    const std::string& filename, 
    const std::string& graphName, 
    const InputGraph& graph, 
    bool append);

void printOutput(
    const std::string& filename, 
    const std::string& graphName, 
    const InputGraph& graph, 
    const Result& result,
    const int verbose);
