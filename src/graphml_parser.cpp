#include "io.h"
#include "common.h"
#include "logging.h"

#include <fstream>
#include <map>

using namespace std;

namespace graphml_parser {

string getAttr(map<string, string>& attrs, const string& name) {
  if (!attrs.count(name)) {
    std::cerr << "attribute '" << name << "' not found\n";
    throw 120;
  }

  return attrs[name];
}

map<string, string> parseAttrs(const vector<string>& content) {
  map<string, string> attrs;

  for (auto line : content) {
    vector<string> tmp = SplitNotNull(line, "=\"/");

    if (tmp.size() == 2) {
      const string key = tmp[0];
      const string value = tmp[1];
      attrs[key] = value;
    }
  }

  return attrs;
}

void parseNode(IOGraph& graph, const std::vector<std::string>& content) {
  map<string, string> attrs = parseAttrs(content);
  const string id = getAttr(attrs, "id");
  auto node = graph.addNode(id);

  node->setAttr("label", id);
}

void parseEdge(IOGraph& graph, const vector<string>& content) {
  map<string, string> attrs = parseAttrs(content);
  const string source = getAttr(attrs, "source");
  const string target = getAttr(attrs, "target");

  // ignore loops
  if (source == target)
    return;

  graph.addEdge(source, target);
}

bool readGraphmlGraphInt(istream& in, IOGraph& graph) {
  bool in_graph = false;
  string line;

  while (getline(in, line)) {
    vector<string> tmp = SplitNotNull(line, " \t\r<>");

    if (tmp[0] == "graph") {
      CHECK(!in_graph, 120);
      in_graph = true;
    } else if (tmp[0] == "/graph") {
      CHECK(in_graph, 120);
      in_graph = false;
    } else if (tmp[0] == "node") {
      parseNode(graph, tmp);
    } else if (tmp[0] == "edge") {
      parseEdge(graph, tmp);
    }
  }

  CHECK(!in_graph, 120);
  return !graph.nodes.empty();
}

} // namespace graphml_parser

bool GraphParser::readGraphmlGraph(istream& in, IOGraph& graph) const {
  try {
    return graphml_parser::readGraphmlGraphInt(in, graph);
  } catch (int code) {
    LOG("graphml file parsing exception: %d", code);
    return false;
  }
}
