#include "io.h"
#include "common.h"
#include "logging.h"

#include <iostream>
#include <fstream>

using namespace std;

bool readLines(istream& in, vector<string>& lines);
bool isSpecial(const string& line);
bool isNode(const string& line);
bool isEdge(const string& line);
IOEdge* parseEdge(IOGraph& graph, const string& line);
IONode* parseNode(IOGraph& graph, const string& line);
string extractId(const string& line);
map<string, string> parseAttr(const string& line);
void extractAttr(const string& line, string& key, string& value);
void splitLine(const string& line, string& beforeBrakets, string& insideBrackets);
vector<string> splitAttr(const string& line, char separator);
string trim(const string& line);

bool readDotGraphInt(istream& in, IOGraph& graph) {
  vector<string> lines;

  if (!readLines(in, lines)) {
    return false;
  }

  for (size_t i = 0; i < lines.size(); i++) {
    if (isSpecial(lines[i])) {
      parseNode(graph, lines[i]);
    } else if (isEdge(lines[i])) {
      parseEdge(graph, lines[i]);
    } else if (isNode(lines[i])) {
      parseNode(graph, lines[i]);
    } else {
      cerr << "Unknown entry: " << lines[i] << "\n";
      return false;
    }
  }

  // make sure all nodes have labels
  for (auto& node : graph.nodes) {
    if (!node.hasAttr("label")) {
      node.setAttr("label", node.id);
    }
  }

  return true;
}

// Parse input *.dot file
// Throws 110 if parsing is not successful
bool GraphParser::readDotGraph(istream& in, IOGraph& graph) const {
  try {
    return readDotGraphInt(in, graph);
  } catch (int code) {
    LOG("dot file parsing exception: %d", code);
    return false;
  }
}

bool readLines(istream& in, vector<string>& res) {
  char c;
  string s = "";
  int state = -1;
  bool insideQuote = false;

  while (in.get(c)) {
    if (c == '{') {
      state = 0;
      continue;
    }

    if (c == '}') {
      state = -1;
      continue;
    }

    if (state == 0) {
      if (c == ';' && !insideQuote) {
        if (s.length() > 0) {
          res.push_back(s);
        }

        s = "";
      } else {
        if (c == '"') {
          insideQuote = !insideQuote;
        }

        s += c;
      }
    }
  }

  return (state == -1 && res.size() > 0);
}

bool isSpecial(const string& line) {
  string beforeBrakets, insideBrackets;
  splitLine(line, beforeBrakets, insideBrackets);
  string s = trim(beforeBrakets);

  if (s == "node") {
    return true;
  }

  if (s == "graph") {
    return true;
  }

  if (s == "edge") {
    return true;
  }

  return false;
}

bool isNode(const string& line) {
  string beforeBrakets, insideBrackets;
  splitLine(line, beforeBrakets, insideBrackets);

  for (int j = 0; j + 1 < (int)beforeBrakets.length(); j++)
    if (beforeBrakets[j] == '-' && beforeBrakets[j + 1] == '-') {
      return false;
    }

  return true;
}

bool isEdge(const string& line) {
  string beforeBrakets, insideBrackets;
  splitLine(line, beforeBrakets, insideBrackets);

  for (int j = 0; j + 1 < (int)beforeBrakets.length(); j++)
    if (beforeBrakets[j] == '-' && beforeBrakets[j + 1] == '-') {
      return true;
    }

  return false;
}

string extractId(const string& line) {
  string s = trim(line);

  if (s[0] == '"') {
    CHECK(s[s.length() - 1] == '"', 110);
    return s.substr(1, s.length() - 2);
  }

  return s;
}

IOEdge* parseEdge(IOGraph& graph, const string& line) {
  string beforeBrakets, insideBrackets;
  splitLine(line, beforeBrakets, insideBrackets);
  IOEdge* e = nullptr;

  //ids
  for (size_t j = 0; j + 1 < beforeBrakets.length(); j++) {
    if (beforeBrakets[j] == '-' && beforeBrakets[j + 1] == '-') {
      string n1 = beforeBrakets.substr(0, j);
      string n2 = beforeBrakets.substr(j + 2, beforeBrakets.length() - j - 2);

      if (n1 == n2) {
        return nullptr;
      }

      e = graph.addEdge(extractId(n1), extractId(n2));
      break;
    }
  }

  assert(e != nullptr);
  //attrs
  e->setAttr(parseAttr(insideBrackets));
  return e;
}

IONode* parseNode(IOGraph& graph, const string& line) {
  string beforeBrakets, insideBrackets;
  splitLine(line, beforeBrakets, insideBrackets);

  string id = extractId(beforeBrakets);
  auto node = graph.addNode(id);

  auto attrs = parseAttr(insideBrackets);
  node->setAttr(attrs);
  if (attrs.count("label") == 0) {
    node->setAttr("label", id);
  }
  return node;
}

map<string, string> parseAttr(const string& line) {
  string s = trim(line);
  vector<string> tmp = splitAttr(s, ',');
  map<string, string> attr;

  for (int i = 0; i < (int)tmp.size(); i++) {
    string key, value;
    extractAttr(tmp[i], key, value);
    attr[key] = value;
  }

  return attr;
}

void extractAttr(const string& line, string& key, string& value) {
  vector<string> tmp = splitAttr(line, '=');
  CHECK(tmp.size() == 2, 110);
  key = trim(tmp[0]);
  value = extractId(tmp[1]);
}

void splitLine(const string& line, string& beforeBrakets, string& insideBrackets) {
  beforeBrakets = insideBrackets = "";
  int i = 0;
  bool insideQuote = false;

  while (i < (int)line.length()) {
    if (line[i] == '[' && !insideQuote) {
      break;
    }

    if (line[i] == '"') {
      insideQuote = !insideQuote;
    }

    i++;
  }

  beforeBrakets = line.substr(0, i);
  CHECK(!insideQuote, 110);
  int j = i;

  while (j < (int)line.length()) {
    if (line[j] == ']' && !insideQuote) {
      break;
    }

    if (line[j] == '"') {
      insideQuote = !insideQuote;
    }

    j++;
  }

  if (j - i - 1 > 0) {
    insideBrackets = line.substr(i + 1, j - i - 1);
  }
}

vector<string> splitAttr(const string& line, char separator) {
  vector<string> res;
  bool insideQuotes = false;
  string s = "";

  for (size_t i = 0; i < line.length(); i++) {
    if (!insideQuotes) {
      if (line[i] == '"') {
        insideQuotes = true;
        s += line[i];
        continue;
      }

      if (line[i] == separator) {
        if (s.length() > 0) {
          res.push_back(s);
        }

        s = "";
        continue;
      }

      s += line[i];
    } else {
      if (line[i] == '"') {
        insideQuotes = false;
      }

      s += line[i];
    }
  }

  if (s.length() > 0) {
    res.push_back(s);
  }

  return res;
}

string trim(const string& line) {
  if (line.length() == 0) {
    return line;
  }

  int i = 0;

  while (i < (int)line.length()) {
    if (line[i] == ' ' || line[i] == '\n' || line[i] == '\t' || line[i] == '\r') {
      i++;
      continue;
    }

    break;
  }

  int j = (int)line.length() - 1;

  while (j >= 0) {
    if (line[j] == ' ' || line[j] == '\n' || line[j] == '\t' || line[j] == '\r') {
      j--;
      continue;
    }

    break;
  }

  CHECK(i <= j, 110);
  return line.substr(i, j - i + 1);
}

void writeStyles(ostream& out, IOGraph& g);
void writeNodes(ostream& out, IOGraph& g);
void writeStyle(ostream& out, const IONode& n, bool useQ);
void writeNode(ostream& out, const IONode& n, bool useQ);
void writeAttr(ostream& out, const map<string, string>& attr);
void writeEdges(ostream& out, IOGraph& g);
void writeEdge(ostream& out, const IOEdge& n);

bool GraphParser::writeDotGraph(ostream& out, IOGraph& graph) const {
  out << "graph {\n";
  writeStyles(out, graph);
  writeNodes(out, graph);
  writeEdges(out, graph);
  out << "}\n";
  return true;
}

void writeStyles(ostream& out, IOGraph& g) {
  for (size_t i = 0; i < g.style.size(); i++) {
    writeStyle(out, g.style[i], false);
  }
}

void writeNodes(ostream& out, IOGraph& g) {
  for (size_t i = 0; i < g.nodes.size(); i++) {
    writeNode(out, g.nodes[i], true);
  }
}

void writeStyle(ostream& out, const IONode& n, bool useQ) {
  if (useQ) {
    out << "  \"" << n.id << "\" ";
  } else {
    out << "  " << n.id << " ";
  }

  writeAttr(out, n.attr);
  out << ";\n";
}

void writeNode(ostream& out, const IONode& n, bool useQ) {
  if (useQ) {
    out << "  \"" << n.id << "\" ";
  } else {
    out << "  " << n.id << " ";
  }

  if (n.attr.size() > 0) {
    writeAttr(out, n.attr);
  }

  out << ";\n";
}

void writeAttr(ostream& out, const map<string, string>& attr) {
  out << "[";

  for (auto iter = attr.begin(); iter != attr.end(); iter++) {
    string key = (*iter).first;
    string value = (*iter).second;

    if (iter != attr.begin()) {
      out << ", ";
    }

    out << key << "=\"" << value << "\"";
  }

  out << "]";
}

void writeEdges(ostream& out, IOGraph& g) {
  for (size_t i = 0; i < g.edges.size(); i++) {
    writeEdge(out, g.edges[i]);
  }
}

void writeEdge(ostream& out, const IOEdge& e) {
  out << "  \"" << e.source << "\" -- \"" << e.target << "\" ";

  if (e.attr.size() > 0) {
    writeAttr(out, e.attr);
  }

  out << ";\n";
}
