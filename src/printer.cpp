#include "io.h"
#include "logging.h"
#include "one_planar.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <queue>
#include <set>
#include <cmath>

using namespace std;

std::string getEdgeColor(int edgePage) {
  if (edgePage == 0) {
    return "#ca0020";
  } else if (edgePage == 1) {
    return "#0571b0";
  } else if (edgePage == 2) {
    return "#008837";
  } else if (edgePage == 3) {
    return "#a6611a";
  } else if (edgePage == 4) {
    return "#7b3294";
  } else if (edgePage == 5) {
    return "#000000";
  } else if (edgePage == 6) {
    return "#999999";
  } else if (edgePage == 7) {
    return "#00ccff";
  } else if (edgePage == 8) {
    return "#ff00ff";
  } else if (edgePage == 9) {
    return "#00ff00";
  } else if (edgePage == 10) {
    return "#ff9900";
  } else if (edgePage == 11) {
    return "#cc99ff";
  }
  return "#000000";
}

std::string getNodeColor(int nodeType) {
  if (nodeType == 0) {
    return "#ccccff";
  } else if (nodeType == 1) {
    return "#ffcc00";
  } else if (nodeType == 2) {
    return "#00ccff";
  } else if (nodeType == 3) {
    return "#999999";
  } else if (nodeType == 4) {
    return "#99ccff";
  } else if (nodeType == 5) {
    return "#ff99cc";
  } else if (nodeType == 6) {
    return "#cc99ff";
  } else if (nodeType == 7) {
    return "#ff6600";
  } else if (nodeType == 8) {
    return "#c0c0c0";
  }

  return "#000000";
}

void printInputDot(const std::string& filename, const std::string& graphName, const InputGraph& graph, bool append=false) {
  std::ofstream out;
  if (append)
    out.open(filename, std::ios::app);
  else
    out.open(filename);
  out << "graph " << graphName << " {\n";

  for (int i = 0; i < graph.n; i++) {
    out << "  " << i << ";\n";
  }

  for (const auto& [u, v] : graph.edges) {
    out << "  " << u << " -- " << v << ";\n";
  }

  out << "}\n";
  out.close();

  LOG(TextColor::blue, "written TXT graph to '%s'", filename.c_str());  
}

void printInputGml(const std::string& filename, const InputGraph& graph) {
  const int n = graph.n;
  const auto& edges = graph.edges;

  IOGraph ioGraph;

  for (int i = 0; i < n; i++) {
    ioGraph.addNode(std::to_string(i));
    ioGraph.nodes[i].setAttr("w", "20");
    ioGraph.nodes[i].setAttr("h", "20");
    ioGraph.nodes[i].setAttr("fill", getNodeColor(0));
  }

  // TODO: draw on a circle!!!

  // linear order of vertices
  const double X_STEP = 60;
  const double Y_STEP = 30;
  double curX = 0;
  map<int, double> coordX;
  map<int, double> coordY;
  for (int i = 0; i < n; i++) {
    auto& node = ioGraph.nodes[i];

    node.setAttr("x", to_string(curX));
    node.setAttr("y", "0");

    coordX[i] = curX;
    coordY[i] = 0;
    curX += X_STEP;
  }

  // TODO: make it possible to draw on a circle with straight lines!

  // edges are arcs
  for (size_t i = 0; i < edges.size(); i++) {
    auto& edge = edges[i];
    IOEdge* ioEdge = nullptr;
    if (graph.directions.empty() || graph.directions[i]) {
      ioEdge = ioGraph.addEdge(to_string(edge.first), to_string(edge.second));
    } else {
      ioEdge = ioGraph.addEdge(to_string(edge.second), to_string(edge.first));
    }
    CHECK(ioEdge != nullptr);

    int l = edge.first;
    int r = edge.second;

    ioEdge->setAttr("fill", getEdgeColor(0));
    ioEdge->setAttr("width", "2");
    // ioEdge->setAttr("type", "arc");
    // ioEdge->setAttr("arcType", "fixedRatio");
    ioEdge->setAttr("arcHeight", to_string(Y_STEP));
    ioEdge->setAttr("arcRatio", "1");
    ioEdge->setAttr("x", to_string(coordX[l]) + "###" + to_string(int(coordX[l] + coordX[r]) / 2) + "###" + to_string(coordX[r]));
    ioEdge->setAttr("y", to_string(coordY[l]) + "###" + to_string(coordY[l] - Y_STEP) + "###" + to_string(coordY[r]));
  }

  GraphParser parser;
  parser.writeGmlGraph(filename, ioGraph);
  LOG(TextColor::blue, "written GML result to '%s'", filename.c_str());  
}

void printInput(const std::string& filename, const std::string& graphName, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;

  if (verbose >= 3) {
    for (size_t i = 0; i < edges.size(); i++) {
      LOG("edge (%d, %d); div = %d", edges[i].first, edges[i].second, n + i);
    }
  }

  if (verbose >= 3 && filename != "") {
    const std::string file = filename.substr(0, filename.find_last_of("."));
    const std::string extension = filename.substr(filename.find_last_of(".") + 1);
    if (extension == "dot") {
      printInputDot(file + "_in.dot", graphName, graph, /* append */ false);
    } else if (extension == "cfg") {
      printInputDot(file + "_in.cfg", graphName, graph, /* append */ true);
    } else if (extension == "gml") {
      printInputGml(file + "_in.gml", graph);
    } else {
      LOG("unsupported output extension '%s' for pritning input graph", extension.c_str());
    }
  }
}

void printResultTxt(const InputGraph& graph, const Result& result) {
  const auto& edges = graph.edges;

  cout << "  \033[90m" << "order    : " << "\033[0m" << "[";
  for (size_t i = 0; i < result.order.size(); i++) {
    if (result.order[i].empty())
      continue;
    CHECK(result.order[i].size() == 1 || result.order[i].size() == 2);

    if (i != 0)
      cout << " ";
    if (result.order[i].size() == 1)
      cout << result.order[i][0];
    else
      cout << result.order[i][0] << "/" << result.order[i][1];
  }
  cout << "]\n";

  const int numPlanar = (int)(edges.size() - 2 * result.crossings.size());
  cout << "  \033[90m" << "planar edges" << " [" << numPlanar << "]:" << "\033[0m";
  for (size_t i = 0; i < edges.size(); i++) {
    if (result.isCrossed[i])
      continue;
    cout << " (" << edges[i].first << ", " << edges[i].second << ")";
  }
  cout << "\n";

  cout << "  \033[90m" << "crossings   " << " [" << result.crossings.size() << "]:" << "\033[0m";
  for (size_t i = 0; i < result.crossings.size(); i++) {
    const int e1 = result.crossings[i].first;
    const int e2 = result.crossings[i].second;
    cout << "  (" << edges[e1].first << ", " << edges[e1].second << ")--(" << edges[e2].first << ", " << edges[e2].second << ")";
  }
  cout << "\n";
}

void printResultGmlCircle(const std::string& filename, const int n, const std::vector<EdgeTy>& edges, const Result& result) {
  IOGraph ioGraph;

  const double r = 500;
  std::vector<std::pair<double, double>> coords;
  for (int i = 0; i < n; i++) {
    auto node = ioGraph.addNode(std::to_string(i));
    node->setAttr("w", "20");
    node->setAttr("h", "20");

    const double ai = 3.14 * i / (n - 1);
    const double x = -r * cos(ai);
    const double y = r * sin(ai);
    node->setAttr("x", std::to_string(x));
    node->setAttr("y", std::to_string(y));
    coords.push_back({x, y});
  }
  for (size_t i = 0; i < result.crossings.size(); i++) {
    const int e1 = result.crossings[i].first;
    const int e2 = result.crossings[i].second;
    const double x = (coords[edges[e1].first].first + coords[edges[e1].second].first + coords[edges[e2].first].first + coords[edges[e2].second].first) / 4;
    const double y = (coords[edges[e1].first].second + coords[edges[e1].second].second + coords[edges[e2].first].second + coords[edges[e2].second].second) / 4;

    const std::string divId = std::to_string(e1) + "-" + std::to_string(e2);
    auto node = ioGraph.addNode(divId);
    node->setAttr("label", "");
    node->setAttr("w", "10");
    node->setAttr("h", "10");
    node->setAttr("fill", getNodeColor(1));
    node->setX(x);
    node->setY(y);
  }
  // planar edges
  for (size_t i = 0; i < edges.size(); i++) {
    if (result.isCrossed[i])
      continue;
    auto edge = ioGraph.addEdge(to_string(edges[i].first), to_string(edges[i].second));
    edge->setAttr("fill", getEdgeColor(0));
    edge->setAttr("width", "2");
  }
  // crosssed edges
  for (size_t i = 0; i < result.crossings.size(); i++) {
    const int e1 = result.crossings[i].first;
    const int e2 = result.crossings[i].second;
    const std::string divId = std::to_string(e1) + "-" + std::to_string(e2);

    auto edge1 = ioGraph.addEdge(to_string(edges[e1].first), divId);
    edge1->setAttr("fill", getEdgeColor(1));
    auto edge2 = ioGraph.addEdge(to_string(edges[e1].second), divId);
    edge2->setAttr("fill", getEdgeColor(1));
    auto edge3 = ioGraph.addEdge(to_string(edges[e2].first), divId);
    edge3->setAttr("fill", getEdgeColor(1));
    auto edge4 = ioGraph.addEdge(to_string(edges[e2].second), divId);
    edge4->setAttr("fill", getEdgeColor(1));
  }

  GraphParser parser;
  parser.writeGmlGraph(filename, ioGraph);
  LOG(TextColor::blue, "written GML result to '%s'", filename.c_str());  
}

void printResultStackArcs(const InputGraph& graph, const Result& result, IOGraph& ioGraph) {
  CHECK(!result.stack.empty(), "this is implemented for stack-planarity only");

  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  const int numSegments = 2 * m;

  const double STEP_X = 50;
  const double STEP_Y = 50;
  double curX = 0;
  std::vector<int> vIndex(numVertices, -1);
  for (int i = 0; i < (int)result.order.size(); i++) {
    if (result.order[i].empty())
      continue;
    if (result.order[i][0] < n) {
      // regular vertex
      CHECK(result.order[i].size() == 1);
      auto node = ioGraph.addNode(std::to_string(i));
      node->setAttr("label", std::to_string(result.order[i][0]));
      node->setAttr("w", "20");
      node->setAttr("h", "20");
      node->setAttr("fill", getNodeColor(0));
      node->setDoubleAttr("x", curX);
      node->setDoubleAttr("y", 0);
    } else {
      // division vertex
      CHECK(result.order[i].size() <= 2);
      auto node = ioGraph.addNode(std::to_string(i));
      // node->setAttr("label", "");
      node->setAttr("label", std::to_string(result.order[i][0]));
      node->setAttr("w", "10");
      node->setAttr("h", "10");
      node->setAttr("fill", getNodeColor(1));
      node->setDoubleAttr("x", curX);
      node->setDoubleAttr("y", 0);
    }
    curX += STEP_X;
    for (int v : result.order[i]) {
      vIndex[v] = i;
    }
  }

  // edge segments
  for (int seg = 0; seg < numSegments; seg++) {
    const auto [e_first, e_second] = graph.seg2edge(seg);
    CHECK(vIndex[e_first] != -1 && vIndex[e_second] != -1 && vIndex[e_first] != vIndex[e_second]);

    auto edge = ioGraph.addEdge(to_string(vIndex[e_first]), to_string(vIndex[e_second]));
    edge->setAttr("width", "2");
    CHECK(!result.stack.empty());

    edge->setAttr("fill", result.stack[seg] ? getEdgeColor(0) : getEdgeColor(1));
    edge->setAttr("type", "arc");
    edge->setAttr("arcType", "fixedRatio");
    edge->setAttr("arcHeight", to_string(STEP_Y));
    if (result.stack[seg] == (vIndex[e_first] < vIndex[e_second])) {
      edge->setAttr("arcRatio", "1");
    } else {
      edge->setAttr("arcRatio", "-1");
    }
  }
}

void printResultStackBiarcs(const InputGraph& graph, const Result& result, IOGraph& ioGraph) {
  CHECK(!result.stack.empty(), "this is implemented for stack-planarity only");

  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  const double STEP_X = 50;
  const double STEP_Y = 50;
  double curX = 0;
  std::vector<int> vIndex(numVertices, -1);
  for (int i = 0; i < (int)result.order.size(); i++) {
    if (result.order[i].empty())
      continue;
    if (result.order[i][0] < n) {
      // regular vertex
      CHECK(result.order[i].size() == 1);
      auto node = ioGraph.addNode(std::to_string(i));
      node->setAttr("label", std::to_string(result.order[i][0]));
      node->setAttr("w", "12");
      node->setAttr("h", "12");
      node->setAttr("fill", getNodeColor(0));
      node->setDoubleAttr("x", curX);
      node->setDoubleAttr("y", 0);
    } else {
      // division vertex
      CHECK(result.order[i].size() <= 2);
      auto node = ioGraph.addNode(std::to_string(i));
      node->setAttr("label", "");
      node->setAttr("w", "2");
      node->setAttr("h", "2");
      // node->setAttr("label", std::to_string(result.order[i][0]));
      // node->setAttr("w", "8");
      // node->setAttr("h", "8");
      node->setAttr("fill", getNodeColor(1));
      node->setDoubleAttr("x", curX);
      node->setDoubleAttr("y", 0);
    }
    curX += STEP_X;
    for (int v : result.order[i]) {
      vIndex[v] = i;
    }
  }

  // splines
  for (int i = 0; i < (int)edges.size(); i++) {
    const auto& [u, v] = edges[i];
    CHECK(vIndex[u] != -1 && vIndex[v] != -1 && vIndex[u] != vIndex[v]);
    const int d = graph.findDivIndex(u, v);
    // LOG("edge %d -> %d -> %d", u, d, v);
    CHECK(d >= n);
    const int uX = ioGraph.getNode(std::to_string(vIndex[u]))->getDoubleAttr("x");
    const int uY = ioGraph.getNode(std::to_string(vIndex[u]))->getDoubleAttr("y");
    const int vX = ioGraph.getNode(std::to_string(vIndex[v]))->getDoubleAttr("x");
    const int vY = ioGraph.getNode(std::to_string(vIndex[v]))->getDoubleAttr("y");
    const int dX = ioGraph.getNode(std::to_string(vIndex[d]))->getDoubleAttr("x");
    const int dY = ioGraph.getNode(std::to_string(vIndex[d]))->getDoubleAttr("y");

    auto edge = ioGraph.addEdge(to_string(vIndex[u]), to_string(vIndex[v]));
    edge->setAttr("width", "2");
    edge->setAttr("fill", getEdgeColor(1));
    edge->setAttr("type", "biarc");
    edge->setAttr("arcHeight", to_string(STEP_Y));

    // Format:
    //  x: "ux ### dx ### vx ###"
    //  y: "-1/1 ### -1/1" (1 for positive arc, -1 for negative)

    std::string xx;
    xx += std::to_string(uX);
    xx += "###";
    xx += std::to_string(dX);
    xx += "###";
    xx += std::to_string(vX);
    edge->setAttr("x", xx);

    CHECK(uY == 0 && vY == 0 && dY == 0);
    const int seg1 = 2 * i;
    CHECK(graph.seg2edge_v2(seg1) == std::make_pair(v, d));
    const int seg2 = 2 * i + 1;
    CHECK(graph.seg2edge_v2(seg2) == std::make_pair(u, d));
    std::string yy;
    if (result.stack[seg2] == (vIndex[u] < vIndex[d])) {    
      yy += std::to_string(-1);
    } else {
      yy += std::to_string(1);
    }
    yy += "###";
    if (result.stack[seg1] == (vIndex[d] < vIndex[v])) {    
      yy += std::to_string(-1);
    } else {
      yy += std::to_string(1);
    }
    edge->setAttr("y", yy);
  }
}

void printResultGmlStack(const std::string& filename, const InputGraph& graph, const Result& result) {
  IOGraph ioGraph;
  printResultStackArcs(graph, result, ioGraph);

  GraphParser parser;
  parser.writeGmlGraph(filename, ioGraph);
  LOG(TextColor::blue, "written GML result to '%s'", filename.c_str());  
}

void printResultSvgStack(const std::string& filename, const InputGraph& graph, const Result& result) {
  IOGraph ioGraph;
  printResultStackBiarcs(graph, result, ioGraph);

  GraphParser parser;
  parser.writeSvgGraph(filename, ioGraph);
  LOG(TextColor::blue, "written SVG result to '%s'", filename.c_str());  
}

void printOutput(const std::string& filename, const std::string& graphName, const InputGraph& graph, const Result& result, const int verbose) {
  const std::string extension = filename.substr(filename.find_last_of(".") + 1);

  if (verbose >= 2 || extension == "txt") {
    printResultTxt(graph, result);
  }

  if (filename != "") {
    if (extension == "txt") {
      // TODO
      // printResultTxt(filename, graph, result);
    } else if (extension == "dot") {
      // pass
    } else if (extension == "gml") {
      // printResultGmlCircle(filename, n, edges, result);
      printResultGmlStack(filename, graph, result);
    } else if (extension == "svg") {
      printResultSvgStack(filename, graph, result);
    } else {
      ERROR("unsupported output extension " + extension);
    }
  }
}

