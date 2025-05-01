#include "io.h"
#include "common.h"
#include "logging.h"

#include <fstream>
#include <map>
#include <cmath>

using namespace std;

namespace svg_parser {

void writeHead(ostream& out, int minX, int minY, int width, int height) {
  out << "<?xml version='1.0' encoding='utf-8' standalone='yes'?>\n";
  // out << "<svg width='" << width << "' height='" << height << "'\n";
  out << "<svg width='100%' height='100%' viewBox='" << minX << " " << minY << " " << width << " " << height << "'\n";
  out << "	version = '1.1'\n";
  out << "	xmlns = 'http://www.w3.org/2000/svg'\n";
  out << "	xmlns:xlink = 'http://www.w3.org/1999/xlink'\n";
  out << "	xmlns:ev = 'http://www.w3.org/2001/xml-events'>\n";
}

void writeTail(ostream& out) {
  out << "</svg>\n";
}

void drawCircle(ostream& out, 
                double x, double y, double r, 
                int width, int dasharray,
                const string& stroke) {
    out << "	<circle cx='" << x << "' cy='" << y 
        << "' r='" << r 
        << "' style='fill:none"
        << "' stroke-width='" << width
        << "' stroke-dasharray='" << dasharray
        << "' stroke='" << stroke	<< "'/>\n";
}

void drawCircle(ostream& out, 
                double x, double y, double r, 
                const string& fill, 
                const string& stroke) {
    out << "	<circle cx='" << x << "' cy='" << y 
            << "' r='" << r 
            << "' fill='" << fill 
            << "' stroke='" << stroke	<< "'/>\n";
}

void drawLine(ostream& out, 
             double x1, double y1, 
             double x2, double y2, 
             int width, int dasharray,
             const string& stroke, const string& custom) {
  out << "	<path d='" << "M" << x1 << " " << y1 
                                            << "L " << x2 << " " << y2
        << "' stroke-width='" << width
        << "' stroke-dasharray='" << dasharray
        << "' fill='" << "none"
        << "' stroke='" << stroke	<< "'";
  if (custom.length() > 0) out << " " << custom;
  out << "/>\n";
}

void drawBezier(ostream& out, 
                double x1, double y1, 
                double x2, double y2, 
                double x3, double y3, 
                double x4, double y4, 
                int width, const std::string& stroke) {
  out << "	<path d='" << "M " << x1 << " " << y1 
      << " C " 
      << x2 << " " << y2 
      << ", " << x3 << " " << y3 
      << ", " << x4 << " " << y4
      << "' stroke-width='" << width
      << "' fill='" << "none"
      << "' stroke='" << stroke	<< "'"
      << "/>\n";
}

void drawSemiarc(ostream& out, 
                 double x1, double y1, 
                 double x2, double y2, 
                 int width, const std::string& stroke) {
  const double r = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) / 2;
  out << "	<path d='" << "M " << x1 << "," << y1 
      << " A " << r << "," << r / 2 << " 0 0,1" 
      << x2 << "," << y2 
      << "' stroke-width='" << width
      << "' fill='" << "none"
      << "' stroke='" << stroke	<< "'"
      << "/>\n";
}

void drawBiarc(ostream& out, 
               double x1, double y1, 
               double x2, double y2, 
               double x3, double y3,
               bool arc1Up, bool arc2Up, 
               int width, const std::string& stroke) {
  const double r1 = sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)) / 2;
  const double r2 = sqrt((x3 - x2) * (x3 - x2) + (y3 - y2) * (y3 - y2)) / 2;

  out << "	<path d='" << "M " << x1 << "," << y1
      << " A " << r1 << "," << 2 * r1 / 3 << (arc1Up ? " 0 0,1 " : " 0 1,0 ") << x2 << " " << y2
      << " A " << r2 << "," << 2 * r2 / 3 << (arc2Up ? " 0 0,1 " : " 0 1,0 ") << x3 << " " << y3
      << "' stroke-width='" << width
      << "' fill='" << "none"
      << "' stroke='" << stroke	<< "'"
      << "/>\n";
}

void drawText(ostream& out, double x, double y, int fontSize, const string& text) {
  out << "	<text x='" << x << "' y='" << y << "'"
      << "  style='font-family:Verdana;font-size:" << fontSize << "' >" 
      << text << "</text>\n";
}

void defineArrowheads(ostream& out, const map<int, string>& colors) {
  out << "<defs>\n";
  for (auto& it : colors) {
    out << "  <marker id='head" << it.first << "' orient='auto' markerWidth='4' markerHeight='4' refX='8.8' refY='2'>\n";
    out << "    <path d='M0,0 V4 L4,2 Z' fill='" << it.second << "' />\n";
    out << "  </marker>\n";
  }
  out << "</defs>\n";
}

} // namespace svg_parser

using namespace svg_parser;

bool GraphParser::writeSvgGraph(ostream& out, IOGraph& graph) const {
  double minX = graph.nodes[0].getDoubleAttr("x");
  double minY = graph.nodes[0].getDoubleAttr("y");
  double maxX = graph.nodes[0].getDoubleAttr("x");
  double maxY = graph.nodes[0].getDoubleAttr("y");
  for (size_t i = 0; i < graph.nodes.size(); i++) {
    const auto& v = graph.nodes[i];
    minX = std::min(minX, v.getDoubleAttr("x"));
    minY = std::min(minY, v.getDoubleAttr("y"));
    maxX = std::max(maxX, v.getDoubleAttr("x"));
    maxY = std::max(maxY, v.getDoubleAttr("y"));
  }

  const double width = (maxX - minX + 5);
  const double height = (maxY - minY + 5);
  writeHead(out, minX - width * 0.1, minY - height * 0.05, width * 1.2, height * 1.1);

  for (size_t i = 0; i < graph.edges.size(); i++) {
    const auto& e = graph.edges[i];
    CHECK(e.hasAttr("x") && e.hasAttr("y"));

    auto xx = SplitNotNullInt(e.getAttr("x"), "###");
    auto yy = SplitNotNullInt(e.getAttr("y"), "###");
    CHECK(xx.size() == 3 && yy.size() == 2);
    CHECK(yy[0] == 1 || yy[0] == -1);
    CHECK(yy[1] == 1 || yy[1] == -1);
    const bool arc1Up = yy[0] == 1;
    const bool arc2Up = yy[1] == 1;

    // white background
    drawBiarc(out, xx[0], 0, xx[1], 0, xx[2], 0, arc1Up, arc2Up, 6, "#FFFFFF");
    // the edge itself
    drawBiarc(out, xx[0], 0, xx[1], 0, xx[2], 0, arc1Up, arc2Up, e.getDoubleAttr("width"), e.getAttr("fill"));
  }

  for (size_t i = 0; i < graph.nodes.size(); i++) {
    const auto& v = graph.nodes[i];
    std::string label = v.hasAttr("label") ? v.getAttr("label") : v.id;
    double x = v.getDoubleAttr("x");
    double y = v.getDoubleAttr("y");
    double r = v.getDoubleAttr("w");
    std::string color = v.getAttr("fill");

    drawCircle(out, x, y, r, color, "#000000");

    if (label.length() == 1) {
      drawText(out, x - 5, y + 5, 15, label);
    } else if (label.length() == 2) {
      drawText(out, x - 8, y + 5, 13, label);
    } else if (label.length() >= 3) {
      drawText(out, x - 9, y + 4, 9, label);
    }
  }

  writeTail(out);
  return true;
}
