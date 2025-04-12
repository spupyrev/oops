#pragma once

#include "planar_graph.h"

using EdgeTy = std::pair<int, int>;

// non-planar graphs
void genByClass(std::vector<EdgeTy>& edges, int n, const std::string& graphClass);
void genPath(std::vector<EdgeTy>& edges, int n);
void genCycle(std::vector<EdgeTy>& edges, int n);
void genStar(std::vector<EdgeTy>& edges, int n);
void genBinaryTree(std::vector<EdgeTy>& edges, int n);
void genTernaryTree(std::vector<EdgeTy>& edges, int n);
void genRegularTree(std::vector<EdgeTy>& edges, int n, int degree);
void genApexTree(std::vector<EdgeTy>& edges, int n, int degree);
void genTree(std::vector<EdgeTy>& edges, int n);
void genXTree(std::vector<EdgeTy>& edges, int n);
void genGrid(std::vector<EdgeTy>& edges, int n, bool addHexEdges);
void genComplete(std::vector<EdgeTy>& edges, int n);
void genCompleteSubdivision(const int N, std::vector<EdgeTy>& edges, int& n);
void genCompleteBipartite(std::vector<EdgeTy>& edges, int n, int m);
void genQueue(std::vector<EdgeTy>& edges, int n, int Q, bool bipartite);
void genStack(std::vector<EdgeTy>& edges, int n, int S, bool bipartite);
void genFullStack(std::vector<EdgeTy>& edges, int n);
void genMergedOuterplanar(std::vector<EdgeTy>& edges, int n, int S);
void genMixed(std::vector<EdgeTy>& edges, int n, int S, int Q, bool bipartite);
void genHypercube(std::vector<EdgeTy>& edges, int d, int& n);
void genRegularBibaprtite(std::vector<EdgeTy>& edges, int n, int d);
void genRegular(std::vector<EdgeTy>& edges, int n, int d);
void genRegularExpander(std::vector<EdgeTy>& edges, int n, int d);
void genCaterpillar(std::vector<EdgeTy>& edges, int n, int legs);
void genPathwidth(std::vector<EdgeTy>& edges, int n, int p);
void genOuterKPlanar(std::vector<EdgeTy>& edges, int n, int k, bool bipartite);
void genOptimal1Planar(std::vector<EdgeTy>& edges, int& n);
void genSubdivision(std::vector<EdgeTy>& edges, int k, int& n);
void genCirculant(std::vector<EdgeTy>& edges, int n, const vector<int>& S);
void genQueuesHard(std::vector<EdgeTy>& edges, int n);
void genQueuesHard2(std::vector<EdgeTy>& edges, int n, int verbose);
void genUpwardPathwidth3(std::vector<EdgeTy>& edges, int &n, int verbose);
void genKTreeHard(const int K, std::vector<EdgeTy>& edges, int& n);

enum class ProductType {
  strong,
  direct,
  cartesian,
  blowup
};

void genProduct(std::vector<EdgeTy>& edges, const std::string& classA, int nA, const std::string& classB, int nB, ProductType ptype);
void genProduct(std::vector<EdgeTy>& edges, int nA, const vector<pair<int, int>>& edgesA, int nB, const vector<pair<int, int>>& edgesB, ProductType ptype);
void genMixedProduct(std::vector<EdgeTy>& edges, int& n, int nA, int nB);
void genCliqueSum(std::vector<EdgeTy>& edges, int nA, const vector<pair<int, int>>& edgesA, int nB, const vector<pair<int, int>>& edgesB);

void reduceMaxDegree(int n, std::vector<EdgeTy>& edges, int maxDegree);

enum class StellationType {
  regular,
  ortho
};

// TODO: move to graph
void genIOGraph(PlanarGraph& g);
Vertex* stellateFace(PlanarGraph& g, const Face& f);
void stellateFaces(PlanarGraph& g);
void stellateFaces(PlanarGraph& g, double prob);
void stellateFaces(PlanarGraph& g, PlanarGraph& gadget);
void stellateFaces(PlanarGraph& g, StellationType stype, const vector<int>& exceptVertices);
void stellateFacesOrtho(PlanarGraph& g);
void stellateFacesNH(PlanarGraph& g);
void stellateFacesNH2(PlanarGraph& g);
void stellateFacesQueue3(PlanarGraph& g);
void addTwinVertices(PlanarGraph& g);

void subdivideEdges(PlanarGraph& g, vector<Edge*>& edges);
void subdivideEdges(PlanarGraph& g, int k);
void subdivideEdges(PlanarGraph& g);

void glueCopies(PlanarGraph& graph, int kNumCopies, int A, int B, bool cyclic);
void addCopies(PlanarGraph& graph, PlanarGraph& gadget, int kNumCopies, int s, int t);

void augmentQuarticBipartite(PlanarGraph& g);

