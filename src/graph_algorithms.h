#pragma once

#include "adjacency.h"

#include <vector>
#include <set>

void printGraphStats(const int n, const std::vector<EdgeTy>& edges, const std::vector<bool>& direction);

bool isPlanar(const int n, const std::vector<std::pair<int, int>>& edges, int verbose);

/// Returns depth for every vertex
std::vector<int> bfs(const int n, const std::vector<EdgeTy>& edges, const std::vector<int>& seeds);

/// Returns the dfs order
std::vector<int> dfs(const std::vector<std::vector<bool>>& adj, int seed);

/// Check if the graph is connected
bool isConnected(int n, const std::vector<EdgeTy>& edges);

/// Check if the graph remains connected after removing some vertices
bool isConnected(const std::vector<int>& vList, const std::vector<std::vector<int>>& adj, std::set<int>& removed);

/// Return sets of nodes, one set for each biconnected component of the graph
std::vector<std::vector<EdgeTy>> biconnectedComponents(const int n, const std::vector<EdgeTy>& edges);

/// Check if the graph is 2-(vertex)-connected
bool is2Connected(int n, const std::vector<EdgeTy>& edges);

/// Check if the graph is bipartite
bool isBipartite(int n, const std::vector<EdgeTy>& edges);

/// Check if the graph is a (possibly disconnected) tree (forest)
bool isTree(const int n, const std::vector<EdgeTy>& edges);

/// Check if the graph is a 3-tree
bool is3Tree(const int n, const std::vector<EdgeTy>& edges);

/// Check if the graph is a 3-tree
bool is3Tree(const std::vector<std::vector<int>>& adj);

/// Check if the provided graph is acyclic
bool isAcyclic(const int n, const std::vector<EdgeTy>& edges, int verbose);

int minDegree(const int n, const std::vector<EdgeTy>& edges);

int maxDegree(const int n, const std::vector<EdgeTy>& edges);

int computeGirth(const int n, const std::vector<EdgeTy>& edges);

bool hasReducibleTriangle(const AdjListTy& adjList);
