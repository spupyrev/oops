#include "graph_generator.h"
#include "logging.h"
#include "common.h"

#include <algorithm>
#include <vector>
#include <queue>
#include <set>

using namespace std;

void genByClass(const std::string& graphClass, int& n, std::vector<EdgeTy>& edges) {
  if (graphClass == "complete" || graphClass == "clique") {
    genComplete(edges, n);
  } else if (graphClass == "complete-bipartite") {
    const int n1 = n / 2;
    const int n2 = n - n1;
    genCompleteBipartite(edges, n1, n2);
  } else if (graphClass == "circulant") {
    genCirculant(edges, n, {1, 2, 3});
  } else {
    ERROR("unknown graph class: '" + graphClass + "'");
  }
}

/// Generate K_n
void genComplete(std::vector<EdgeTy>& edges, int n) {
  CHECK(n >= 2);

  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      edges.push_back({i, j});
    }
  }
}

/// Generate K_n,m
void genCompleteBipartite(std::vector<EdgeTy>& edges, int n, int m) {
  CHECK(n >= 1 && m >= 1);

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      edges.push_back({i, n + j});
    }
  }
}

// regular bipartite non-planar graph
void genRegularBipartite(std::vector<EdgeTy>& edges, int n, int d) {
  CHECK(n >= 2);
  CHECK(d >= 2);

  //int n1 = Rand::next(MAX_DEGREE, n - MAX_DEGREE - 1);
  int n1 = n / 2;
  int n2 = n - n1;
  CHECK(n1 >= d && n2 >= d);
  //cerr << n1 << " " << n2 << "\n";

  auto degree = vector<int>(n, 0);
  auto adj = vector<vector<int>>(n, vector<int>(n, 0));

  bool progress = true;
  //int iter = 0;
  while (progress) {
    progress = false;
    //cerr << "iteration " << iter++ << "\n";

    auto p1 = Rand::permutation(n1);
    auto p2 = Rand::permutation(n2);
    for (size_t i1 = 0; i1 < p1.size(); i1++) {
      for (size_t i2 = 0; i2 < p2.size(); i2++) {
        int r1 = p1[i1];
        int r2 = n1 + p2[i2];
        if (r1 == r2 || adj[r1][r2]) continue;
        if (degree[r1] >= d || degree[r2] >= d) continue;
        adj[r1][r2] = adj[r2][r1] = 1;
        degree[r1]++;
        degree[r2]++;
        edges.push_back(std::make_pair(r1, r2));
        progress = true;
        break;
      }
      if (progress) break;
    }
  }

  CHECK((int)edges.size() >= d);
  for (int i = 0; i < n; i++) {
    CHECK(degree[i] <= d);
  }
}

// regular non-planar graph
void genRegular(std::vector<EdgeTy>& edges, int n, int d) {
  CHECK(n >= 2);
  CHECK(d >= 2);

  auto degree = vector<int>(n, 0);
  auto adj = vector<vector<int>>(n, vector<int>(n, 0));

  bool progress = true;
  while (progress) {
    progress = false;

    auto p1 = Rand::permutation(n);
    auto p2 = Rand::permutation(n);
    for (size_t i1 = 0; i1 < p1.size(); i1++) {
      for (size_t i2 = 0; i2 < p2.size(); i2++) {
        int r1 = p1[i1];
        int r2 = p2[i2];
        if (r1 > r2) swap(r1, r2);
        if (r1 == r2 || adj[r1][r2]) continue;
        if (degree[r1] >= d || degree[r2] >= d) continue;
        adj[r1][r2] = adj[r2][r1] = 1;
        degree[r1]++;
        degree[r2]++;
        edges.push_back(std::make_pair(r1, r2));
        progress = true;
        break;
      }
      if (progress) break;
    }
  }

  CHECK((int)edges.size() >= d);
  for (int i = 0; i < n; i++) {
    CHECK(degree[i] <= d);
  }
}

/// circulant
void genCirculant(std::vector<EdgeTy>& edges, int n, const vector<int>& S) {
  CHECK(n >= 3);
  CHECK(S.size() >= 1);

  for (int j : S) {
    CHECK(1 <= j && j <= n / 2);
  }

  for (int j : S) {
    for (int i = 0; i < n; i++) {
      int k = (i + j) % n;
      if (i != k) {
        edges.push_back({min(i, k), max(i, k)});
      }
    }
  }

  sort_unique(edges);
}
