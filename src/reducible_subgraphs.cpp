#include "graph_algorithms.h"

#include <queue>
#include <unordered_set>

bool hasReducibleTriangle(const AdjListTy& adjList);

bool hasReducibleCutPair(const AdjListTy& adjList);

bool hasReducibleTriangle(const AdjListTy& adjList) {
  auto connected = [&](const int u, const int v) -> bool {
    return adjList[u][0] == v || adjList[u][1] == v || adjList[u][2] == v;
  };
  auto remaining = [&](const int u, const int v1, const int v2) -> int {
    for (int x : adjList[u]) {
      if (x != v1 && x != v2)
        return x;
    }
    ERROR("error");
    return -1;
  };

  for (size_t i = 0; i < adjList.size(); i++) {
    for (size_t i1 = 0; i1 < adjList[i].size(); i1++) {
      for (size_t i2 = i1 + 1; i2 < adjList[i].size(); i2++) {
        const int j = adjList[i][i1];
        const int k = adjList[i][i2];
        if (!connected(j, k))
          continue;

        const int x = remaining(i, j, k);
        const int remJ = remaining(j, i, k);
        const int remK = remaining(k, i, j);

        // found a triangle connected to three distinct outer-vertices
        if (x != remJ && x != remK && remJ != remK)
          return true;

        // found two attached triangles
        if (x == remJ) {
          const int remX = remaining(x, i, j);
          if (remX != remK)
            return true;
        }
      }
    }
  }

  return false;
}

int countComponentsWithoutPair(const AdjListTy& adjList, const int u, const int v) {
  const int n = (int)adjList.size();
  std::vector<bool> removed(n, false);
  std::vector<bool> used(n, false);
  removed[u] = true;
  removed[v] = true;

  int numComponents = 0;
  std::queue<int> q;
  for (int start = 0; start < n; start++) {
    if (removed[start] || used[start]) {
      continue;
    }

    numComponents++;
    used[start] = true;
    q.push(start);
    while (!q.empty()) {
      const int now = q.front();
      q.pop();
      for (int next : adjList[now]) {
        if (removed[next] || used[next]) {
          continue;
        }
        used[next] = true;
        q.push(next);
      }
    }
  }

  return numComponents;
}

bool isComponentWithTriplePlanar(const AdjListTy& adjList, const std::vector<int>& component,
                                 const int u, const int v, const int w) {
  const int n = (int)adjList.size();
  std::vector<bool> inSubgraph(n, false);
  for (int node : component) {
    inSubgraph[node] = true;
  }
  inSubgraph[u] = true;
  inSubgraph[v] = true;
  inSubgraph[w] = true;

  std::vector<int> remap(n, -1);
  int subgraphN = 0;
  for (int node = 0; node < n; node++) {
    if (!inSubgraph[node]) {
      continue;
    }
    remap[node] = subgraphN;
    subgraphN++;
  }

  std::vector<EdgeTy> subgraphEdges;
  std::unordered_set<uint64_t> usedEdges;
  auto addEdge = [&](const int aRaw, const int bRaw) -> void {
    CHECK(aRaw != bRaw);
    const int a = std::min(aRaw, bRaw);
    const int b = std::max(aRaw, bRaw);
    const uint64_t key = (uint64_t(uint32_t(a)) << 32) | uint64_t(uint32_t(b));
    if (!usedEdges.insert(key).second) {
      return;
    }
    subgraphEdges.push_back({a, b});
  };

  for (int node = 0; node < n; node++) {
    if (remap[node] == -1) {
      continue;
    }
    for (int next : adjList[node]) {
      if (remap[next] == -1 || node >= next) {
        continue;
      }
      addEdge(remap[node], remap[next]);
    }
  }

  const int uu = remap[u];
  const int vv = remap[v];
  const int ww = remap[w];
  addEdge(uu, vv);
  addEdge(uu, ww);
  addEdge(vv, ww);

  return isPlanar(subgraphN, subgraphEdges, 0);
}

int countComponentsWithoutTriple(const AdjListTy& adjList, const int u, const int v, const int w) {
  const int n = (int)adjList.size();
  std::vector<bool> removed(n, false);
  std::vector<bool> used(n, false);
  removed[u] = true;
  removed[v] = true;
  removed[w] = true;

  int numComponents = 0;
  std::queue<int> q;
  for (int start = 0; start < n; start++) {
    if (removed[start] || used[start]) {
      continue;
    }

    std::vector<int> component;
    used[start] = true;
    q.push(start);
    while (!q.empty()) {
      const int now = q.front();
      q.pop();
      component.push_back(now);
      for (int next : adjList[now]) {
        if (removed[next] || used[next]) {
          continue;
        }
        used[next] = true;
        q.push(next);
      }
    }

    if (component.size() < 2) {
      continue;
    }

    if (isComponentWithTriplePlanar(adjList, component, u, v, w)) {
      numComponents++;
    }
  }

  return numComponents;
}

bool hasReducibleCutPair(const AdjListTy& adjList) {
  const int n = (int)adjList.size();
  for (int u = 0; u < n; u++) {
    for (int v = u + 1; v < n; v++) {
      if (countComponentsWithoutPair(adjList, u, v) > 1) {
        return true;
      }
    }
  }
  return false;
}

bool hasReducibleCutTriples(const AdjListTy& adjList) {
  const int n = (int)adjList.size();
  for (int u = 0; u < n; u++) {
    for (int v = u + 1; v < n; v++) {
      for (int w = v + 1; w < n; w++) {
        if (countComponentsWithoutTriple(adjList, u, v, w) > 0) {
          LOG("found separating triple: (%d, %d, %d)", u, v, w);
          //CHECK(false);
          return true;
        }
      }
    }
  }
  return false;
}

bool hasReducibleSubgraph(const AdjListTy& adjList) {
  // this is working only for cubic graph
  for (size_t i = 0; i < adjList.size(); i++) {
    CHECK(adjList[i].size() == 3, "not supported for non-cubic graphs");
  }

  if (hasReducibleTriangle(adjList))
    return true;

  // if (hasReducibleCutPair(adjList))
  //   return true;

  // if (hasReducibleCutTriples(adjList))
  //   return true;

  return false;
}
