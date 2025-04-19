#include "graph_algorithms.h"
#include "adjacency.h"
#include "common.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <set>
#include <queue>
#include <map>

using namespace std;

/// Print graph stats: n, m, degree, girth
void printGraphStats(ostream& out, const int n, const std::vector<EdgeTy>& edges, const std::vector<bool>& direction) {
  int numSinks = 0, numSources = 0;
  if (!direction.empty()) {
    CHECK(edges.size() == direction.size());
    auto inDegree = vector<int>(n, 0);
    auto outDegree = vector<int>(n, 0);
    for (size_t i = 0; i < edges.size(); i++) {
      if (direction[i]) {
        outDegree[edges[i].first]++;
        inDegree[edges[i].second]++;
      } else {
        outDegree[edges[i].second]++;
        inDegree[edges[i].first]++;
      }
    }
    for (int i = 0; i < n; i++) {
      if (inDegree[i] == 0) {
        numSources++;
      }
      if (outDegree[i] == 0) {
        numSinks++;
      }
    }
  }

  const int minDeg = minDegree(n, edges);
  const int maxDeg = maxDegree(n, edges);
  const bool bipartite = isBipartite(n, edges);
  const bool tree = isTree(n, edges);
  const int girth = computeGirth(n, edges);

  out << "|V| = " << n;
  out << "   |E| = " << edges.size();
  out << "   " << minDeg << " <= deg <= " << maxDeg;
  out << "   bipartite = " << bipartite;
  out << "   tree = " << tree;
  out << "   girth = " << girth;
  if (!direction.empty()) {
    out << "   |sources| = " << numSources;
    out << "   |sinks| = " << numSinks;
  }
}

void printGraphStats(const int n, const std::vector<EdgeTy>& edges, const std::vector<bool>& direction) {
  stringstream ss;
  printGraphStats(ss, n, edges, direction);
  LOG(ss.str());
}

std::vector<int> bfs(const int n, const std::vector<EdgeTy>& edges, const std::vector<int>& seeds) {
  CHECK(!seeds.empty());
  vector<int> depth(n, -1);
  queue<int> q;
  for (int seed : seeds) {
    depth[seed] = 0;
    q.push(seed);
  }

  Adjacency adj(n);
  adj.from_edges(edges);

  while (!q.empty()) {
    int now = q.front();
    q.pop();

    adj.forEach(now, [&](int next) {
      CHECK(now != next);
      if (depth[next] == -1) {
        depth[next] = depth[now] + 1;
        q.push(next);
      }
    });
  }

  return depth;
}

void dfs(int now, const vector<vector<bool>>& adj, vector<int>& order, vector<bool>& used) {
  CHECK(!used[now]);
  order.push_back(now);
  used[now] = true;

  for (int i = 0; i < (int)adj[now].size(); i++) {
    if (!adj[now][i])
      continue;

    if (used[i])
      continue;

    dfs(i, adj, order, used);
  }
}

std::vector<int> dfs(const std::vector<std::vector<bool>>& adj, int seed) {
  vector<int> order;
  vector<bool> used(adj.size(), false);
  dfs(seed, adj, order, used);
  return order;
}

bool isBipartite(int n, const std::vector<EdgeTy>& edges, vector<int>& part1, vector<int>& part2) {
  vector<int> depth(n, -1);
  Adjacency adj(n);
  adj.from_edges(edges);

  for (int i = 0; i < n; i++) {
    if (depth[i] != -1) continue;
    int now = i;
    depth[now] = 0;
    queue<int> q;
    q.push(now);

    while (!q.empty()) {
      now = q.front();
      q.pop();

      adj.forEach(now, [&](int next) {
        CHECK(now != next);
        if (depth[next] == -1) {
          depth[next] = depth[now] + 1;
          q.push(next);
        }
      });
    }
  }

  for (auto& edge : edges) {
    int color1 = depth[edge.first] % 2;
    int color2 = depth[edge.second] % 2;

    if (color1 == color2)
      return false;
  }

  for (int i = 0; i < n; i++) {
    if (depth[i] % 2 == 0)
      part1.push_back(i);
    else
      part2.push_back(i);
  }

  return true;
}

bool isBipartite(int n, const std::vector<EdgeTy>& edges) {
  vector<int> part1;
  vector<int> part2;
  return isBipartite(n, edges, part1, part2);
}

bool isConnected(const vector<int>& vList, const vector<vector<int>>& adj, set<int>& removed) {
  queue<int> q;
  set<int> processed;

  for (size_t i = 0; i < vList.size(); i++) {
    if (!removed.count(vList[i])) {
      q.push(vList[i]);
      processed.insert(vList[i]);
      break;
    }
  }

  while (!q.empty()) {
    int now = q.front();
    q.pop();

    for (size_t i = 0; i < adj[now].size(); i++) {
      int next = adj[now][i];

      if (removed.count(next)) {
        continue;
      }

      if (processed.find(next) != processed.end()) {
        continue;
      }

      q.push(next);
      processed.insert(next);
    }
  }

  return (processed.size() + removed.size() == vList.size());
}

bool isConnected(int n, const std::vector<EdgeTy>& edges) {
  auto color = bfs(n, edges, {0});
  for (int i = 0; i < n; i++) {
    if (color[i] == -1)
      return false;
  }
  return true;
}

void bicompDFS(
    const Adjacency& adj,
    const int now, 
    int& cur_time,
    std::vector<int>& discovery, 
    std::vector<int>& low, 
    std::vector<int>& parent,
    std::vector<EdgeTy>& edge_stack,
    std::vector<std::vector<EdgeTy>>& bicomponents) {
  // Initialize discovery time and low value
  discovery[now] = low[now] = ++cur_time;
  int numChildren = 0;

  // Go through all vertices adjacent to now
  adj.forEach(now, [&](const int next) {
    CHECK(now != next);

    if (discovery[next] == -1) {
      // If v is not visited yet, then recurse for it
      numChildren++;
      parent[next] = now;
      // store the edge in stack
      edge_stack.push_back({now, next});
      bicompDFS(adj, next, cur_time, discovery, low, parent, edge_stack, bicomponents);

      // Check if the subtree rooted with 'v' has a connection to one of the ancestors of 'now'
      low[now] = std::min(low[now], low[next]);

      // If u is an articulation point, pop all edges from stack till (now, v)
      if ((discovery[now] == 1 && numChildren > 1) || (discovery[now] > 1 && low[next] >= discovery[now])) {
        while (edge_stack.back().first != now || edge_stack.back().second != next) {
          const int u = std::min(edge_stack.back().first, edge_stack.back().second);
          const int v = std::max(edge_stack.back().first, edge_stack.back().second);
          //std::cout << u << "--" << v << " ";
          bicomponents.back().push_back({u, v});
          edge_stack.pop_back();
        }
        const int u = std::min(edge_stack.back().first, edge_stack.back().second);
        const int v = std::max(edge_stack.back().first, edge_stack.back().second);
        //std::cout << u << "--" << v << " ";
        bicomponents.back().push_back({u, v});
        edge_stack.pop_back();
        //std::cout << endl;
        bicomponents.push_back({});
      }
    } else if (next != parent[now]) {
      // Update low value of 'u' only of 'v' is still in stack
      // (i.e. it's a back edge, not cross edge).
      low[now] = std::min(low[now], discovery[next]);
      if (discovery[next] < discovery[now]) {
        edge_stack.push_back({now, next});
      }
    }
  });
}

/// Depth-first search algorithm to generate articulation points and biconnected components
std::vector<std::vector<EdgeTy>> biconnectedComponents(const int n, const std::vector<EdgeTy>& edges) {
  // Initialize disc and low, and parent arrays
  std::vector<int> discovery(n, -1);
  std::vector<int> low(n, -1);
  std::vector<int> parent(n, -1);
  std::vector<EdgeTy> edge_stack;
  int cur_time = 0;

  Adjacency adj(n);
  adj.from_edges(edges);

  std::vector<std::vector<EdgeTy>> bicomponents;
  bicomponents.push_back({});

  for (int i = 0; i < n; i++) {
    if (discovery[i] == -1) {
      bicompDFS(adj, i, cur_time, discovery, low, parent, edge_stack, bicomponents);
    }

    // If the edge stack is not empty, pop all edges from the stack
    while (!edge_stack.empty()) {
      const int u = std::min(edge_stack.back().first, edge_stack.back().second);
      const int v = std::max(edge_stack.back().first, edge_stack.back().second);
      bicomponents.back().push_back({u, v});
      edge_stack.pop_back();
    }
    if (!bicomponents.back().empty()) {
      bicomponents.push_back({});
    }
  }

  CHECK(bicomponents.back().empty());
  bicomponents.pop_back();

  // Verify the result
  size_t numEdges = 0;
  for (auto& comp : bicomponents) {
    CHECK(comp.size() >= 1);
    numEdges += comp.size();
  }
  CHECK(numEdges == edges.size());

  return bicomponents;
}

bool is2Connected(int n, const std::vector<EdgeTy>& edges) {
  CHECK(n > 0);
  auto biComponents = biconnectedComponents(n, edges);
  return biComponents.size() == 1;
}

bool isTree(const int n, const std::vector<EdgeTy>& edges) {
  Adjacency adj(n);
  adj.from_edges(edges);
  std::vector<int> degree(n, 0);
  for (int i = 0; i < n; i++) {
    degree[i] = adj.degree(i);
  }

  std::vector<bool> alive(n, true);
  int numAlive = n;
  while (numAlive > 0) {
    bool progress = false;
    for (int i = 0; i < n; i++) {
      if (!alive[i]) continue;
      if (degree[i] > 1) continue;
      
      progress = true;
      alive[i] = false;
      numAlive--;
      adj.forEach(i, [&](int j) {
        if (alive[j]) {
          CHECK(degree[j] >= 1);
          degree[j]--;
        }
      });
    }

    if (!progress) 
      return false;
  }

  return true;
}

bool is3Tree(const int n, const std::vector<EdgeTy>& edges) {
  AdjListTy adj(n, std::vector<int>());
  for (auto& [src, dst] : edges) {
    CHECK(0 <= src && src < n);
    CHECK(0 <= dst && dst < n);
    CHECK(src < dst);
    adj[src].push_back(dst);
    adj[dst].push_back(src);
  }

  int N = n;
  while (N > 3) {
    bool removed = false;

    for (int i = 0; i < n; i++) {
      if (adj[i].empty())
        continue;
      CHECK(adj[i].size() >= 3);
      if (adj[i].size() != 3)
        continue;

      // remove i
      for (size_t j = 0; j < adj[i].size(); j++) {
        int r = adj[i][j];
        CHECK(i != r);
        adj[r].erase(remove(adj[r].begin(), adj[r].end(), i), adj[r].end());
      }

      adj[i].clear();
      removed = true;
      N--;
      break;
    }

    if (!removed)
      return false;
  }

  return true;
}

bool is3Tree(const std::vector<vector<int>>& adj0) {
  vector<vector<int>> adj = adj0;
  int n = adj.size();

  while (n > 3) {
    bool removed = false;

    for (size_t i = 0; i < adj.size(); i++) {
      if (adj[i].size() == 0)
        continue;

      CHECK(adj[i].size() >= 3);

      if (adj[i].size() != 3)
        continue;

      // remove i
      for (size_t j = 0; j < adj[i].size(); j++) {
        int r = adj[i][j];
        CHECK((int)i != r);
        adj[r].erase(remove(adj[r].begin(), adj[r].end(), i), adj[r].end());
      }

      adj[i].clear();
      removed = true;
      n--;
      //cout << "removing " << i << "\n";
      break;
    }

    if (!removed)
      return false;
  }

  std::vector<int> clique;
  for (size_t i = 0; i < adj.size(); i++) {
    if (adj[i].size() > 0)
      clique.push_back(i);
  }
  LOG("parent clique for 3-tree: %s", to_string(clique).c_str());
  return true;
}

bool isAcyclic(const int n, const std::vector<EdgeTy>& edges, int verbose) {
  vector<vector<bool>> adj(n, vector<bool>(n, false));
  for (const auto& [s, t] : edges) {
    adj[s][t] = true;
  }

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (adj[i][j] && adj[j][i]) {
        LOG_IF(verbose >= 2, "found a cycle for (%d, %d)", i, j);
        return false;
      }
    }
  }

  for (int k = 0; k < n; k++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (i != j && i != k && j != k) {
          if (adj[i][k] && adj[k][j]) {
            if (adj[j][i]) {
              LOG_IF(verbose >= 2, "found a cycle for (%d, %d)", i, j);
              return false;
            }
            adj[i][j] = true;
          }
        }
      }
    }
  }
  return true;
}

int minDegree(const int n, const std::vector<EdgeTy>& edges) {
  std::vector<int> degree(n, 0);
  for (auto& [u, v] : edges) {
    degree[u]++;
    degree[v]++;
  }
  return *std::min_element(degree.begin(), degree.end());
}

int maxDegree(const int n, const std::vector<EdgeTy>& edges) {
  std::vector<int> degree(n, 0);
  for (auto& [u, v] : edges) {
    degree[u]++;
    degree[v]++;
  }
  return *std::max_element(degree.begin(), degree.end());
}

int computeGirth(const int n, const std::vector<EdgeTy>& edges) {
  Adjacency adj(n);
  adj.from_edges(edges);

  using Node = std::pair<int, int>; // vertex, depth

  std::vector<int> labels(n);
  int best = n - 1;

  std::queue<Node> queue;

  // Start a BFS from every vertex except the last two (not needed)
  int root = 0;
  while (root < n - 2 && best > 3) {
    std::fill(labels.begin(), labels.end(), -1);

    // Add initial node to the queue
    labels[root] = 0;
    queue.push({root, 0});

    while (!queue.empty() && best > 3) {
      Node node = queue.front();
      queue.pop();
      if ((node.second + 1) * 2 - 1 >= best)
        break;
      const int depth = node.second + 1;

      // Check all neighbours
      adj.forEach(node.first, [&](int neighbour) {
        CHECK(node.first != neighbour);
        // We haven't seen this neighbour before
        if (labels[neighbour] < 0) {
          queue.push({neighbour, depth});
          labels[neighbour] = depth;
        // Cycle with odd number of edges
        } else if (labels[neighbour] == depth - 1) {
          if (depth * 2 - 1 < best)
            best = depth * 2 - 1;
        // Cycle with even number of edges
        } else if (labels[neighbour] == depth) {
            if (depth * 2 < best)
              best = depth * 2;
        }
      });
    }

    // Clear the queue and prepare to start a BFS from a next vertex
    queue = {};
    root++;
  }

  return best > 0 ? best : 1;
}  
