#include "common.h"
#include "logging.h"

#include <algorithm>
#include <cassert>
#include <vector>
#include <deque>

using namespace std;

namespace planarity_lr {

struct EdgeTy {
  explicit EdgeTy() : first(-1), second(-1) {}
  explicit EdgeTy(int first, int second) : first(first), second(second) {}

  bool operator == (const EdgeTy& other) const {
    return first == other.first && second == other.second;
  }

  bool is_none() const {
    return first == -1 && second == -1;
  }

  int first = -1;
  int second = -1;
};

struct Graph {
public:
  explicit Graph() {}

  void init(bool needAM_, int n_, const std::vector<std::pair<int, int>>& edges_) {
    needAM = needAM_;
    n = n_;
    edges.reserve(edges_.size());
    for (const auto& [u, v] : edges_) {
      edges.push_back(EdgeTy(u, v));
    }
    if (needAM) {
      CHECK(edges_.empty());
      adj_init();
    }
  }

  void adj_init() {
    // adj.assign(n, std::vector<bool>(n, false));
    if ((int)adj.size() < n) {
      adj.assign(n, std::vector<bool>(n, false));
    }
    for (int i = 0; i < n; i++) {
      adj[i].assign(n, false);
    }
  }

  bool has_unordered_edge(int u, int v) const {
    CHECK(needAM);
    return adj[u][v] || adj[v][u];
  }

  void add_edge(int u, int v) {
    edges.push_back(EdgeTy(u, v));
    if (needAM)
      adj[u][v] = true;
  }

  void init_adjacency_list(std::vector<std::vector<int>>& adj_list, bool symmetric) const {
    std::vector<int> degree(n, 0);
    for (const auto& [u, v] : edges) {
      degree[u]++;
      if (symmetric)
        degree[v]++;
    }
    adj_list.assign(n, std::vector<int>());
    for (int i = 0; i < n; i++) {
      adj_list[i].reserve(degree[i]);
    }
    for (const auto& [u, v] : edges) {
      adj_list[u].push_back(v);
      if (symmetric)
        adj_list[v].push_back(u);
    }
  }

private:
  bool needAM = false;
  int n = 0;
  std::vector<EdgeTy> edges;
  static std::vector<std::vector<bool>> adj;
};

std::vector<std::vector<bool>> Graph::adj;

struct LRPlanarity;

/// Represents a set of return edges.
/// 
/// All return edges in an interval induce a same constraint on the contained
/// edges, which means that all edges must either have a left orientation or
/// all edges must have a right orientation.
struct Interval {
  explicit Interval() {}
  explicit Interval(const EdgeTy& low, const EdgeTy& high) : low(low), high(high) {}

  bool empty() const {
    // Check if the interval is empty
    return low.is_none() && high.is_none();
  }

  Interval copy() const {
    // Returns a copy of this interval
    return Interval(low, high);
  }

  bool conflicting(const EdgeTy& b, LRPlanarity *planarity_state) const;

public:
  EdgeTy low;
  EdgeTy high;  
};

/// Represents a different constraint between two intervals.
/// 
/// The edges in the left interval must have a different orientation than
/// the one in the right interval.
struct ConflictPair {
  explicit ConflictPair() {}
  explicit ConflictPair(const Interval& left, const Interval& right) : left(left), right(right) {}

  void swap() {
    // Swap left && right intervals
    std::swap(left, right);
  }

  int lowest(LRPlanarity *planarity_state) const;

  bool operator == (const ConflictPair& other) const {
    return left.low == other.left.low && 
           left.high == other.left.high && 
           right.low == other.right.low && 
           right.high == other.right.high;
  }

public:
  Interval left;
  Interval right;  
};        

struct LRPlanarity {
  explicit LRPlanarity(const int n, const std::vector<std::pair<int, int>>& edges) 
      : n(n), m(edges.size()) {
    CHECK(n >= 5 && m >= 9);
    CHECK(m <= 3 * n - 6);
    // init data
    G.init(false, n, edges);
    // init2
    DG.init(true, n, std::vector<std::pair<int, int>>());
    // vertex init
    height.resize(n, -1);
    parent_edge.resize(n);
    next_index.resize(n);
    // edge init
    edgeData.assign(n, std::vector<EdgeData>(n));
  }

  bool lr_planarity() {
    // make adjacency lists for dfs
    G.init_adjacency_list(adjs, true);

    // orientation of the graph by depth first search traversal
    for (int v = 0; v < n; v++) {
      if (height[v] == -1) {
        height[v] = 0;
        roots.push_back(v);
        dfs_orientation(v);
      }
    }

    DG.init_adjacency_list(ordered_adjs, false);

    // sort the adjacency lists by nesting depth
    // note: this sorting leads to non linear time
    for (int v = 0; v < n; v++) {
      if (ordered_adjs[v].size() <= 1)
        continue;
      std::sort(ordered_adjs[v].begin(), ordered_adjs[v].end(), [&](int l, int r) {
        const auto& lData = edgeData[v][l];
        const auto& rData = edgeData[v][r];
        return lData.nesting_depth < rData.nesting_depth;
      });
    }

    // planarity testing
    for (int v : roots) {
      if (!dfs_testing(v)) {
        return false;
      }
    }

    return true;
  }

  /// Orient the graph by DFS, compute lowpoints && nesting order
  void dfs_orientation(const int root) { 
    // the recursion stack
    dfs_stack = {root};
    // index of next edge to handle in adjacency list of each node
    next_index.assign(n, 0);
    // boolean to indicate whether to skip the initial work for an edge
    clear_skip_init();

    while (!dfs_stack.empty()) {
      const int v = dfs_stack.back();
      dfs_stack.pop_back();
      const EdgeTy& e = parent_edge[v];
      EdgeData& eData = edgeData[e.first][e.second];

      for (size_t i = next_index[v]; i < adjs[v].size(); i++) {
        const int w = adjs[v][i];
        const EdgeTy vw(v, w);
        EdgeData& vwData = edgeData[v][w];

        if (!vwData.skip_init) {
          if (DG.has_unordered_edge(v, w)) {
            next_index[v] += 1;
            continue;  // the edge was already oriented
          }

          DG.add_edge(v, w);  // orient the edge

          vwData.lowpt = height[v];
          vwData.lowpt2 = height[v];
          if (height[w] == -1) {  // (v, w) is a tree edge
            parent_edge[w] = vw;
            height[w] = height[v] + 1;

            dfs_stack.push_back(v);  // revisit v after finishing w
            dfs_stack.push_back(w);  // visit w next
            vwData.skip_init = true; // don't redo this block
            break;  // handle next node in dfs_stack (i.e. w)
          } else {  // (v, w) is a back edge
            vwData.lowpt = height[w];
          }
        }

        // determine nesting graph
        vwData.nesting_depth = 2 * vwData.lowpt;
        if (vwData.lowpt2 < height[v])  // chordal
          vwData.nesting_depth += 1;

        // update lowpoints of parent edge e
        if (!e.is_none()) {
          if (vwData.lowpt < eData.lowpt) {
            eData.lowpt2 = std::min(eData.lowpt, vwData.lowpt2);
            eData.lowpt = vwData.lowpt;
          } else if (vwData.lowpt > eData.lowpt) {
            eData.lowpt2 = std::min(eData.lowpt2, vwData.lowpt);
          } else {
            eData.lowpt2 = std::min(eData.lowpt2, vwData.lowpt2);
          }
        }       

        next_index[v] += 1;
      }
    }
  }

  /// Test for LR partition
  bool dfs_testing(const int root) {
    // the recursion stack
    dfs_stack = {root};
    // index of next edge to handle in adjacency list of each node
    next_index.assign(n, 0);
    // boolean to indicate whether to skip the initial work for an edge
    clear_skip_init();

    while (!dfs_stack.empty()) {
      const int v = dfs_stack.back();
      dfs_stack.pop_back();
      const EdgeTy& e = parent_edge[v];
      EdgeData& eData = edgeData[e.first][e.second];
      // to indicate whether to skip the final block after the for loop
      bool skip_final = false;

      for (size_t i = next_index[v]; i < ordered_adjs[v].size(); i++) {
        const int w = ordered_adjs[v][i];
        const EdgeTy vw = EdgeTy(v, w);
        EdgeData& vwData = edgeData[v][w];

        if (!vwData.skip_init) {          
          vwData.stack_bottom = top_of_stack(S);

          if (vw == parent_edge[w]) {  // tree edge
            dfs_stack.push_back(v);  // revisit v after finishing w
            dfs_stack.push_back(w);  // visit w next
            vwData.skip_init = true;  // don't redo this block
            skip_final = true;  // skip final work after breaking
            break;  // handle next node in dfs_stack (i.e. w)
          } else {  // back edge
            vwData.lowpt_edge = vw;
            S.push_back(ConflictPair(Interval(), Interval(vw, vw)));
          }
        }

        // integrate new return edges
        if (vwData.lowpt < height[v]) {
          if (w == ordered_adjs[v][0]) {  // e_i has return edge
            eData.lowpt_edge = vwData.lowpt_edge;
          } else {  // add constraints of e_i
            if (!add_constraints(vw, e)) {
              // graph is not planar
              return false;
            }
          }
        }

        next_index[v] += 1;
      }

      if (!skip_final) {
        // remove back edges returning to parent
        if (!e.is_none()) {  // v isn't root
          remove_back_edges(e);
        }
      }
    }

    return true;
  }

  bool add_constraints(const EdgeTy& vw, const EdgeTy& e) {
    const auto& eData = edgeData[e.first][e.second];
    const auto& vwData = edgeData[vw.first][vw.second];
    auto P = ConflictPair(Interval(), Interval());
    // merge return edges of e_i into P.right
    while (true) {
      auto Q = S.back();
      S.pop_back();

      if (!Q.left.empty()) {
        Q.swap();
      }
      if (!Q.left.empty()) { // not planar
        return false;
      }

      auto& qrlData = edgeData[Q.right.low.first][Q.right.low.second];
      if (qrlData.lowpt > eData.lowpt) {
        // merge intervals
        if (P.right.empty()) {  // topmost interval
          P.right = Q.right.copy();
        } else {
          auto& prlData = edgeData[P.right.low.first][P.right.low.second];
          prlData.ref = Q.right.high;
        }
        P.right.low = Q.right.low;
      } else {  // align
        qrlData.ref = eData.lowpt_edge;
      }      
      if (top_of_stack(S) == vwData.stack_bottom) {
        break;
      }
    }

    // merge conflicting return edges of e_1,...,e_i-1 into P.L
    while (true) {
      const auto topS = top_of_stack(S);
      if (!topS.left.conflicting(vw, this) && !topS.right.conflicting(vw, this))
        break;

      auto Q = S.back();
      S.pop_back();

      if (Q.right.conflicting(vw, this)) {
        Q.swap();
      }
      if (Q.right.conflicting(vw, this)) { // not planar
        return false;
      }

      // merge interval below lowpt(e_i) into P.R
      auto& prlData = edgeData[P.right.low.first][P.right.low.second];
      prlData.ref = Q.right.high;
      if (!Q.right.low.is_none()) {
        P.right.low = Q.right.low;
      }

      if (P.left.empty()) {  // topmost interval
        P.left = Q.left.copy();
      } else {
        auto& pllData = edgeData[P.left.low.first][P.left.low.second];
        pllData.ref = Q.left.high;
      }
      P.left.low = Q.left.low;
    }

    if (!P.left.empty() || !P.right.empty()) {
      S.push_back(P);
    }
    return true;
  }

  void remove_back_edges(const EdgeTy& e) {
    auto& eData = edgeData[e.first][e.second];
    const int u = e.first;
    // trim back edges ending at parent u
    // drop entire conflict pairs
    while (!S.empty() && top_of_stack(S).lowest(this) == height[u]) {
      S.pop_back();
    }

    if (!S.empty()) {  // one more conflict pair to consider
      auto P = S.back();
      S.pop_back();

      // trim left interval
      while (!P.left.high.is_none() && P.left.high.second == u) {
        const auto& plhData = edgeData[P.left.high.first][P.left.high.second];
        P.left.high = plhData.ref;
      }
      if (P.left.high.is_none() && !P.left.low.is_none()) {
        // just emptied
        auto& pllData = edgeData[P.left.low.first][P.left.low.second];
        pllData.ref = P.right.low;
        P.left.low = EdgeTy();
      }
      // trim right interval
      while (!P.right.high.is_none() && P.right.high.second == u) {
        const auto& prhData = edgeData[P.right.high.first][P.right.high.second];
        P.right.high = prhData.ref;
      }
      if (P.right.high.is_none() && !P.right.low.is_none()) {
        // just emptied
        auto& prlData = edgeData[P.right.low.first][P.right.low.second];
        prlData.ref = P.left.low;
        P.right.low = EdgeTy();
      }
      S.push_back(P);
    }

    // side of e is side of a highest return edge
    if (eData.lowpt < height[u]) {  // e has return edge
      EdgeTy hl = top_of_stack(S).left.high;
      EdgeTy hr = top_of_stack(S).right.high;

      if (!hl.is_none() && (hr.is_none() || edgeData[hl.first][hl.second].lowpt > edgeData[hr.first][hr.second].lowpt)) {
        eData.ref = hl;
      } else {
        eData.ref = hr;
      }
    }
  }

  /// Returns the element on top of the stack
  ConflictPair top_of_stack(const std::deque<ConflictPair>& l) const {  
    return !l.empty() ? l.back() : ConflictPair();
  }

  void clear_skip_init() {
    for (int i = 0; i < n; i++) {
      for (int j : adjs[i]) {
        edgeData[i][j].skip_init = false;
      }
    }
  }

public:
  const int n;
  const int m;
  // input graph
  Graph G;
  // oriented DFS graph
  Graph DG;

  // adjacency lists
  std::vector<vector<int>> adjs;
  std::vector<vector<int>> ordered_adjs;

  // Vertex data
  std::vector<int> roots;
  std::vector<int> height;
  std::vector<EdgeTy> parent_edge;
  std::vector<size_t> next_index;

  // Edge data
  struct EdgeData {
    // height of lowest return point of an edge
    int lowpt = 0;
    // height of second lowest return point
    int lowpt2 = 0;
    // for nesting order
    int nesting_depth = 0;
    // 
    EdgeTy ref;
    // conflict pairs
    ConflictPair stack_bottom;
    EdgeTy lowpt_edge;
    // 
    bool skip_init;
  };

  // making it static to avoid re-allocation of the data
  static std::vector<std::vector<EdgeData>> edgeData;

  // dfs stack
  std::deque<int> dfs_stack;
  // stack of conflict pairs
  std::deque<ConflictPair> S;
};

std::vector<std::vector<LRPlanarity::EdgeData>> LRPlanarity::edgeData;

/// Returns True if interval I conflicts with edge b
bool Interval::conflicting(const EdgeTy& b, LRPlanarity *planarity_state) const { 
  if (empty()) 
    return false;

  const auto& highData = planarity_state->edgeData[high.first][high.second];
  const auto& bData = planarity_state->edgeData[b.first][b.second];
  return highData.lowpt > bData.lowpt;
}

/// Returns the lowest lowpoint of a conflict pair
int ConflictPair::lowest(LRPlanarity *planarity_state) const {
  if (left.empty()) {
    const auto& data = planarity_state->edgeData[right.low.first][right.low.second];
    return data.lowpt;
  }
  if (right.empty()) {
    const auto& data = planarity_state->edgeData[left.low.first][left.low.second];
    return data.lowpt;
  }
  const auto& dataL = planarity_state->edgeData[left.low.first][left.low.second];
  const auto& dataR = planarity_state->edgeData[right.low.first][right.low.second];
  return std::min(dataL.lowpt, dataR.lowpt);
}

} // namespace planarity_lr

bool isPlanar(const int n, const std::vector<std::pair<int, int>>& edges, int verbose) {  
  // Check if the input is correct
  for (size_t i = 0; i < edges.size(); i++) {
    const auto& ei = edges[i];
    CHECK(ei.first < ei.second);
    CHECK(0 <= ei.first && ei.first < n);
    CHECK(0 <= ei.second && ei.second < n);
    for (size_t j = i + 1; j < edges.size(); j++) {
      const auto& ej = edges[j];
      CHECK(ei.first != ej.first || ei.second != ej.second);
    }
  }

  // Stop early, if possible
  const int m = (int)edges.size();
  if (n < 5 || m < 9)
    return true;
  if (m > 3 * n - 6)
    return false;

  // Planarity testing
  auto alg = planarity_lr::LRPlanarity(n, edges);
  const bool res = alg.lr_planarity();

  if (res)
    LOG_IF(verbose >= 2, "passed planarity test for |V| = %d and |E| = %d", n, edges.size());
  else
    LOG_IF(verbose >= 2, "non-planar graph with |V| = %d and |E| = %d", n, edges.size());

  return res;
}
