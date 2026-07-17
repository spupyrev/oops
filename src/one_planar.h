#pragma once

#include "common.h"
#include "logging.h"
#include "forbidden_crossings.h"
#include "glucose/SolverSimp21.h"
#include "breakid/sat_symmetry.h"
#include "satsuma/sat_symmetry.h"
#include "satsuma/tsl.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <zlib.h>

using namespace Simp21;

using EdgeTy = std::pair<int, int>;

constexpr size_t MAX_CROSSINGS_FOR_2GT = 3840;
constexpr size_t MAX_CROSSINGS_FOR_3EQ = 1400;
constexpr size_t MAX_CROSSINGS_FOR_3GT = 1280;
constexpr int MAX_CYCLE_LEN_FOR_SEPCYCLES = 9;

inline EdgeTy make_edge(const int u, const int v) {
  CHECK(u != v);
  return EdgeTy(std::min(u, v), std::max(u, v));
}

/// The object is a representation of the input graph that contains n vertices and m edges.
///
/// InputGraph has n "regular" vertices and m "division" vertices. (so every edge is subdivided once)
/// There are 2*m edge segments; i-th edge is split into (2*i)-th and (2*i+1)-th segments
/// Regular vertices are in [0, n); division vertices are in [n, n + m)
struct InputGraph {
  int n = 0;
  // vertices are in [0..n); first < second
  std::vector<EdgeTy> edges;
  // adjacency list
  std::vector<std::vector<int>> adj;
  // edge directions: empty for undirected, (true => first->second, false => second->first)
  std::vector<bool> directions;
  // Maps an original (parent-graph) vertex index to this graph's local index, or
  // -1 if the vertex is absent here. Identity for a top-level graph; a proper
  // remap when this graph is a re-indexed biconnected component. Lets external
  // (parent-label) inputs such as -fix-cross1 resolve to the right local edge.
  std::vector<int> origToLocal;

  InputGraph() {}

  /// Create a graph for a given number of vertices and given edges
  InputGraph(int n, const std::vector<EdgeTy>& edges, const std::vector<bool>& directions) :
      n(n), edges(edges), directions(directions) {
    origToLocal.resize(n);
    for (int i = 0; i < n; i++) origToLocal[i] = i;
    init();
  }

  /// Create a graph for a given set of edges (possibly by re-indexing its vertices)
  InputGraph(const std::vector<EdgeTy>& edgesC) {
    int maxN = 0;
    for (const auto& [u, v] : edgesC) {
      maxN = std::max(maxN, u);
      maxN = std::max(maxN, v);
    }
    std::vector<bool> used(maxN + 1, false);
    for (const auto& [u, v] : edgesC) {
      used[u] = true;
      used[v] = true;
    }

    n = 0;
    std::vector<int> index(maxN + 1, -1);
    for (int i = 0; i < maxN + 1; i++) {
      if (!used[i]) continue;
      index[i] = n;
      n++;
    }
    for (const auto& [u, v] : edgesC) {
      CHECK(index[u] != -1 && index[v] != -1);
      CHECK(index[u] < index[v], "u = %d; v = %d", u, v);
      edges.push_back({index[u], index[v]});
    }
    origToLocal = index;

    init();
  }

  void init() {
    CHECK(edges.size() > 0, "empty input graph");
    CHECK(directions.empty() || directions.size() == edges.size());
    for (const auto& [u, v] : edges) {
      CHECK(0 <= u && u < n);
      CHECK(0 <= v && v < n);
      CHECK(u < v);
    }

    adj.resize(n);
    for (const auto& [u, v] : edges) {
      adj[u].push_back(v);
      adj[v].push_back(u);
    }
  }

  /// Check if the input graph is directed
  bool isDirected() const {
    return !directions.empty();
  }

  /// Minimum vertex degree
  int minDegree() const {
    int res = (int)adj[0].size();
    for (int i = 1; i < n; i++) {
      res = std::min(res, (int)adj[i].size());
    }
    return res;
  }

  /// Vertex degree
  int degree(int v) const {
    return (int)adj[v].size();
  }

  /// Find the index of edge (u, v) or "-1" if it doesn't exist
  int findEdgeIndex(const int u, const int v) const {
    for (int i = 0; i < (int)edges.size(); i++) {
      const auto& [s, t] = edges[i];
      if ((s == u && t == v) || (s == v && t == u))
        return i;
    }
    return -1;
  }

  /// Find the index of the division vertex for edge (u, v)
  int findDivIndex(int u, int v) const {
    const int edgeIndex = findEdgeIndex(u, v);
    CHECK(edgeIndex != -1);
    return edgeIndex + n;
  }

  /// Check if the two vertices are connected
  bool adjacentVV(const int u, const int v) const {
    CHECK(0 <= u && u < n);
    CHECK(0 <= v && v < n);
    return contains(adj[u], v);
  }

  /// Check if the edge adjacent to the vertex
  bool adjacentEV(const int edgeIndex, const int nodeIndex) const {
    CHECK(0 <= edgeIndex && edgeIndex < (int)edges.size());
    CHECK(0 <= nodeIndex && nodeIndex < n + (int)edges.size());
    return (edges[edgeIndex].first == nodeIndex || edges[edgeIndex].second == nodeIndex);
  }

  /// Check if the two edges are adjacent
  bool adjacentEE(const EdgeTy& edgeI, const EdgeTy& edgeJ) const {
    if (edgeI.first == edgeJ.first)
      return true;
    if (edgeI.first == edgeJ.second)
      return true;
    if (edgeI.second == edgeJ.first)
      return true;
    if (edgeI.second == edgeJ.second)
      return true;
    return false;
  }

  /// Return endpoints of a segment
  std::pair<int, int> seg2edge(const int seg) const {
    CHECK(0 <= seg && seg < 2 * int(edges.size()));
    const int src = seg % 2;
    const int edge_idx = seg / 2;
    const auto& ee = edges[edge_idx];
    const int e_first = src ? ee.first : edge_idx + n;
    const int e_second = src ? edge_idx + n : ee.second;
    return std::make_pair(e_first, e_second);
  }

  /// Return endpoints of a segment
  std::pair<int, int> seg2edge_v2(const int seg) const {
    CHECK(0 <= seg && seg < 2 * int(edges.size()));
    const int src = seg % 2;
    const int edge_idx = seg / 2;
    const auto& ee = edges[edge_idx];
    if (src)
      return std::make_pair(ee.first, edge_idx + n);
    else
      return std::make_pair(ee.second, edge_idx + n);
  }
};

/// SAT solver parameters
enum SolverType { STACK, MOVE, BRUTE_FORCE };
struct Params {
  Params() = default;

  Params(const Params&) = delete;
  Params& operator=(const Params&) = delete;
  Params(Params&&) = delete;
  Params& operator=(Params&&) = delete;  

  // Verbosity level
  int verbose = 0;
  // Timeout (in seconds)
  int timeout = 0;
  // Force generating a solution
  bool alwaysCreateSolution = false;

  // Whether to apply symmetry-breaking constraints
  bool applyBreakID = false;
  bool applySatsuma = false;

  // Dimacs input/output
  std::string modelFile = "";
  std::string resultFile = "";

  SolverType solverType = SolverType::STACK;
  bool useSATConstraints = false;
  bool useUNSATConstraints = false;
  bool useIC = false;
  bool useNIC = false;

  int strict = 0;
  bool forbidCrossings = false;
  std::string swapConstraints = "";
  std::string partialConstraints = "";
  std::string sepCycleConstraints = "";
  int unsatLevel = 0;

  std::string custom = "";
  bool ignoreTransitiveRels = false;
  bool useSepCycleUP = false;
  // Bias CDCL toward branching on cross1/cross2 vars first, with cross1=true
  // preferred. Strong on g.31.9-class instances; mixed on rome-style ones.
  bool crossPriority = false;

  // Fixed cross1 units: ';'-separated 'u:v+/-' (edge {u,v} crossed/uncrossed).
  // Negative units also prune crossablePairs pre-encoding.
  std::string fixCross1 = "";
  bool cubicVerification = false;

  bool hasCross1Restrictions() const {
    return !fixCross1.empty() || cubicVerification;
  }

  std::string to_string() const {
    std::ostringstream ss;
    if (solverType == SolverType::STACK)
      ss << "stack-planar ";
    else if (solverType == SolverType::MOVE)
      ss << "move-planar ";
    else
      ss << "brute-force ";

    ss << "[";
    ss << "sat=" << int(useSATConstraints);
    ss << "; unsat=" << int(useUNSATConstraints);
    ss << "; IC=" << int(useIC);
    ss << "; NIC=" << int(useNIC);
    ss << "; strict=" << strict;
    ss << "]";

    return ss.str();
  }
};

enum class ResultCodeTy { SAT, UNSAT, TIMEOUT, ERROR };

struct Result {
  Result() : code(ResultCodeTy::ERROR) {}
  Result(ResultCodeTy code_) : code(code_) {}

  std::string getCodeDesc() const {
    if (code == ResultCodeTy::SAT)
      return "SAT";
    else if (code == ResultCodeTy::UNSAT)
      return "UNSAT";
    else if (code == ResultCodeTy::TIMEOUT)
      return "TIMEOUT";
    else 
      return "ERROR";
  }

  // result code
  ResultCodeTy code;
  // pairs of crossed edges
  std::vector<std::pair<int, int>> crossings;
  // whether the i-th edges is crossed
  std::vector<bool> isCrossed;

  // the order of the original and division (dummy) vertices
  std::vector<std::vector<int>> order;
  // whether i-th segment on stack 0 or 1
  std::vector<bool> stack;
};

struct MVar {
  int id;
  bool positive;

  MVar(): id(0), positive(false) {}
  MVar(int id, bool positive): id(id), positive(positive) {}

  bool operator != (const MVar& other) const {
    if (id != other.id)
      return true;
    return positive != other.positive;
  }

  bool operator < (const MVar& other) const {
    if (id != other.id)
      return id < other.id;
    return positive < other.positive;
  }
};

// Small-buffer vector of MVars. Clauses almost always have 2-3 literals, so an
// inline buffer avoids a heap allocation per clause (the dominant allocation in
// the encoder). Larger clauses spill to the heap. MVar is trivially copyable, so
// growth is a plain memcpy-style copy. Only the API used by MClause is provided.
class MVarVec {
  static constexpr int kInline = 4;
  MVar inln_[kInline];
  MVar* data_ = inln_;
  int size_ = 0;
  int cap_ = kInline;

  void grow(int need) {
    int nc = cap_ * 2;
    while (nc < need) nc *= 2;
    MVar* nd = new MVar[nc];
    for (int i = 0; i < size_; i++) nd[i] = data_[i];
    if (data_ != inln_) delete[] data_;
    data_ = nd;
    cap_ = nc;
  }
  void assign(const MVar* src, int n) {
    if (n > cap_) grow(n);
    for (int i = 0; i < n; i++) data_[i] = src[i];
    size_ = n;
  }

public:
  MVarVec() = default;
  MVarVec(std::initializer_list<MVar> l) { assign(l.begin(), (int)l.size()); }
  MVarVec(const std::vector<MVar>& v) { assign(v.data(), (int)v.size()); }
  MVarVec(const MVarVec& o) { assign(o.data_, o.size_); }
  MVarVec& operator=(const MVarVec& o) { if (this != &o) assign(o.data_, o.size_); return *this; }
  MVarVec(MVarVec&& o) noexcept { *this = std::move(o); }
  MVarVec& operator=(MVarVec&& o) noexcept {
    if (this == &o) return *this;
    if (data_ != inln_) delete[] data_;
    if (o.data_ == o.inln_) {           // other is inline: copy elements
      for (int i = 0; i < o.size_; i++) inln_[i] = o.inln_[i];
      data_ = inln_; cap_ = kInline;
    } else {                            // steal heap buffer
      data_ = o.data_; cap_ = o.cap_;
      o.data_ = o.inln_; o.cap_ = kInline;
    }
    size_ = o.size_; o.size_ = 0;
    return *this;
  }
  ~MVarVec() { if (data_ != inln_) delete[] data_; }

  int size() const { return size_; }
  MVar& operator[](int i) { return data_[i]; }
  const MVar& operator[](int i) const { return data_[i]; }
  MVar* begin() { return data_; }
  MVar* end() { return data_ + size_; }
  const MVar* begin() const { return data_; }
  const MVar* end() const { return data_ + size_; }
  void push_back(const MVar& v) { if (size_ == cap_) grow(size_ + 1); data_[size_++] = v; }
};

struct MClause {
  MVarVec vars;

  MClause() {
  }

  MClause(const MVar& v1) {
    vars.push_back(v1);
  }

  MClause(std::initializer_list<MVar> mvars): vars(mvars) {
  }

  MClause(const std::vector<MVar>& mvars): vars(mvars) {
  }

  void addVar(const MVar& v1) {
    vars.push_back(v1);
  }

  bool operator == (const MClause& other) const {
    return equals(other);
  }

  bool operator < (const MClause& other) const {
    const int size_min = std::min(vars.size(), other.vars.size());
    for (int i = 0; i < size_min; i++) {
      if (vars[i] != other.vars[i]) {
        return vars[i] < other.vars[i];
      }
    }
    return vars.size() < other.vars.size();
  }

  void sort() {
    std::sort(vars.begin(), vars.end());
  }

  // Compute a 64-bit FNV-1a style hash over the clause variables
  size_t hash() const noexcept {
    uint64_t h = 1469598103934665603ULL;
    for (const auto& v : vars) {
      uint64_t x = (uint64_t(v.id) << 1) | (v.positive ? 1ULL : 0ULL);
      h ^= x;
      h *= 1099511628211ULL;
    }
    return (size_t)h;
  }

  // Equality check (same as operator== but as a member)
  bool equals(const MClause& other) const noexcept {
    if (vars.size() != other.vars.size())
      return false;
    for (int i = 0; i < vars.size(); ++i) {
      if (vars[i].id != other.vars[i].id || vars[i].positive != other.vars[i].positive)
        return false;
    }
    return true;
  }
};

// Remove duplicate clauses in a single hash pass, compacting survivors in place.
// The set stores pointers into `vec`, so no clause vector is copied. Only binary
// clauses are ever emitted more than once; larger clauses come from
// already-deduplicated sources, so we neither sort nor hash them (their literal
// order does not affect the CNF). Duplicates are purely a performance concern --
// a repeated clause is still correct -- so any that slip through are harmless.
inline void make_vec_unique_fast(std::vector<MClause>& vec) {
  struct PtrHash {
    size_t operator()(const MClause* c) const noexcept { return c->hash(); }
  };
  struct PtrEq {
    bool operator()(const MClause* a, const MClause* b) const noexcept { return a->equals(*b); }
  };

  const size_t n = vec.size();
  tsl::robin_set<const MClause*, PtrHash, PtrEq> seen;
  seen.reserve(n);
  size_t w = 0;
  for (size_t r = 0; r < n; r++) {
    if (w != r) vec[w] = std::move(vec[r]);
    MClause& c = vec[w];
    if (c.vars.size() != 2) {
      w++;
      continue;
    }
    if (c.vars[1] < c.vars[0])
      std::swap(c.vars[0], c.vars[1]);
    if (seen.insert(&c).second)
      w++;
  }
  vec.resize(w);
}

class SATModel {
  struct pair_hash {
    inline std::size_t operator()(const std::pair<int, int>& v) const {
      std::size_t h1 = std::hash<int>{}(v.first);
      std::size_t h2 = std::hash<int>{}(v.second);
      return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
    }
  };

  std::unordered_map<int, Simp21::Var> vars;
  std::vector<MClause> clauses;
  int curId = 0;

  // relative order node-variables [node_i][node_j]
  std::vector<std::vector<int>> relVars;
  // cross (merged) variables [division_node_index][division_node_index].
  // Flat matrix (like relVars) for O(1) lookup: getCross2Var is called in the
  // O(numVertices^3) merge/K4 loops, so hashing a pair key showed up in profiles.
  // Stored canonically with row < col; -1 means "no variable".
  std::vector<std::vector<int>> cross2Vars;
  // cross variables [division_node_index]
  std::unordered_map<int, int> cross1Vars;
  // move variables [node_index][segment_index]
  std::unordered_map<std::pair<int, int>, int, pair_hash> moveVars;
  // cover variables [edge_index][node_index]
  std::vector<std::vector<int>> coverVars;
  // dir variables [edge_index]
  std::vector<int> dirVars;
  // aux variable
  std::unordered_map<int, int> auxVars;

  // solution (provided by an external solver)
  std::unordered_map<int, bool> externalVars;

  // added 2- and 3-clauses
  ForbiddenTuples forbiddenTuples;

 public:
  SATModel() {
    vars.clear();
    clauses.clear();
    curId = 0;
  }

  int addVar() {
    curId++;
    return curId - 1;
  }

  void addClause(std::initializer_list<MVar> mvars) {
    clauses.emplace_back(mvars);
  }

  void addClause(MClause c) {
    clauses.push_back(std::move(c));
  }

  void reserveClauses(size_t numClauses) {
    clauses.reserve(clauses.size() + numClauses);
  }

  void reserveRelVars(size_t numVertices) {
    relVars = std::vector<std::vector<int>>(numVertices, std::vector<int>(numVertices, -1));
  }

  int addRelVar(int i, int j, bool isSymmetric) {
    CHECK(i != j);
    if (isSymmetric && i > j)
      return -1;
    CHECK(relVars[i][j] == -1);
    relVars[i][j] = addVar();
    return relVars[i][j];
  }

  MVar getRelVar(int i, int j, bool positive) const {
    CHECK(i != j);
    CHECK(relVars[i][j] != -1 || relVars[j][i] != -1, "relVar for tuple (%d, %d) not added", i, j);

    if (relVars[i][j] != -1) {
      return MVar(relVars[i][j], positive);
    } else {
      return MVar(relVars[j][i], !positive);
    }
  }

  MVar getMoveVar(int node, int edge, bool positive) const {
    const auto pair = std::make_pair(node, edge);
    CHECK(moveVars.count(pair));
    const int index = (*moveVars.find(pair)).second;
    return MVar(index, positive);
  }

  void addMoveVar(int node, int edge) {
    const int var = addVar();
    const auto pair = std::make_pair(node, edge);
    CHECK(moveVars.count(pair) == 0);
    moveVars[pair] = var;
  }

  void reserveCross2Vars(size_t numVertices) {
    cross2Vars = std::vector<std::vector<int>>(numVertices, std::vector<int>(numVertices, -1));
  }

  int addCross2Var(int i, int j) {
    CHECK(i != j);
    if (i > j) std::swap(i, j);
    CHECK(cross2Vars[i][j] == -1);
    cross2Vars[i][j] = addVar();
    return cross2Vars[i][j];
  }

  MVar getCross2Var(int i, int j, bool positive) const {
    CHECK(i != j);
    if (i > j) std::swap(i, j);
    CHECK(cross2Vars[i][j] != -1);
    return MVar(cross2Vars[i][j], positive);
  }

  MVar getCross1Var(int edge, bool positive) const {
    CHECK(cross1Vars.count(edge), "cross1Var for edge %d not found", edge);
    return MVar(cross1Vars.at(edge), positive);
  }

  void addCross1Var(int edge) {
    CHECK(!cross1Vars.count(edge));
    cross1Vars[edge] = addVar();
  }

  void reserveCoverVars(size_t numEdges, size_t numVertices) {
    coverVars = std::vector<std::vector<int>>(numEdges, std::vector<int>(numVertices, -1));
    dirVars = std::vector<int>(numEdges, -1);
  }

  int addCoverVar(int e, int v) {
    CHECK(0 <= e && e < (int)coverVars.size());
    CHECK(0 <= v && v < (int)coverVars[e].size());
    CHECK(coverVars[e][v] == -1);
    coverVars[e][v] = addVar();
    return coverVars[e][v];
  }

  bool hasCoverVar(int e, int v) const {
    CHECK(0 <= e && e < (int)coverVars.size());
    CHECK(0 <= v && v < (int)coverVars[e].size());
    return coverVars[e][v] != -1;
  }

  MVar getCoverVar(int e, int v, bool positive) const {
    CHECK(coverVars[e][v] != -1, "coverVars for tuple (%d, %d) not added", e, v);
    return MVar(coverVars[e][v], positive);
  }

  int addDirVar(int e) {
    CHECK(0 <= e && e < (int)dirVars.size());
    CHECK(dirVars[e] == -1);
    dirVars[e] = addVar();
    return dirVars[e];
  }

  MVar getDirVar(int e, bool positive) const {
    CHECK(dirVars[e] != -1, "dirVars for edge %d not added", e);
    return MVar(dirVars[e], positive);
  }

  int addAuxVar() {
    const int var = addVar();
    auxVars[var] = var;
    return var;
  }

  MVar getAuxVar(int id, bool positive) const {
    CHECK(auxVars.find(id) != auxVars.end(), "auxVars for id = %d not added", id);
    return MVar(auxVars.at(id), positive);
  }

  void initVars(Solver& solver) {
    // create vars
    for (int i = 0; i < curId; i++) {
      const auto var = solver.newVar();
      vars[i] = var;
      CHECK(0 <= vars[i] && vars[i] < curId);
    }
    // compress clauses (removes duplicate binary clauses; see make_vec_unique_fast)
    make_vec_unique_fast(clauses);
  }

  Simp21::Lit getSolverLit(const MVar& v) const {
    CHECK(vars.count(v.id));
    const auto var = vars.at(v.id);
    return v.positive ? Simp21::mkLit(var) : ~Simp21::mkLit(var);
  }

  void initClauses(Solver& solver) {
    Simp21::vec<Simp21::Lit> clause;
    for (auto& c : clauses) {
      clause.clear();

      for (const auto& l : c.vars) {
        CHECK(vars.count(l.id));
        const auto var = vars[l.id];

        if (l.positive) {
          clause.push(Simp21::mkLit(var));
        } else {
          clause.push(~Simp21::mkLit(var));
        }
      }

      solver.addClause_(clause);
    }
  }

  void toDimacs(const std::string& filename) {
    std::string ext = filename.substr(filename.find_last_of(".") + 1);

    if (ext == "gz") {
      // zlib
      std::ostringstream oss;
      toDimacs(oss);
      const auto& content = oss.str(); // make it "const auto content = " if there are issues
      gzFile out = gzopen(filename.c_str(), "wb9");
      gzwrite(out, content.c_str(), content.length());
      gzclose(out);
    } else {
      // usual route
      std::ofstream out;
      out.open(filename);
      toDimacs(out);
      out.close();
    }
  }

  void toDimacs(std::ostream& out) {
    const int nvars = varCount();
    out << "p cnf " << nvars << " " << clauseCount() << "\n";

    for (auto& c : clauses) {
      for (auto& l : c.vars) {
        CHECK(vars.count(l.id));
        auto var = vars[l.id] + 1;
        CHECK(1 <= var && var <= nvars);

        if (l.positive) {
          out << var << " ";
        } else {
          out << "-" << var << " ";
        }
      }

      out << "0\n";
    }
  }

  std::string fromDimacs(const std::string& filename) {
    std::ifstream in;
    in.open(filename);
    std::string line;
    std::string externalResult = "";

    while (std::getline(in, line)) {
      std::istringstream iss(line);
      char mode;

      if (!(iss >> mode)) {
        break; // hmm
      }

      if (mode == 's') {
        // result
        iss >> externalResult;
      } else if (mode == 'v') {
        // var
        int vv;

        while (iss >> vv) {
          if (vv > 0) {
            int id = vv - 1;
            CHECK(externalVars.find(id) == externalVars.end());
            externalVars[id] = true;
          } else if (vv < 0) {
            int id = -vv - 1;
            CHECK(externalVars.find(id) == externalVars.end());
            externalVars[id] = false;
          }
        }
      }
    }

    in.close();
    CHECK(externalResult != "");

    if (externalResult == "SATISFIABLE" && externalVars.size() != vars.size()) {
      ERROR("incorrect number of variables in '" + filename + "': " + std::to_string(vars.size()) + " != " + std::to_string(externalVars.size()));
    }

    return externalResult;
  }

  void applyBreakID(int verbose, Solver& solver) {
    // variables are in [1..nvars]; negations are negative
    const int nvars = varCount();
    std::vector<std::vector<int>> cl;

    for (auto& c : clauses) {
      std::vector<int> clause;
      for (auto& l : c.vars) {
        CHECK(vars.count(l.id));
        auto var = vars[l.id] + 1;
        CHECK(1 <= var && var <= nvars);

        if (l.positive) {
          clause.push_back(var);
        } else {
          clause.push_back(-var);
        }
      }
      cl.push_back(clause);
    }

    auto result = breakSymmetries(verbose, nvars, cl);
    LOG_IF(verbose, "introduced %d new variables and %d symmetry-breaking clauses", result.first - nvars, result.second.size());

    // adding vars
    for (int i = nvars; i < result.first; i++) {
      curId++;
      auto var = solver.newVar();
      vars[i] = var;
      CHECK(0 <= vars[i] && vars[i] < curId);
    }
    // adding clauses
    for (auto& c : result.second) {
      MClause clause;
      for (int l : c) {
        CHECK(l != 0);
        if (l > 0) {
          clause.addVar(MVar(l - 1, true));
        } else {
          clause.addVar(MVar(-l - 1, false));
        }
      }
      addClause(clause);

      Simp21::vec<Simp21::Lit> satClause;
      for (auto& l : clause.vars) {
        CHECK(vars.count(l.id));
        auto var = vars[l.id];

        if (l.positive) {
          satClause.push(Simp21::mkLit(var));
        } else {
          satClause.push(~Simp21::mkLit(var));
        }
      }
      solver.addClause_(satClause);
    }
  }

  void applySatsuma(int verbose, Solver& solver) {
    // variables are in [1..nvars]; negations are negative
    const int nvars = varCount();
    std::vector<std::vector<int>> cl;
    cl.reserve(clauses.size());

    std::vector<int> clause;
    for (auto& c : clauses) {
      clause.clear();
      for (auto& l : c.vars) {
        CHECK(vars.count(l.id));
        auto var = vars[l.id] + 1;
        CHECK(1 <= var && var <= nvars);

        if (l.positive) {
          clause.push_back(var);
        } else {
          clause.push_back(-var);
        }
      }
      cl.push_back(clause);
    }

    auto result = applySatsumaSymmetry(verbose, nvars, cl);
    LOG_IF(verbose, "introduced %'d new variables and %'d symmetry-breaking clauses", 
        result.first - nvars, (int)result.second.size() - (int)clauses.size());

    // nothing to add
    if (result.first == nvars && result.second.size() == clauses.size())
      return;

    // adding vars
    for (int i = nvars; i < result.first; i++) {
      curId++;
      auto var = solver.newVar();
      vars[i] = var;
      CHECK(0 <= vars[i] && vars[i] < curId);
    }
    // adding clauses
    clauses.clear();
    for (const auto& c : result.second) {
      MClause clause;
      for (int l : c) {
        CHECK(l != 0);
        if (l > 0) {
          clause.addVar(MVar(l - 1, true));
        } else {
          clause.addVar(MVar(-l - 1, false));
        }
      }
      clauses.push_back(clause);
    }
  }

  bool value(Solver& solver, int id) const {
    CHECK(vars.count(id));
    Simp21::Var var = vars.at(id);

    if (!externalVars.empty()) {
      return externalVars.at(id);
    }

    return (solver.model[var] == l_True ? true : false);
  }

  bool value(Solver& solver, MVar v) const {
    CHECK(vars.count(v.id));
    Simp21::Var var = vars.at(v.id);
    bool positive = v.positive;

    if (!externalVars.empty()) {
      return externalVars.at(v.id) ? positive : !positive;
    }

    return (solver.model[var] == l_True ? positive : !positive);
  }

  size_t varCount() const {
    return vars.size();
  }

  size_t clauseCount() const {
    return clauses.size();
  }

  void fillClauseCount(size_t &clause2Count, size_t &clause3Count, size_t &clause4PlusCount) const {
    clause2Count = 0;
    clause3Count = 0; 
    clause4PlusCount = 0;
    for (const auto& clause : clauses) {
      if (clause.vars.size() == 2)
        clause2Count++;
      else if (clause.vars.size() == 3)
        clause3Count++;
      else 
        clause4PlusCount++;
    }
  }

  ForbiddenTuples& getForbiddenTuples() {
    return forbiddenTuples;
  }
};

bool canBeMerged(int u, int v, const int n, const std::vector<EdgeTy>& edges);
bool canBeMerged(int u, int v, const InputGraph& graph);
void initCrossablePairs(const Params& params, const InputGraph& graph);
void encodeStackPlanar(
    SATModel& model, const InputGraph& graph, const Params& params);
void initSATSolver(
    const Params& params, const InputGraph& graph,
    SATModel& model, Solver& solver);
lbool solveSATModel(
    const Params& params, SATModel& model, Solver& solver);
void fillResultStack(
    const SATModel& model, Solver& solver, const InputGraph& graph,
    const Params& params, Result& result);
Result runSolver(const Params& params, const InputGraph& graph);

void minimizeCrossings(const InputGraph& graph, Result& result, int verbose);
bool isPlanarWithCrossings(const InputGraph& graph, const std::vector<std::pair<int, int>>& crossings);
int computeSkewness(const InputGraph& graph, const int verbose, const int max_skewnees);
Result bruteForce(const Params& params, const InputGraph& graph);

void encodeICConstraints(SATModel& model, const InputGraph& graph, const int verbose, const int C);
void encodeSwapConstraints(SATModel& model, const InputGraph& graph, const Params& params);
void encodeK4Constraints(SATModel& model, const InputGraph& graph, const Params& params);
void encodeCoverConstraints(SATModel& model, const InputGraph& graph, const Params& params);
void encodeStrictConstraints(SATModel& model, const InputGraph& graph, const Params& params);
void encodeStrictConstraintsOption4b(SATModel& model, const InputGraph& graph, const Params& params);
void encodePartialConstraints(SATModel& model, const InputGraph& graph, const Params& params);
void encodeSepCyclesConstraints(SATModel& model, const InputGraph& graph, const Params& params);
std::unique_ptr<Simp21::UserPropagator> createSepCyclesDynamic(
    const SATModel& model,
    const InputGraph& graph,
    const Params& params);
