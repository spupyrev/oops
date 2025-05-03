#include "one_planar.h"
#include "logging.h"

using namespace std;

/// Check if a crossing plus the corresponding K4-induced planar edges violate 4*n-8 density
bool violateDensity(int u, int v, const InputGraph& graph) {
  const int n = graph.n;
  const int m = (int)graph.edges.size();
  CHECK(n >= 4);
  CHECK(u >= n && v >= n);

  const int e_first = graph.edges[u - n].first;
  const int e_second = graph.edges[u - n].second;
  const int f_first = graph.edges[v - n].first;
  const int f_second = graph.edges[v - n].second;
  const std::vector<EdgeTy> k4 = {
    {e_first, f_first},
    {e_first, f_second},
    {e_second, f_first},
    {e_second, f_second}
  };
  int extraEdges = 0;
  for (const auto& [s, t] : k4) {
    if (!graph.hasEdge(s, t)) {
      extraEdges++;
    }
  }
  if (m + extraEdges > 4 * n - 8) {
    return true;
  }
  return false;
}

/// Check if vertices x and a can be swapped so as to eliminate a crossing between edges (x, y) and (a, b)
bool canSwapToReduceCrossings(int x, int y, int a, int b, const InputGraph& graph) {
  const auto& adj = graph.adj;

  std::vector<int> adjX = adj[x];
  CHECK(contains(adj[x], y));
  remove_value(adjX, y);
  std::vector<int> adjA = adj[a];
  CHECK(contains(adj[a], b));
  remove_value(adjA, b);

  // ignore edge (x, a)
  if (contains(adj[x], a)) {
    remove_value(adjX, a);
    remove_value(adjA, x);
  }

  if (equal_unsorted(adjX, adjA))
    return true;

  // one extra edge is possible
  const bool hasEdgeAY = contains(adj[a], y);
  const bool hasEdgeXB = contains(adj[x], b);

  // If one of the edges is missing, then we can always reduce the crossings
  if (hasEdgeAY && !hasEdgeXB) {
    remove_value(adjA, y);
    if (equal_unsorted(adjX, adjA)) {
      return true;    
    }
  }
  if (!hasEdgeAY && hasEdgeXB) {
    remove_value(adjX, b);
    if (equal_unsorted(adjX, adjA)) {
      return true;    
    }
  }

  // If both edges are present, then we require the smaller pair to be crossing-free
  // TODO: disabled for now to avoid interactions with twin ordering
  // if (hasEdgeAY && hasEdgeXB) {
  //   remove_value(adjA, y);
  //   remove_value(adjX, b);
  //   if (equal_unsorted(adjX, adjA)) {
  //     // possible crossings: (x, y)--(a, b) and (x, b)--(a, y)
  //     const int e1 = graph.findDivIndex(x, y);
  //     const int e2 = graph.findDivIndex(a, b);
  //     const int e3 = graph.findDivIndex(x, b);
  //     const int e4 = graph.findDivIndex(a, y);
  //     if (std::minmax(e1, e2) < std::minmax(e3, e4)) {
  //       // LOG("disable a crossing for (x=%d, y=%d) and (a=%d, b=%d); adjX = %s; adjA = %s", x, y, a, b, 
  //       //   to_string(adjX).c_str(), to_string(adjA).c_str());
  //       return true;
  //     }
  //   }
  // }

  return false;
}

// Pairs of edges that could be crossed; equivalently, pairs of division vertices that can be merged
std::vector<std::vector<bool>> crossablePairs;

/// Find all pairs of edges that can be crossed
void initCrossablePairs(const Params& params, const InputGraph& graph) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  crossablePairs = std::vector<std::vector<bool>>(numVertices, std::vector<bool>(numVertices, false));

  if (params.forbidCrossings) {
    CHECK(graph.isDirected());
    return;
  }

  // Skip adjacent edges and density violations
  int numDensitySkipped = 0;
  for (int u = n; u < numVertices; u++) {
    for (int v = u + 1; v < numVertices; v++) {
      // two (non-adjacent) division vertices can be non-comparable
      if (graph.adjacent(edges[u - n], edges[v - n])) {
        continue;
      }
      // some crossings violate density
      if (violateDensity(u, v, graph)) {
        numDensitySkipped++;
        continue;
      }
      crossablePairs[u][v] = true;
      crossablePairs[v][u] = true;
    }
  }

  const auto& adj = graph.adj;
  
  // Disable certain degree-3 crossings
  int numDegree3Skipped = 0;
  if (!graph.isDirected()) {
    for (int v = 0; v < n; v++) {
      if (adj[v].size() > 3)
        continue;
      for (int x : adj[v]) {
        const int e1 = graph.findDivIndex(v, x);
        for (int e2 = n; e2 < numVertices; e2++) {
          if (!crossablePairs[e1][e2])
            continue;
          const int y = edges[e2 - n].first;
          const int z = edges[e2 - n].second;
          const bool allVEdgesToXYZ = std::all_of(adj[v].begin(), adj[v].end(), [&](int t) { return t == x || t == y || t == z; });
          if (!allVEdgesToXYZ)
            continue;

          CHECK(e1 != e2);
          crossablePairs[e1][e2] = false;
          crossablePairs[e2][e1] = false;
          numDegree3Skipped++;
        }
      }
    }
  }

  // Disable crossings for almost-twins
  int numAlmostTwinsSkipped = 0;
  if (!graph.isDirected()) {
    for (int u = n; u < numVertices; u++) {
      for (int v = u + 1; v < numVertices; v++) {
        if (!crossablePairs[u][v])
          continue;        
        // edges are (x, y) and (a, b)
        const int x = edges[u - n].first;
        const int y = edges[u - n].second;
        const int a = edges[v - n].first;
        const int b = edges[v - n].second;

        if (canSwapToReduceCrossings(x, y, a, b, graph) || canSwapToReduceCrossings(x, y, b, a, graph) || 
            canSwapToReduceCrossings(y, x, a, b, graph) || canSwapToReduceCrossings(y, x, b, a, graph)) {
          crossablePairs[u][v] = false;
          crossablePairs[v][u] = false;
          numAlmostTwinsSkipped++;
        }
      }
    }
  }


  // //////////////////////////////////////////////////////////
  // // Count distances
  // const int INF = 123456789;
  // std::vector<std::vector<int>> dist(n, std::vector<int>(n, INF));
  // for (int i = 0; i < n; i++) {
  //   dist[i][i] = 0;
  // }
  // for (const auto& [u, v] : edges) {
  //   dist[u][v] = 1;
  //   dist[v][u] = 1;
  // }
  // for (int k = 0; k < n; k++) {
  //   for (int i = 0; i < n; i++) {
  //     for (int j = 0; j < n; j++) {
  //       if (dist[i][j] > dist[i][k] + dist[k][j])
  //         dist[i][j] = dist[i][k] + dist[k][j];
  //     }
  //   }
  // }
  // int num0 = 0;
  // int num1 = 0;
  // int num2 = 0;
  // int num3 = 0;
  // int num4 = 0;
  // int numSP = 0;
  // for (int e1 = 0; e1 < m; e1++) {
  //   for (int e2 = e1 + 1; e2 < m; e2++) {
  //     const int u1 = edges[e1].first;
  //     const int v1 = edges[e1].second;
  //     const int u2 = edges[e2].first;
  //     const int v2 = edges[e2].second;
  //     const int d = std::min(
  //       std::min(dist[u1][u2], dist[u1][v2]), 
  //       std::min(dist[v1][u2], dist[v1][v2])
  //     );
  //     const int D = std::max(
  //       std::max(dist[u1][u2], dist[u1][v2]), 
  //       std::max(dist[v1][u2], dist[v1][v2])
  //     );
  //     if (d == 0)
  //       num0++;
  //     if (d == 1)
  //       num1++;
  //     if (d == 2)
  //       num2++;
  //     if (d == 3)
  //       num3++;
  //     if (d == 4)
  //       num4++;
  //     if (d > 0 && d + 2 == D) {
  //       numSP++;
  //       // crossablePairs[e1 + n][e2 + n] = false;
  //       // crossablePairs[e2 + n][e1 + n] = false;
  //     }
  //   }
  // }
  // // LOG("num0 = %d", num0);
  // // LOG("num1 = %d", num1);
  // // LOG("num2 = %d", num2);
  // // LOG("num3 = %d", num3);
  // // LOG("num4 = %d", num4);
  // LOG_IF(verbose, "numSP = %d", numSP);

  // //////////////////////////////////////////////////////////


  // Count pairs
  int mergablePairs = 0;
  int possiblePairs = 0;
  for (int u = n; u < numVertices; u++) {
    for (int v = u + 1; v < numVertices; v++) {
      possiblePairs++;
      if (crossablePairs[u][v])
        mergablePairs++;
    }
  }
  LOG_IF(params.verbose, "created %d (%.2lf%% out of %d) crossing pairs; "
                         "filtered out %d (%.2lf%%) due to density, %d (%.2lf%%) due to degree-3, and %d (%.2lf%%) almost twins", 
    mergablePairs, 100.0 * mergablePairs / double(possiblePairs), possiblePairs,
    numDensitySkipped, 100.0 * numDensitySkipped / double(possiblePairs),
    numDegree3Skipped, 100.0 * numDegree3Skipped / double(possiblePairs),
    numAlmostTwinsSkipped, 100.0 * numAlmostTwinsSkipped / double(possiblePairs)
  );
}

/// Return true iff the two vertices can be merged
bool canBeMerged(int u, int v, const int n, const std::vector<EdgeTy>& edges) {
  CHECK(u != v);
  CHECK(u >= n && v >= n);
  CHECK((int)crossablePairs.size() == n + (int)edges.size());
  return crossablePairs[u][v];
}

/// Encode relative position constraints for pairs of vertices
void encodeRelativeVariables(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  // Create variables for every pair of vertices
  model.reserveRelVars(numVertices);
  for (int i = 0; i < numVertices; i++) {
    for (int j = 0; j < numVertices; j++) {
      if (i == j)
        continue;
      
      if (i < n || j < n)
        model.addRelVar(i, j, true);
      else
        model.addRelVar(i, j, !canBeMerged(i, j, n, edges));
    }
  }

  // Ensure associativity for non-mergable vertices
  model.reserveClauses(numVertices * numVertices * numVertices + numVertices * numVertices);
  for (int i = n; i < numVertices; i++) {
    for (int j = n + 1; j < numVertices; j++) {
      if (i == j)
        continue;
      if (canBeMerged(i, j, n, edges)) {
        // !(i < j) or !(j < i)
        model.addClause({
            model.getRelVar(i, j, false), 
            model.getRelVar(j, i, false)
        });
      }
    }
  }

  // Ensure transitivity
  for (int i = 0; i < numVertices; i++) {
    for (int j = 0; j < numVertices; j++) {
      for (int k = 0; k < numVertices; k++) {
        if (i == j || i == k || j == k)
          continue;

        // i <= j && j <= k => i < k
        model.addClause(MClause({
            model.getRelVar(j, i, true), 
            model.getRelVar(k, j, true), 
            model.getRelVar(i, k, true)
        }));
      }
    }
  }

  // Ensure edge directions
  if (graph.isDirected()) {
    for (size_t i = 0; i < edges.size(); i++) {
      const int s = graph.directions[i] ? edges[i].first : edges[i].second;
      const int t = graph.directions[i] ? edges[i].second : edges[i].first;
      const int div = n + i;
      model.addClause(MClause(model.getRelVar(s, t, true)));
      model.addClause(MClause(model.getRelVar(s, div, true)));
      model.addClause(MClause(model.getRelVar(div, t, true)));
    }
  }
}

void encodeICConstraints(SATModel& model, const InputGraph& graph, const int verbose, const int C);
void encodeSwapConstraints(SATModel& model, const InputGraph& graph, const int verbose);
void encodeK4Constraints(SATModel& model, const InputGraph& graph, const int verbose);
void encodeCoverConstraints(SATModel& model, const InputGraph& graph, const int verbose);

/// Add pairwise crossing (merge) variables
void encodeCross2Variables(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  // Pairwise "merged" variables 
  for (int i = n; i < numVertices; i++) {
    for (int j = i + 1; j < numVertices; j++) {
      if (canBeMerged(i, j, n, edges)) {
        model.addCross2Var(i, j);
        // Sync cross variables with relative variables
        model.addClause({
            model.getRelVar(i, j, true), 
            model.getRelVar(j, i, true), 
            model.getCross2Var(i, j, true) 
        });
        model.addClause({
            model.getRelVar(i, j, false), 
            model.getCross2Var(i, j, false) 
        });
        model.addClause({
            model.getRelVar(j, i, false), 
            model.getCross2Var(i, j, false) 
        });
      }
    }
  }
  // At most two division vertices are merged together
  for (int i = n; i < numVertices; i++) {
    for (int j = i + 1; j < numVertices; j++) {
      for (int k = j + 1; k < numVertices; k++) {        
        if (canBeMerged(i, j, n, edges) && canBeMerged(i, k, n, edges)) {
          model.addClause({
              model.getCross2Var(i, j, false), 
              model.getCross2Var(i, k, false)
          });
        }
        if (canBeMerged(i, j, n, edges) && canBeMerged(j, k, n, edges)) {
          model.addClause({
              model.getCross2Var(i, j, false), 
              model.getCross2Var(j, k, false)
          });
        }
        if (canBeMerged(i, k, n, edges) && canBeMerged(j, k, n, edges)) {
          model.addClause({
              model.getCross2Var(i, k, false), 
              model.getCross2Var(j, k, false)
          });
        }
      }
    }
  }

  if (!graph.isDirected()) {
    // TODO: always enable for UNSAT?
    encodeICConstraints(model, graph, verbose, 2);

    // TODO: always enable for UNSAT?
    encodeSwapConstraints(model, graph, verbose);

    // TODO: always enable for UNSAT?
    encodeK4Constraints(model, graph, verbose);
  }
}

/// Encode edge-crossing variables and constraints
void encodeCross1Constraints(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  // Variables
  for (int i = n; i < numVertices; i++) {
    model.addCross1Var(i);
  }
  // A pair crosses => the two edges cross
  for (int i = n; i < numVertices; i++) {
    for (int j = i + 1; j < numVertices; j++) {
      if (!canBeMerged(i, j, n, edges)) 
        continue;

      model.addClause({model.getCross2Var(i, j, false), model.getCross1Var(i, true)});
      model.addClause({model.getCross2Var(i, j, false), model.getCross1Var(j, true)});
    }
  }
  // An edge crosses => at least one pairwise crossing
  for (int i = n; i < numVertices; i++) {
    MClause clause(model.getCross1Var(i, false));
    for (int j = n; j < numVertices; j++) {
      if ((i < j && canBeMerged(i, j, n, edges)) || (j < i && canBeMerged(j, i, n, edges))) {
        clause.addVar(model.getCross2Var(i, j, true));
      }
    }
    model.addClause(clause);
  }

  // Lexicographic crossings on degree-2
  const auto& adj = graph.adj;
  int numDegree2Lex = 0;
  for (int v = 0; v < n; v++) {
    if (adj[v].size() != 2)
      continue;
    int x = adj[v][0];
    int y = adj[v][1];
    if (x > y) std::swap(x, y);
    const int e1 = graph.findDivIndex(v, x);
    const int e2 = graph.findDivIndex(v, y);

    model.addClause({
      model.getCross1Var(e1, false), 
      model.getCross1Var(e2, true)
    });
    numDegree2Lex++;
  }

  // A pair cross => K4-edges do not cross
  int numPlanarK4 = 0;
  for (int i = n; i < numVertices; i++) {
    for (int j = i + 1; j < numVertices; j++) {
      if (!canBeMerged(i, j, n, edges)) 
        continue;

      const int e_first = edges[i - n].first;
      const int e_second = edges[i - n].second;
      const int f_first = edges[j - n].first;
      const int f_second = edges[j - n].second;
      std::vector<EdgeTy> k4 = {
        {e_first, f_first},
        {e_first, f_second},
        {e_second, f_first},
        {e_second, f_second}
      };
      for (const auto& [u, v] : k4) {
        const int x = graph.findEdgeIndex(u, v);
        if (x != -1) {
          model.addClause({model.getCross2Var(i, j, false), model.getCross1Var(x + n, false)});
          numPlanarK4++;
        }
      }
    }
  }

  LOG_IF(verbose, "cross-1 constraints: %d planar-K4; %d degree-2-lex", numPlanarK4, numDegree2Lex);

  if (!graph.isDirected()) {
    encodeCoverConstraints(model, graph, verbose);
  }
}

/// Encode cover variables and constraints
void encodeCoverConstraints(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  // Variables
  model.reserveCoverVars(m, numVertices);
  for (int i = n; i < numVertices; i++) {
    const int e = i - n;
    model.addCoverVar(e, i);

    const int u = edges[e].first;
    const int v = edges[e].second;

    // Cover is true => vertex is in between
    model.addClause({
        model.getCoverVar(e, i, false), 
        model.getRelVar(i, u, true),
        model.getRelVar(i, v, true)
    });
    model.addClause({
        model.getCoverVar(e, i, false), 
        model.getRelVar(u, i, true),
        model.getRelVar(v, i, true)
    });

    // // Cover is false => vertex is not in between (optional)
    // model.addClause({
    //     model.getCoverVar(e, i, true), 
    //     model.getRelVar(u, v, false),
    //     model.getRelVar(i, u, true),
    //     model.getRelVar(v, i, true)
    // });
    // model.addClause({
    //     model.getCoverVar(e, i, true), 
    //     model.getRelVar(u, v, true),
    //     model.getRelVar(i, v, true),
    //     model.getRelVar(u, i, true)
    // });
  }

  // Constraints
  // An edge doesn't cross => the division is in between
  int numMid = 0;
  for (int i = n; i < numVertices; i++) {
    const int e = i - n;
    model.addClause({
      model.getCross1Var(i, true), 
      model.getCoverVar(e, i, true)
    });
    numMid++;
  }
  // Every crossed division vertex is in between
  int numMid2 = 0;
  for (int i = n; i < numVertices; i++) {
    for (int j = i + 1; j < numVertices; j++) {
      if (!canBeMerged(i, j, n, edges)) 
        continue;

      const int e1 = i - n;
      const int e2 = j - n;

      model.addClause({
        model.getCross2Var(i, j, false), 
        model.getCoverVar(e1, i, true),
        model.getCoverVar(e2, j, true)
      });
      numMid2++;
    }
  }

  LOG_IF(verbose, "cover constraints: %d mid-1; %d mid-2", numMid, numMid2);
}

/// C=0 <=> IC; C=1 <=> NIC; C=2 <=> 1-planar
void encodeICConstraints(SATModel& model, const InputGraph& graph, const int verbose, const int C) {
  CHECK(0 <= C && C <= 3);

  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  // Two pairs of crossings have at most C vertices in common
  int numIC = 0;
  for (int i1 = n; i1 < numVertices; i1++) {
    for (int j1 = i1 + 1; j1 < numVertices; j1++) {
      if (!canBeMerged(i1, j1, n, edges))
        continue;
      // a crossing for edges (u, v) and (x, y)
      const int u1 = edges[i1 - n].first;
      const int v1 = edges[i1 - n].second;
      const int x1 = edges[j1 - n].first;
      const int y1 = edges[j1 - n].second;
      const std::vector<int> k4_vertices = {u1, v1, x1, y1};

      for (int i2 = i1 + 1; i2 < numVertices; i2++) {
        for (int j2 = i2 + 1; j2 < numVertices; j2++) {
          if (!canBeMerged(i2, j2, n, edges))
            continue;
          if (i1 == i2 || i1 == j2 || j1 == i2 || j1 == j2)
            continue;
          const int u2 = edges[i2 - n].first;
          const int v2 = edges[i2 - n].second;
          const int x2 = edges[j2 - n].first;
          const int y2 = edges[j2 - n].second;

          int numCommon = 0;
          if (contains(k4_vertices, u2)) numCommon++;
          if (contains(k4_vertices, v2)) numCommon++;
          if (contains(k4_vertices, x2)) numCommon++;
          if (contains(k4_vertices, y2)) numCommon++;

          CHECK(i1 < j1 && i2 < j2 && i1 < i2);
          if (numCommon > C) {
            model.addClause({
                model.getCross2Var(i1, j1, false),
                model.getCross2Var(i2, j2, false)
            });
            numIC++;
            // LOG("i1 = %2d; j1 = %2d; i2 = %2d; j2 = %2d  [numCommon = %d]", i1, j1, i2, j2, numCommon);
          }
        }
      }
    }
  }

  LOG_IF(verbose, "added %2d IC-%d constraints", numIC, C);
}

/// Disable pairs of crossings that can be eliminated by swapping pairs of vertices
void encodeSwapConstraints(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const auto& adj = graph.adj;

  int numConstraints = 0;
  // try to eliminate crossings (x, b)--(y, a) and (x, v)--(y, u) by exchanging x and y
  for (int x = 0; x < n; x++) {
    for (int y = x + 1; y < n; y++) {
      for (size_t ib = 0; ib < adj[x].size(); ib++) {
        for (size_t iv = 0; iv < adj[x].size(); iv++) {
          if (ib == iv) continue;
          const int b = adj[x][ib];
          const int v = adj[x][iv];
          if (b == y || v == y) continue;
          CHECK(b != v);

          const int be = graph.findDivIndex(b, x);
          const int ve = graph.findDivIndex(v, x);
          for (size_t ia = 0; ia < adj[y].size(); ia++) {
            for (size_t iu = 0; iu < adj[y].size(); iu++) {
              if (ia == iu) continue;
              const int a = adj[y][ia];
              const int u = adj[y][iu];
              if (a == x || u == x) continue;
              CHECK(a != u);

              const int ae = graph.findDivIndex(a, y);
              const int ue = graph.findDivIndex(u, y);

              if (b == a || b == u || v == u || v == a)
                continue;

              CHECK(be != ae && ve != ue);
              if (!canBeMerged(be, ae, n, edges))
                continue;
              if (!canBeMerged(ve, ue, n, edges))
                continue;
              CHECK(b != a && b != u && v != u && v != a);

              // vertices that can be reached from x w/o violating 1-planarity
              std::vector<int> reachX = adj[x];
              CHECK(contains(reachX, b) && contains(reachX, b));
              reachX.push_back(x);
              reachX.push_back(y);
              reachX.push_back(a);
              reachX.push_back(u);
              // vertices that can be reached from y w/o violating 1-planarity
              std::vector<int> reachY = adj[y];
              CHECK(contains(reachY, a) && contains(reachY, u));
              reachY.push_back(y);
              reachY.push_back(x);
              reachY.push_back(b);
              reachY.push_back(v);

              if (!equal_unsorted(reachX, reachY))
                continue;

              if (std::minmax(ae, be) < std::minmax(ve, ue)) {
                //LOG("swap crossing0: x = %d; y = %d; a = %d; b = %d; u = %d; v = %d", x, y, a, b, u, v);
                model.addClause({
                    model.getCross2Var(ae, be, false), 
                    model.getCross2Var(ve, ue, false)
                });
                numConstraints++;
              }
            }
          }
        }
      }
    }
  }

  LOG_IF(verbose, "added %3d swap constraints", numConstraints);

  int num3Constraints = 0;
  // forbid crossings between every pair of edges for two degree-3 vertices
  for (int u = 0; u < n; u++) {
    if (graph.degree(u) != 3) 
      continue;
    for (int v = u + 1; v < n; v++) {
      if (graph.degree(v) != 3) 
        continue;
      if (graph.hasEdge(u, v))
        continue;

      // vertex u is adjacent to (a, b, c) via edges (ua, ub, uc)
      const int a = adj[u][0];
      const int ua = graph.findDivIndex(u, a);
      const int b = adj[u][1];
      const int ub = graph.findDivIndex(u, b);
      const int c = adj[u][2];
      const int uc = graph.findDivIndex(u, c);
      vector<int> perm = identity(3);
        // vertex v is adjacent to (x, y, z) via edges (vx, vy, vz)
        const int x = adj[v][perm[0]];
        const int vx = graph.findDivIndex(v, x);
        const int y = adj[v][perm[1]];
        const int vy = graph.findDivIndex(v, y);
        const int z = adj[v][perm[2]];
        const int vz = graph.findDivIndex(v, z);

        if (!canBeMerged(ua, vx, n, edges))
          continue;
        if (!canBeMerged(ub, vy, n, edges))
          continue;
        if (!canBeMerged(uc, vz, n, edges))
          continue;

        LOG_IF(verbose >= 3, "swap-3 crossing: (%d--%d); (%d--%d); (%d--%d)", ua, vx, ub, vy, uc, vz);
        model.addClause({
            model.getCross2Var(ua, vx, false), 
            model.getCross2Var(ub, vy, false),
            model.getCross2Var(uc, vz, false)
        });
        num3Constraints++;
      do {
      } while (std::next_permutation(perm.begin(), perm.end()));
    }
  }

  LOG_IF(verbose, "added %3d 3-swap constraints", num3Constraints);
}

/// A pair cross => K4-edges do not cross
void encodeK4Constraints(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  const auto& adj = graph.adj;

  int numConstraints = 0;
  for (int divUV = n; divUV < numVertices; divUV++) {
    for (int divXY = n; divXY < numVertices; divXY++) {
      if (divUV == divXY)
        continue;
      if (!canBeMerged(divUV, divXY, n, edges))
        continue;
      // a crossing for edges (u, v) and (x, y)
      const int u = edges[divUV - n].first;
      const int v = edges[divUV - n].second;

      for (size_t iu = 0; iu < adj[u].size(); iu++) {
        const int a = adj[u][iu];
        if (a == v) continue;
        const int divUA = graph.findDivIndex(u, a);
        for (size_t iv = 0; iv < adj[v].size(); iv++) {
          const int b = adj[v][iv];
          if (b == u) continue;
          const int divVB = graph.findDivIndex(v, b);
          if (!canBeMerged(divUA, divVB, n, edges))
            continue;
          CHECK(divUV != divVB && divUV != divUA);
          if (std::minmax(divUV, divXY) < std::minmax(divVB, divUA)) {
            model.addClause({
                model.getCross2Var(divUV, divXY, false),
                model.getCross2Var(divVB, divUA, false)
            });
            numConstraints++;
          }
        }
      }
    }
  }

  LOG_IF(verbose, "added %2d K4 constraints", numConstraints);
}

void encodeMoveVariables(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  const int numSegments = 2 * m;

  // Create "move" variables for vertex->segment
  for (int v = 0; v < numVertices; v++) {
    for (int seg = 0; seg < numSegments; seg++) {
      // skip adjacent
      const auto [e_first, e_second] = graph.seg2edge(seg);
      if (e_first == v || e_second == v)
        continue;

      model.addMoveVar(v, seg);
    }
  }

  // Sync move variables for merged vertices
  for (int seg = 0; seg < numSegments; seg++) {
    const auto [e_first, e_second] = graph.seg2edge(seg);
    for (int v1 = n; v1 < numVertices; v1++) {
      for (int v2 = v1 + 1; v2 < numVertices; v2++) {
        if (e_first == v1 || e_second == v1 || e_first == v2 || e_second == v2)
          continue;
        if (!canBeMerged(v1, v2, n, edges))
          continue;

        // merged(i, j) => move(v1, seg) == move(v2, seg)
        model.addClause(MClause({
            model.getCross2Var(v1, v2, false), 
            model.getMoveVar(v1, seg, true),
            model.getMoveVar(v2, seg, false)
        }));
        model.addClause(MClause({
            model.getCross2Var(v1, v2, false), 
            model.getMoveVar(v1, seg, false),
            model.getMoveVar(v2, seg, true)
        }));
      }
    }
  }

}

void encodeMoveConstraints(SATModel& model, const InputGraph& graph, const int verbose) {
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numSegments = 2 * m;

  auto rm0Clause = [&](int seg_e, int seg_f, int se, int te, int sf, int tf, bool moveFlag) {
    return MClause({
        model.getRelVar(se, sf, false),
        model.getRelVar(sf, te, false),
        model.getRelVar(te, tf, false),
        model.getMoveVar(sf, seg_e, moveFlag),
        model.getMoveVar(te, seg_f, moveFlag)
    });
  };
  auto rm1Clause = [&](int seg_e, int seg_f, int se, int te, int sf, int tf, bool moveFlag) {
    return MClause({
        model.getRelVar(se, sf, false),
        model.getRelVar(sf, tf, false),
        model.getRelVar(tf, te, false),
        model.getMoveVar(sf, seg_e, !moveFlag),
        model.getMoveVar(tf, seg_e, moveFlag)
    });
  };

  for (int seg_e = 0; seg_e < numSegments; seg_e++) {
    for (int seg_f = 0; seg_f < numSegments; seg_f++) {
      if (seg_e == seg_f) 
        continue;

      // edge e
      const auto [e_first, e_second] = graph.seg2edge(seg_e);
      // edge f
      const auto [f_first, f_second] = graph.seg2edge(seg_f);

      // skip adjacent edge
      if (graph.adjacent(EdgeTy(e_first, e_second), EdgeTy(f_first, f_second)))
        continue;

      // R^m_0
      // Assume (SE=e.first < TE=e.second) && (SF=f.first < TF=f.second)
      model.addClause(rm0Clause(seg_e, seg_f, e_first, e_second, f_first, f_second, false));
      model.addClause(rm0Clause(seg_e, seg_f, e_first, e_second, f_first, f_second, true));
      // Assume (SE=e.second < TE=e.first) && (SF=f.first < TF=f.second)
      model.addClause(rm0Clause(seg_e, seg_f, e_second, e_first, f_first, f_second, false));
      model.addClause(rm0Clause(seg_e, seg_f, e_second, e_first, f_first, f_second, true));
      // Assume (SE=e.first < TE=e.second) && (SF=f.second < TF=f.first)
      model.addClause(rm0Clause(seg_e, seg_f, e_first, e_second, f_second, f_first, false));
      model.addClause(rm0Clause(seg_e, seg_f, e_first, e_second, f_second, f_first, true));
      // Assume (SE=e.second < TE=e.first) && (SF=f.second < TF=f.first)
      model.addClause(rm0Clause(seg_e, seg_f, e_second, e_first, f_second, f_first, false));
      model.addClause(rm0Clause(seg_e, seg_f, e_second, e_first, f_second, f_first, true));

      // R^m_1
      // Assume (SE=e.first < TE=e.second) && (SF=f.first < TF=f.second)
      model.addClause(rm1Clause(seg_e, seg_f, e_first, e_second, f_first, f_second, false));
      model.addClause(rm1Clause(seg_e, seg_f, e_first, e_second, f_first, f_second, true));
      // Assume (SE=e.second < TE=e.first) && (SF=f.first < TF=f.second)
      model.addClause(rm1Clause(seg_e, seg_f, e_second, e_first, f_first, f_second, false));
      model.addClause(rm1Clause(seg_e, seg_f, e_second, e_first, f_first, f_second, true));
      // Assume (SE=e.first < TE=e.second) && (SF=f.second < TF=f.first)
      model.addClause(rm1Clause(seg_e, seg_f, e_first, e_second, f_second, f_first, false));
      model.addClause(rm1Clause(seg_e, seg_f, e_first, e_second, f_second, f_first, true));
      // Assume (SE=e.second < TE=e.first) && (SF=f.second < TF=f.first)
      model.addClause(rm1Clause(seg_e, seg_f, e_second, e_first, f_second, f_first, false));
      model.addClause(rm1Clause(seg_e, seg_f, e_second, e_first, f_second, f_first, true));
    }
  }
}

void encodeMoveSymmetry(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  int numExtraConstraints = 0;

  if (graph.isDirected())
    return;

  // the order starts at 0
  for (int i = 1; i < n; i++) {
    model.addClause(MClause(model.getRelVar(0, i, true)));
    numExtraConstraints++;
  }

  // fix one move variable
  model.addClause(MClause(model.getMoveVar(0, 0, false)));
  numExtraConstraints++;

  // vertex twins
  const auto& adj = graph.adj;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      std::vector<int> adjI = adj[i];
      adjI.push_back(i);
      adjI.push_back(j);
      sort_unique(adjI);
      std::vector<int> adjJ = adj[j];
      adjJ.push_back(i);
      adjJ.push_back(j);
      sort_unique(adjJ);
      if (adjI == adjJ) {
        LOG_IF(verbose >= 3, "identified equal adjacencies for vertices %d and %d: (%s) vs (%s)", i, j, to_string(adj[i]).c_str(), to_string(adj[j]).c_str());
        model.addClause(MClause(model.getRelVar(i, j, true)));
        numExtraConstraints++;
      }
    }
  }

  LOG_IF(verbose && numExtraConstraints > 0, "added %d symmetry-breaking constraints", numExtraConstraints);
}

/// Encode move-planarity
void encodeMovePlanar(
    SATModel& model, 
    const InputGraph& graph,
    const Params& params) {
  LOG_IF(params.verbose, "encoding %s", params.to_string().c_str()); 
  CHECK(params.useSATConstraints, "move-planarity should be used with -cross2");

  // Main encoding
  encodeRelativeVariables(model, graph, params.verbose);
  encodeCross2Variables(model, graph, params.verbose);
  encodeMoveVariables(model, graph, params.verbose);
  encodeMoveConstraints(model, graph, params.verbose);

  // Optional encoding
  if (params.useUNSATConstraints) {
    CHECK(params.useSATConstraints);
    encodeCross1Constraints(model, graph, params.verbose);
  }
  if (params.useIC) {
    CHECK(params.useSATConstraints);
    encodeICConstraints(model, graph, params.verbose, 0);
  }
  if (params.useNIC) {
    CHECK(params.useSATConstraints);
    encodeICConstraints(model, graph, params.verbose, 1);
  }

  // Symmetry
  encodeMoveSymmetry(model, graph, params.verbose);
}

void encodeStackConstraints(SATModel& model, const InputGraph& graph, const int verbose) {
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numSegments = 2 * m;

  // forbid order a < b < c < d
  auto crossClause = [&](int a, int b, int c, int d) {
    return MClause({
        model.getRelVar(a, b, false),
        model.getRelVar(b, c, false),
        model.getRelVar(c, d, false)
    });
  };

  for (int seg_e = 0; seg_e < numSegments; seg_e++) {
    for (int seg_f = seg_e + 1; seg_f < numSegments; seg_f++) {
      const auto [e_org, e_div] = graph.seg2edge_v2(seg_e);
      const auto [f_org, f_div] = graph.seg2edge_v2(seg_f);
      // skip adjacent edges
      if (graph.adjacent(EdgeTy(e_org, e_div), EdgeTy(f_org, f_div)))
        continue;

      // e_org < f_org < e_div < f_div [stack 0]
      model.addClause(crossClause(e_org, f_org, e_div, f_div));
      // f_org < e_org < f_div < e_div [stack 0]
      model.addClause(crossClause(f_org, e_org, f_div, e_div));
      // e_div < f_div < e_org < f_org [stack 1]
      model.addClause(crossClause(e_div, f_div, e_org, f_org));
      // f_div < e_div < f_org < e_org [stack 1]
      model.addClause(crossClause(f_div, e_div, f_org, e_org));
    }
  }
}

void encodeStackSymmetry(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  int numExtraConstraints = 0;

  // the order starts at 0
  for (int i = 1; i < numVertices; i++) {
    model.addClause(MClause(model.getRelVar(0, i, true)));
    numExtraConstraints++;
  }

  /////////////////////////////////////////////////////////////////////////////////
  // // TMP: custom constraints
  // for (int i = 0; i < numVertices; i++) {
  //   if (6 != i)
  //     model.addClause(MClause(model.getRelVar(6, i, true)));
  //   numExtraConstraints++;
  // }
  // model.addClause(MClause(model.getCross2Var(graph.findDivIndex(4, 10), graph.findDivIndex(5, 9), true)));
  // model.addClause(MClause(model.getCross2Var(graph.findDivIndex(1, 7), graph.findDivIndex(2, 6), true)));
  // LOG(TextColor::red, "encoding contains custom constraints!!!");
  // for (int i = 0; i < m; i++) {
  //   LOG("edge (%d, %d): %d", edges[i].first, edges[i].second, graph.findDivIndex(edges[i].first, edges[i].second));
  // }
  /////////////////////////////////////////////////////////////////////////////////

  // vertex twins
  const auto& adj = graph.adj;
  for (int i = 0; i < n; i++) {
    for (int j = i + 1; j < n; j++) {
      std::vector<int> adjI = adj[i];
      adjI.push_back(i);
      adjI.push_back(j);
      std::vector<int> adjJ = adj[j];
      adjJ.push_back(i);
      adjJ.push_back(j);
      if (equal_unsorted(adjI, adjJ)) {
        LOG_IF(verbose >= 3, "identified equal adjacencies for vertices %d and %d: (%s) vs (%s)", i, j, to_string(adj[i]).c_str(), to_string(adj[j]).c_str());
        model.addClause(MClause(model.getRelVar(i, j, true)));
        numExtraConstraints++;
      }
    }
  }

  LOG_IF(verbose && numExtraConstraints > 0, "added %d symmetry-breaking constraints", numExtraConstraints);
}

/// Encode stack-planarity
void encodeStackPlanar(
    SATModel& model, 
    const InputGraph& graph,
    const Params& params) {
  LOG_IF(params.verbose, "encoding %s", params.to_string().c_str()); 
  CHECK(!graph.isDirected(), "directed edges should be used with move-planarity");

  // Main encoding
  encodeRelativeVariables(model, graph, params.verbose);
  encodeStackConstraints(model, graph, params.verbose);

  // Optional encoding
  if (params.useSATConstraints) {
    encodeCross2Variables(model, graph, params.verbose);
  }
  if (params.useUNSATConstraints) {
    CHECK(params.useSATConstraints);
    encodeCross1Constraints(model, graph, params.verbose);
  }
  if (params.useIC) {
    CHECK(params.useSATConstraints);
    encodeICConstraints(model, graph, params.verbose, 0);
  }
  if (params.useNIC) {
    CHECK(params.useSATConstraints);
    encodeICConstraints(model, graph, params.verbose, 1);
  }

  // Symmetry
  encodeStackSymmetry(model, graph, params.verbose);
}

/// Extract a result from the SAT encoding for move-based planarity
void fillResultMove(
    const SATModel& model, 
    Solver& solver, 
    const InputGraph& graph,
    const Params& params,
    Result& result) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  const int numSegments = 2 * m;

  // Verify relative orders
  for (int i = n; i < numVertices; i++) {
    for (int j = n; j < numVertices; j++) {
      if (i == j)
        continue;
      if (canBeMerged(i, j, n, edges))
        continue;
      CHECK(model.value(solver, model.getRelVar(i, j, true)) != model.value(solver, model.getRelVar(j, i, true)));
    }
  }

  // Extract the order and crossings
  result.order = std::vector<std::vector<int>>(numVertices);
  result.isCrossed = std::vector<bool>(m, false);
  for (int i = 0; i < numVertices; i++) {
    std::vector<int> smaller;
    std::vector<int> equal;

    for (int j = 0; j < numVertices; j++) {
      if (i == j)
        continue;
      if (model.value(solver, model.getRelVar(i, j, true))) {
        CHECK(!model.value(solver, model.getRelVar(j, i, true)), "problem with vertices %d and %d", i, j);
        smaller.push_back(j);
      }
      if (model.value(solver, model.getRelVar(i, j, false)) && model.value(solver, model.getRelVar(j, i, false))) {
        equal.push_back(j);
      }
    }
    // LOG("vertex %d; |smaller| = %d; |equal| = %d (%s)", i, smaller.size(), equal.size(), to_string(equal).c_str());

    CHECK((i < n && equal.size() == 0) || (i >= n && equal.size() <= 1));
    CHECK(numVertices >= int(smaller.size() + 1));
    const int index = numVertices - int(smaller.size() + 1);
    CHECK(result.order[index].size() <= 1);
    result.order[index].push_back(i);

    if (equal.size() == 1 && i < equal[0]) {
      const int e1 = i - n;
      const int e2 = equal[0] - n;
      CHECK(e1 != e2);
      CHECK(!graph.adjacent(edges[e1], edges[e2]));
      result.crossings.push_back({e1, e2});

      CHECK(!result.isCrossed[e1] && !result.isCrossed[e2]);
      result.isCrossed[e1] = true;
      result.isCrossed[e2] = true;
    }
  }
  // Verify the order
  for (int i = 0; i < numVertices; i++) {
    if (result.order[i].empty())
      continue;
    for (int j = 0; j < numVertices; j++) {
      if (result.order[j].empty())
        continue;
      for (int x : result.order[i]) {
        for (int y : result.order[j]) {
          if (x == y)
            continue;
          CHECK(i != j || model.value(solver, model.getRelVar(x, y, false)));
          CHECK(i >= j || model.value(solver, model.getRelVar(x, y, true)));
          CHECK(i <= j || model.value(solver, model.getRelVar(x, y, false)));
        }
      }
    }
  }
  // Verify move variables
  for (int i = 0; i < numVertices; i++) {
    if (result.order[i].size() != 2)
      continue;
    const int v1 = result.order[i][0];
    const int v2 = result.order[i][1];
    for (int seg = 0; seg < numSegments; seg++) {
      // skip adjacent
      const auto [e_first, e_second] = graph.seg2edge(seg);
      if (e_first == v1 || e_second == v1 || e_first == v2 || e_second == v2)
        continue;
      CHECK(model.value(solver, model.getMoveVar(v1, seg, true)) == model.value(solver, model.getMoveVar(v2, seg, true)), 
        "move constraints not synced for merged vertices %d and %d", v1, v2);
    }
  }
}

/// Extract a result from the SAT encoding for stack-based planarity
void fillResultStack(
    const SATModel& model, 
    Solver& solver, 
    const InputGraph& graph,
    const Params& params,
    Result& result) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  const int numSegments = 2 * m;

  // Verify relative orders
  for (int i = n; i < numVertices; i++) {
    for (int j = n; j < numVertices; j++) {
      if (i == j)
        continue;
      if (canBeMerged(i, j, n, edges))
        continue;
      CHECK(model.value(solver, model.getRelVar(i, j, true)) != model.value(solver, model.getRelVar(j, i, true)));
    }
  }

  // Extract the order and crossings
  result.order = std::vector<std::vector<int>>(numVertices);
  std::vector<int> invOrder(numVertices, -1);
  result.isCrossed = std::vector<bool>(m, false);
  for (int i = 0; i < numVertices; i++) {
    std::vector<int> smaller;
    std::vector<int> equal;

    for (int j = 0; j < numVertices; j++) {
      if (i == j)
        continue;
      if (model.value(solver, model.getRelVar(i, j, true))) {
        CHECK(!model.value(solver, model.getRelVar(j, i, true)), "problem with vertices %d and %d", i, j);
        smaller.push_back(j);
      }
      if (model.value(solver, model.getRelVar(i, j, false)) && model.value(solver, model.getRelVar(j, i, false))) {
        equal.push_back(j);
      }
    }

    CHECK((i < n && equal.size() == 0) || (i >= n && equal.size() <= 1));
    CHECK(numVertices >= int(smaller.size() + 1));
    const int index = numVertices - int(smaller.size() + 1);
    CHECK(result.order[index].size() <= 1);
    result.order[index].push_back(i);
    invOrder[i] = index;

    if (equal.size() == 1 && i < equal[0]) {
      const int e1 = i - n;
      const int e2 = equal[0] - n;
      CHECK(e1 != e2);
      CHECK(!graph.adjacent(edges[e1], edges[e2]));
      CHECK(canBeMerged(i, equal[0], n, edges));
      result.crossings.push_back({e1, e2});

      CHECK(!result.isCrossed[e1] && !result.isCrossed[e2]);
      result.isCrossed[e1] = true;
      result.isCrossed[e2] = true;
      // verify corresponding cross-2 variable
      if (params.useSATConstraints) {
        CHECK(model.value(solver, model.getCross2Var(e1 + n, e2 + n, true)), "problem cross2Var for edges %d and %d", e1, e2);
      }
    }
  }
  // Verify cross-1 variables
  if (params.useUNSATConstraints) {
    for (int e = 0; e < m; e++) {
      if (result.isCrossed[e]) {
        CHECK(model.value(solver, model.getCross1Var(e + n, true)), "problem cross1Var for edge %d (%d, %d)", e, edges[e].first, edges[e].second);
      } else {
        CHECK(model.value(solver, model.getCross1Var(e + n, false)), "problem cross1Var for edge %d (%d, %d)", e, edges[e].first, edges[e].second);
      }
    }
  }
  // Verify the order
  for (int i = 0; i < numVertices; i++) {
    if (result.order[i].empty())
      continue;
    for (int j = 0; j < numVertices; j++) {
      if (result.order[j].empty())
        continue;
      for (int x : result.order[i]) {
        for (int y : result.order[j]) {
          if (x == y)
            continue;
          CHECK(i != j || model.value(solver, model.getRelVar(x, y, false)));
          CHECK(i >= j || model.value(solver, model.getRelVar(x, y, true)));
          CHECK(i <= j || model.value(solver, model.getRelVar(x, y, false)));
        }
      }
    }
  }

  // Fill in stacks / pages
  result.stack = std::vector<bool>(numSegments);
  for (int seg = 0; seg < numSegments; seg++) {
    const auto [e_org, e_div] = graph.seg2edge_v2(seg);
    CHECK(e_org < e_div && e_div >= n);
    result.stack[seg] = (invOrder[e_org] < invOrder[e_div]);
  }
  // Verify crossings on stacks
  std::vector<int> vIndex(numVertices, -1);
  for (int i = 0; i < (int)result.order.size(); i++) {
    if (result.order[i].empty())
      continue;
    CHECK(1 <= result.order[i].size() && result.order[i].size() <= 2, "%d", result.order[i].size());
    for (int v : result.order[i]) {
      CHECK(vIndex[v] == -1);
      vIndex[v] = i;
    }
  }

  for (int seg_e = 0; seg_e < numSegments; seg_e++) {
    for (int seg_f = seg_e + 1; seg_f < numSegments; seg_f++) {
      const auto [e_first, e_second] = graph.seg2edge(seg_e);
      const auto [f_first, f_second] = graph.seg2edge(seg_f);

      if (result.stack[seg_e] != result.stack[seg_f])
        continue;

      const int el = std::min(vIndex[e_first], vIndex[e_second]);
      const int er = std::max(vIndex[e_first], vIndex[e_second]);
      const int fl = std::min(vIndex[f_first], vIndex[f_second]);
      const int fr = std::max(vIndex[f_first], vIndex[f_second]);
      CHECK(el != -1 && el < er);
      CHECK(fl != -1 && fl < fr);

      CHECK(!cross(el, er, fl, fr), "edges (%d, %d) and (%d, %d) cross", e_first, e_second, f_first, f_second);
    }
  }
}

///
Result runSolver(const Params& params, const InputGraph& graph) {
  const int verbose = params.verbose;

  // Init the model
  SATModel model;
  if (params.useMovePlanarity)
    encodeMovePlanar(model, graph, params);
  else
    encodeStackPlanar(model, graph, params);

  // Create a solver
  Solver solver;
  solver.verbosity = verbose;
  if (params.timeout > 0) {
    solver.timeout_ms = 1000 * params.timeout;
  }

  // Init a model
  model.initVars(solver);

  // Symmetry-breaking
  if (params.applySatsuma) {
    LOG_IF(verbose, "applying Satsuma for %'d variables and %'d constraints", model.varCount(), model.clauseCount());
    model.applySatsuma(verbose, solver);
  }

  // Init a model
  model.initClauses(solver);

  if (params.applyBreakID) {
    LOG_IF(verbose, "applying BreakID for %'d variables and %'d constraints", model.varCount(), model.clauseCount());
    model.applyBreakID(verbose, solver);
  }

  if (params.modelFile != "") {
    LOG_IF(verbose, "encoded %'d variables and %'d constraints", model.varCount(), model.clauseCount());
    model.toDimacs(params.modelFile);
    LOG("SAT model in dimacs format saved to '%s'", params.modelFile.c_str());
    exit(0);
  }

  lbool ret;
  if (params.resultFile != "") {
    LOG("parsing SAT model with %'d variables and %'d clauses from '%s'", model.varCount(), model.clauseCount(), params.resultFile.c_str());
    auto externalResult = model.fromDimacs(params.resultFile);

    if (externalResult == "SATISFIABLE") {
      ret = l_True;
    } else if (externalResult == "UNSATISFIABLE") {
      ret = l_False;
    } else {
      ret = l_Undef;
    }
  } else {
    solver.simplify();
    LOG_IF(verbose, "solving SAT model with %'d variables and %'d clauses...", solver.nVars(), solver.nClauses());

    if (!solver.okay()) {
      ret = l_False;
    } else {
      ret = solver.solve();
    }
  }
  CHECK(ret == l_True || ret == l_False || ret == l_Undef, "an error within SAT solver");

  Result result;
  if (ret == l_True) {
    result.code = ResultCodeTy::SAT;
    if (params.useMovePlanarity)
      fillResultMove(model, solver, graph, params, result);
    else
      fillResultStack(model, solver, graph, params, result);
  } else if (ret == l_False) {
    result.code = ResultCodeTy::UNSAT;
  } else if (ret == l_Undef) {
    result.code = ResultCodeTy::TIMEOUT;
    LOG_IF(verbose, "  time limit exceeded");
  } else {
    result.code = ResultCodeTy::ERROR;
  }

  LOG_IF(verbose, "SAT solver completed with return code: %d (%s)", result.code, result.getCodeDesc().c_str());
  return result;
}
