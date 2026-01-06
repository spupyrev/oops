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

/// Return true iff the two (division) vertices can be merged
bool canBeMerged(int u, int v, const int n, const std::vector<EdgeTy>& edges) {
  CHECK(u != v);
  CHECK(u >= n && v >= n);
  CHECK((int)crossablePairs.size() == n + (int)edges.size());
  return crossablePairs[u][v];
}

bool canBeMerged(int u, int v, const InputGraph& graph) {
  return canBeMerged(u, v, graph.n, graph.edges);
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

  // Ensure associativity for mergable vertices
  model.reserveClauses(numVertices * numVertices * numVertices + numVertices * numVertices);
  for (int i = n; i < numVertices; i++) {
    for (int j = n; j < numVertices; j++) {
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

  LOG_IF(verbose, "cross-2 constraints:");

  if (!graph.isDirected()) {
    encodeK4Constraints(model, graph, verbose);
    encodeICConstraints(model, graph, verbose, 2);
  }
}

/// Encode edge-crossing variables and constraints
void encodeCross1Constraints(SATModel& model, const InputGraph& graph, const Params& params) {
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

  LOG_IF(params.verbose, "cross-1 constraints: %d planar-K4; %d degree-2-lex", numPlanarK4, numDegree2Lex);

  if (!graph.isDirected()) {
    encodeCoverConstraints(model, graph, params);
  }
}

/// Encode cover variables and constraints
void encodeCoverConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;

  // Variables
  model.reserveCoverVars(m, numVertices);
  for (int div = n; div < numVertices; div++) {
    const int e = div - n;
    const auto [s, t] = edges[e];

    // edge e covers its division vertex, div
    model.addCoverVar(e, div);

    // Cover is true => vertex is in between
    model.addClause({
        model.getCoverVar(e, div, false), 
        model.getRelVar(div, s, true),
        model.getRelVar(div, t, true)
    });
    model.addClause({
        model.getCoverVar(e, div, false), 
        model.getRelVar(s, div, true),
        model.getRelVar(t, div, true)
    });

    // Between => cover
    model.addClause({
      model.getRelVar(div, s, true),
      model.getRelVar(t, div, true),
      model.getCoverVar(e, div, true)
    });
    model.addClause({
      model.getRelVar(div, t, true),
      model.getRelVar(s, div, true),
      model.getCoverVar(e, div, true)
    });
  }

  // Constraints
  // An edge doesn't cross => the division is in between
  int numMid1 = 0;
  for (int d = n; d < numVertices; d++) {
    const int e = d - n;
    model.addClause({
      model.getCross1Var(d, true), 
      model.getCoverVar(e, d, true)
    });
    numMid1++;
  }
  // Every crossed division vertex is in between one of the two edges
  int numMid2 = 0;
  for (int d1 = n; d1 < numVertices; d1++) {
    for (int d2 = d1 + 1; d2 < numVertices; d2++) {
      if (!canBeMerged(d1, d2, n, edges)) 
        continue;

      const int e1 = d1 - n;
      const int e2 = d2 - n;

      model.addClause({
        model.getCross2Var(d1, d2, false), 
        model.getCoverVar(e1, d1, true),
        model.getCoverVar(e2, d2, true)
      });
      numMid2++;
    }
  }

  // Two edges cross => edge intervals overlap
  int numStrict1 = 0;
  if (params.strict) {
    for (int div = n; div < numVertices; div++) {
      const int e = div - n;
      const auto [s, t] = edges[e];
      // edge e is directed s < t
      model.addDirVar(e);
      model.addClause({ model.getRelVar(s, t, true),  model.getDirVar(e, false) });
      model.addClause({ model.getRelVar(s, t, false), model.getDirVar(e, true) });
      model.addClause({ model.getRelVar(t, s, true),  model.getDirVar(e, true) });
      model.addClause({ model.getRelVar(t, s, false), model.getDirVar(e, false) });
    }

    // forbid relative order a < b < c < d  with d1+d2 crossed
    // auto addStrictClause = [&](int a, int b, int c, int d, int d1, int d2) {
    //   model.addClause({
    //     model.getCross2Var(d1, d2, false), 
    //     model.getRelVar(a, b, false),
    //     model.getRelVar(b, c, false),
    //     model.getRelVar(c, d, false)
    //   });
    // };

    auto addStrictClause2 = [&](int e1, bool guard1_is_dir,
                                int e2, bool guard2_is_dir,
                                int sepA, int sepB,
                                int d1, int d2) {
      model.addClause({
        model.getCross2Var(d1, d2, false),
        model.getDirVar(e1, !guard1_is_dir),
        model.getDirVar(e2, !guard2_is_dir),
        model.getRelVar(sepA, sepB, false)
      });
    };

    for (int d1 = n; d1 < numVertices; d1++) {
      for (int d2 = d1 + 1; d2 < numVertices; d2++) {
        if (!canBeMerged(d1, d2, n, edges)) 
          continue;

        const int e1 = d1 - n;
        const int e2 = d2 - n;
        const auto [u1, v1] = edges[e1];
        const auto [u2, v2] = edges[e2];

        // addStrictClause(u1, v1, u2, v2,  d1, d2);
        // addStrictClause(u1, v1, v2, u2,  d1, d2);
        // addStrictClause(v1, u1, u2, v2,  d1, d2);
        // addStrictClause(v1, u1, v2, u2,  d1, d2);
        // addStrictClause(u2, v2, u1, v1,  d1, d2);
        // addStrictClause(u2, v2, v1, u1,  d1, d2);
        // addStrictClause(v2, u2, u1, v1,  d1, d2);
        // addStrictClause(v2, u2, v1, u1,  d1, d2);

        // e1 before e2: forbid right(e1) < left(e2)
        addStrictClause2(e1, true,  e2, true,   v1, u2,  d1, d2);
        addStrictClause2(e1, true,  e2, false,  v1, v2,  d1, d2);
        addStrictClause2(e1, false, e2, true,   u1, u2,  d1, d2);
        addStrictClause2(e1, false, e2, false,  u1, v2,  d1, d2);
        // e2 before e1: forbid right(e2) < left(e1)
        addStrictClause2(e2, true,  e1, true,   v2, u1,  d1, d2);
        addStrictClause2(e2, true,  e1, false,  v2, v1,  d1, d2);
        addStrictClause2(e2, false, e1, true,   u2, u1,  d1, d2);
        addStrictClause2(e2, false, e1, false,  u2, v1,  d1, d2);

        numStrict1++;
      }
    }
  }

  LOG_IF(params.verbose, "cover constraints: %d mid-1; %d mid-2; %d strict1", numMid1, numMid2, numStrict1);
}

/// C=0 <=> IC; C=1 <=> NIC; C=2 <=> 1-planar
void encodeICConstraints(SATModel& model, const InputGraph& graph, const int verbose, const int C) {
  CHECK(0 <= C && C <= 3);

  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  ForbiddenTuples& forbiddenTuples = model.getForbiddenTuples();

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
            const CrossingPair crossPair(m, i1 - n, j1 - n, i2 - n, j2 - n);
            if (forbiddenTuples.contains(crossPair))
              continue;

            forbiddenTuples.insert(crossPair);
            numIC++;
          }
        }
      }
    }
  }

  LOG_IF(verbose, "  found %'9d 2-clauses (IC-%d)", numIC, C);
}

/// A pair cross => K4-edges do not cross
void encodeK4Constraints(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const auto& edges = graph.edges;
  const int m = (int)edges.size();
  const int numVertices = n + m;
  const auto& adj = graph.adj;
  ForbiddenTuples& forbiddenTuples = model.getForbiddenTuples();

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

      for (const int a: adj[u]) {
        if (a == v) 
          continue;
        const int divUA = graph.findDivIndex(u, a);
        for (const int b: adj[v]) {
          if (b == u) 
            continue;
          const int divVB = graph.findDivIndex(v, b);
          if (!canBeMerged(divUA, divVB, n, edges))
            continue;
          CHECK(divUV != divVB && divUV != divUA);

          const CrossingPair crossPair(m, divUV - n, divXY - n, divVB - n, divUA - n);
          if (forbiddenTuples.contains(crossPair))
            continue;

          forbiddenTuples.insert(crossPair);
          numConstraints++;
        }
      }
    }
  }

  LOG_IF(verbose, "  found %'9d 2-clauses (K4)", numConstraints);
}

/// Add identified 2- and 3-clauses to the model 
void encodeForbiddenTuples(SATModel& model, const InputGraph& graph, const int verbose) {
  const int n = graph.n;
  const ForbiddenTuples& forbiddenTuples = model.getForbiddenTuples();

  // 2-clauses
  for (const auto& clause : forbiddenTuples.pairs()) {
    const int divE1 = clause.e1 + n;
    const int divE2 = clause.e2 + n;
    CHECK(canBeMerged(divE1, divE2, graph));
    const int divP1 = clause.p1 + n;
    const int divP2 = clause.p2 + n;
    CHECK(canBeMerged(divP1, divP2, graph));

    model.addClause({
        model.getCross2Var(divE1, divE2, false),
        model.getCross2Var(divP1, divP2, false)
    });
  }
  LOG_IF(verbose, "encoded %'9d forbidden 2-clauses", forbiddenTuples.pairs().size());

  // 3-clauses
  for (const auto& clause : forbiddenTuples.triples()) {
    const int divA1 = clause.a1 + n;
    const int divA2 = clause.a2 + n;
    CHECK(canBeMerged(divA1, divA2, graph));
    const int divB1 = clause.b1 + n;
    const int divB2 = clause.b2 + n;
    CHECK(canBeMerged(divB1, divB2, graph));
    const int divC1 = clause.c1 + n;
    const int divC2 = clause.c2 + n;
    CHECK(canBeMerged(divC1, divC2, graph));

    model.addClause({
        model.getCross2Var(divA1, divA2, false),
        model.getCross2Var(divB1, divB2, false),
        model.getCross2Var(divC1, divC2, false),
    });
  }
  LOG_IF(verbose, "encoded %'9d forbidden 3-clauses", forbiddenTuples.triples().size());
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
  CHECK(params.useSATConstraints, "move-planarity should be used with -sat=1");

  // Main encoding
  encodeRelativeVariables(model, graph, params.verbose);
  encodeCross2Variables(model, graph, params.verbose);
  encodeMoveVariables(model, graph, params.verbose);
  encodeMoveConstraints(model, graph, params.verbose);

  // Optional encoding
  if (params.useUNSATConstraints) {
    CHECK(params.useSATConstraints);
    encodeCross1Constraints(model, graph, params);
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
  // auto addOrder = [&](const std::vector<int>& order) {
  //   for (size_t i = 0; i < order.size(); i++) {
  //     for (size_t j = i + 1; j < order.size(); j++) {
  //       model.addClause(MClause(model.getRelVar(order[i], order[j], true)));
  //     }
  //   }
  // };
  // auto addGroup = [&](const std::vector<int>& group) {
  //   for (int x = 0; x < n; x++) {
  //     if (contains(group, x)) 
  //       continue;

  //     for (size_t i = 0; i < group.size(); i++) {
  //       for (size_t j = 0; j < group.size(); j++) {
  //         if (i == j)
  //           continue;
  //         // forbid group[i] < x < group[j]
  //         model.addClause(MClause({model.getRelVar(group[i], x, true), model.getRelVar(x, group[j], true)}));
  //       }
  //     }
  //   }
  // };
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
    encodeCross1Constraints(model, graph, params);
  }
  if (params.useIC) {
    CHECK(params.useSATConstraints);
    encodeICConstraints(model, graph, params.verbose, 0);
  }
  if (params.useNIC) {
    CHECK(params.useSATConstraints);
    encodeICConstraints(model, graph, params.verbose, 1);
  }

  // Misc constraints containing 2- and 3-clauses
  if (params.swapConstraints != "") {
    CHECK(params.useSATConstraints);
    encodeSwapConstraints(model, graph, params);
  }
  if (params.partialConstraints != "") {
    CHECK(params.useSATConstraints);
    encodePartialConstraints(model, graph, params);
  }
  if (params.sepCycleConstraints != "") {
    CHECK(params.useSATConstraints);
    encodeSepCyclesConstraints(model, graph, params);
  }

  encodeForbiddenTuples(model, graph, params.verbose);

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

/// Run a SAT solver for a given graph
Result runSolver(const Params& params, const InputGraph& graph) {
  CHECK(params.solverType == SolverType::STACK || params.solverType == SolverType::MOVE);
  const int verbose = params.verbose;

  // Init the model
  SATModel model;
  if (params.solverType == SolverType::MOVE)
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
    if (params.solverType == SolverType::MOVE)
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
