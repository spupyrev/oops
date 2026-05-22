--- a/src/forbidden_crossings.h
+++ b/src/forbidden_crossings.h
@@ -92,6 +92,13 @@
     return clauses3;
   }

+  void clear() {
+    forbiddenCrossingPairs.clear();
+    forbiddenCrossingTriples.clear();
+    clauses2.clear();
+    clauses3.clear();
+  }
+
 private:
   // std::unordered_set<int64_t> forbiddenCrossingPairs;
   // std::unordered_set<int64_t> forbiddenCrossingTriples;
--- a/src/glucose/SolverSimp21.h
+++ b/src/glucose/SolverSimp21.h
@@ -369,6 +369,9 @@

 public:
   int level(Var x) const;
+  // Public wrapper exposing varBumpActivity for callers that want to
+  // prioritise domain-specific decision variables.
+  void publicBumpActivity(Var v, double mult) { varBumpActivity(v, mult); }

 protected:
   template <class V> int computeLBD(const V &c) {
--- a/src/main.cpp
+++ b/src/main.cpp
@@ -68,6 +68,8 @@
   options.addAllowedOption("-partial-constraints", "", "Add partial constraints: num_pairs");
   options.addAllowedOption("-sep-cycles", "", "Add constraints based on separating cycles: num_pairs");
   options.addAllowedOption("-up-sepcycles", "false", "Use sep-cycle user propagation");
+  options.addAllowedOption("-cross-priority", "false", "Bias CDCL to cross1/cross2 first; cross1=true preferred; phase_saving=0");
+  options.addAllowedOption("-witness-rule", "false", "Apply witness-rule filter inside initCrossablePairs to prune crossing pairs");

   // External SAT solver
   options.addAllowedOption("-dimacs", "", "Output dimacs file");
@@ -427,6 +429,8 @@
   params.custom = options.getStr("-custom");
   params.ignoreTransitiveRels = options.hasCustomOption("no-transitive");
   params.useSepCycleUP = options.getBool("-up-sepcycles");
+  params.crossPriority = options.getBool("-cross-priority");
+  params.useWitnessRule = options.getBool("-witness-rule");

   params.useSATConstraints = options.getBool("-sat");
   const int unsatLevel = options.getInt("-unsat");
--- a/src/one_planar.cpp
+++ b/src/one_planar.cpp
@@ -77,6 +77,11 @@
 // Pairs of edges that could be crossed; equivalently, pairs of division vertices that can be merged
 std::vector<std::vector<bool>> crossablePairs;

+// Cache of (cross2, cross2) pair-pair incompatibilities populated by
+// initCrossablePairs's witness-rule filter and consumed later by
+// encodeSepCyclesConstraints. Lifecycle parallels crossablePairs.
+ForbiddenTuples sepCycleCache;
+
 /// Find all pairs of edges that can be crossed
 void initCrossablePairs(const Params& params, const InputGraph& graph) {
   const int n = graph.n;
@@ -218,6 +223,20 @@
     }
   }

+  // Witness-rule filter: forbid pair (A, B) if there exists an edge e (not
+  // adjacent to A or B) such that EVERY admissible crossing partner X of e
+  // has the pair-pair (A, B), (e, X) flagged as incompatible by the static
+  // sep-cycle pair-pair pass. Such (A, B) is implied UNSAT by the existing
+  // 2-clauses combined with the at-least-one-cross clause for e.
+  int numWitnessSkipped = 0;
+  if (params.useWitnessRule) {
+    numWitnessSkipped = applyWitnessRuleFilter(graph, params,
+                                               crossablePairs, sepCycleCache);
+  } else {
+    sepCycleCache.clear();
+  }
+  (void)numWitnessSkipped;
+
   // Count pairs
   int mergablePairs = 0;
   int possiblePairs = 0;
@@ -1381,6 +1400,11 @@

   // Init the model
   SATModel model;
+  // Seed the model's forbiddenTuples cache with pair-pair incompatibilities
+  // discovered by the witness-rule filter inside initCrossablePairs.
+  for (const auto& cp : sepCycleCache.pairs()) {
+    model.getForbiddenTuples().insert(cp);
+  }
   if (params.solverType == SolverType::MOVE)
     encodeMovePlanar(model, graph, params);
   else
@@ -1413,6 +1437,33 @@
   if (params.useSepCycleUP && params.modelFile == "" && params.resultFile == "") {
     solver.setUserPropagator(createSepCyclesDynamic(model, graph, params));
   }
+
+  // Cross-priority: bump activity on cross1/cross2 vars (CDCL branches on
+  // them first), set cross1 polarity to TRUE-first (try `edge crosses`),
+  // and disable phase saving so the polarity bias persists across restarts.
+  if (params.crossPriority) {
+    solver.phase_saving = 0;
+    const double crossBoost = 1e6;
+    int bumpedCross1 = 0, bumpedCross2 = 0;
+    for (size_t e = 0; e < graph.edges.size(); e++) {
+      const int divE = (int)e + graph.n;
+      if (model.getCross1Vars().count(divE)) {
+        const Simp21::Lit lit = model.getSolverLit(model.getCross1Var(divE, true));
+        const Simp21::Var v = Simp21::var(lit);
+        solver.publicBumpActivity(v, crossBoost);
+        // Glucose: polarity[v]=false means try `v=true` first.
+        solver.setPolarity(v, false);
+        bumpedCross1++;
+      }
+    }
+    for (const auto& kv : model.getCross2Vars()) {
+      const Simp21::Lit lit = model.getSolverLit(model.getCross2Var(kv.first.first, kv.first.second, true));
+      solver.publicBumpActivity(Simp21::var(lit), crossBoost);
+      bumpedCross2++;
+    }
+    LOG_IF(verbose, "cross-priority: bumped %d cross1 + %d cross2 vars; phase_saving=0",
+           bumpedCross1, bumpedCross2);
+  }
   const auto endTimeEncoding = chrono::steady_clock::now();
   LOG_IF(verbose >= 1, "SAT encoding took %s", ms_to_str(startTimeEncoding, endTimeEncoding).c_str());

   --- a/src/one_planar.h
+++ b/src/one_planar.h
@@ -220,6 +220,14 @@
   std::string custom = "";
   bool ignoreTransitiveRels = false;
   bool useSepCycleUP = false;
+  // Bias CDCL toward branching on cross1/cross2 vars first, with cross1=true
+  // preferred. Disables phase saving so the polarity bias persists across
+  // restarts. Strong on g.31.9-class instances; mixed on rome-style ones.
+  bool crossPriority = false;
+  // Apply the witness-rule filter inside initCrossablePairs to forbid
+  // additional (cross2, cross2) pairs whose constraints alone make the
+  // pair unreachable. Strong on rome-style UNSAT instances.
+  bool useWitnessRule = false;

   std::string to_string() const {
     std::ostringstream ss;
@@ -499,6 +507,9 @@
     return MVar(cross1Vars.at(edge), positive);
   }

+  const std::unordered_map<int, int>& getCross1Vars() const { return cross1Vars; }
+  const std::unordered_map<std::pair<int, int>, int, pair_hash>& getCross2Vars() const { return cross2Vars; }
+
   void addCross1Var(int edge) {
     CHECK(!cross1Vars.count(edge));
     cross1Vars[edge] = addVar();
@@ -855,3 +866,11 @@
     const SATModel& model,
     const InputGraph& graph,
     const Params& params);
+
+// Witness-rule filter: forbid (A, B) pairs where some non-incident edge e
+// has all admissible partners X already known to be incompatible with (A, B).
+// Mutates crossablePairs in place; populates `cache` with discovered
+// (cross2, cross2) pair-pair incompatibilities.
+int applyWitnessRuleFilter(const InputGraph& graph, const Params& params,
+                           std::vector<std::vector<bool>>& crossablePairs,
+                           ForbiddenTuples& cache);
--- a/src/sepcycles_static.cpp
+++ b/src/sepcycles_static.cpp
@@ -7,8 +7,8 @@
   static constexpr size_t MAX_CROSSINGS_FOR_3GT = 1280;
   static constexpr size_t MAX_CROSSINGS_FOR_2GT = 3840;

   -  SepCyclesStatic(SATModel& model, const InputGraph& graph, const Params& params, ForbiddenTuples& tuples)
-    : model(model), graph(graph), params(params), n(graph.n), m((int)graph.edges.size()),
+  SepCyclesStatic(const InputGraph& graph, const Params& params, ForbiddenTuples& tuples)
+    : graph(graph), params(params), n(graph.n), m((int)graph.edges.size()),
       verbose(params.verbose), forbiddenTuples(tuples)
     {}

@@ -42,76 +42,101 @@
     return true;
   }

-  /// Find pairs of crossings that cannot happen in a 1-planar drawing
-  size_t build2Clauses() {
+  /// Find pairs of crossings that cannot happen in a 1-planar drawing.
+  /// Records the (cross2, cross2) pair-pair incompatibilities into
+  /// forbiddenTuples. Does not emit anything to a SAT model.
+  size_t findForbiddenCrossingPairs() {
     if (possibleCrossings.size() > MAX_CROSSINGS_FOR_2GT) {
       LOG_IF(verbose, "  skipped building sep-cycles 2-clauses due to too many possible_crossings");
       return 0;
     }

     size_t numClauses2 = 0;
-    size_t numEqFlowClauses2 = 0;
-    size_t numEqFlowClauses3 = 0;

-    // the main crossing
     for (size_t cross0 = 0; cross0 < possibleCrossings.size(); cross0++) {
       const auto [e1_cross0, e2_cross0] = possibleCrossings[cross0];
       const auto [u, v] = graph.edges[e1_cross0];
       const auto [x, y] = graph.edges[e2_cross0];
       CHECK(all_unique({x, y, u, v}));

-      // circular order: u--x--v--y
       markCrossing(u, v, x, y);

-      std::vector<std::pair<int, int>> takenCrossings0 = {{e1_cross0, e2_cross0}};
-      numEqFlowClauses2 += findEqualFlow(x, y, u, v, takenCrossings0);
-      numEqFlowClauses2 += findEqualFlow(u, v, x, y, takenCrossings0);
-
-      // secondary crossings
       for (size_t cross1 = 0; cross1 < possibleCrossings.size(); cross1++) {
         const auto [e1_cross1, e2_cross1] = possibleCrossings[cross1];
         if (cross0 == cross1)
           continue;
-        // circular order: u2--x2--v2--y2
         const auto [u2, v2] = graph.edges[e1_cross1];
         const auto [x2, y2] = graph.edges[e2_cross1];

         -        // check that the selected crossings are compatibe
         if (!canMarkCrossing(u2, v2, x2, y2))
           continue;
         CHECK(all_unique({e1_cross0, e2_cross0, e1_cross1, e2_cross1}));

-        // stop early, if the crossing pair is already processed
         const CrossingPair crossPair(m, e1_cross0, e2_cross0, e1_cross1, e2_cross1);
         if (forbiddenTuples.contains(crossPair))
           continue;

         markCrossing(u2, v2, x2, y2);

-        // need to search cycles both via (x, y) and via (u, v)
         std::vector<std::pair<int, int>> takenCrossings = {{e1_cross0, e2_cross0}, {e1_cross1, e2_cross1}};
         if (findSepCycle(x, y, u, v, takenCrossings) || findSepCycle(u, v, x, y, takenCrossings)) {
-          // save the 2-clause
           numClauses2++;
           forbiddenTuples.insert(crossPair);
         }

-        if (possibleCrossings.size() < MAX_CROSSINGS_FOR_3EQ) {
-          numEqFlowClauses3 += findEqualFlow(x, y, u, v, takenCrossings);
-          numEqFlowClauses3 += findEqualFlow(u, v, x, y, takenCrossings);
-        }
-
         unmarkCrossing(u2, v2, x2, y2);
       }

       unmarkCrossing(u, v, x, y);
     }

+    return numClauses2;
+  }
+
+  /// Emit equal-flow constraints to the model.
+  /// 2-clauses: cross2[A,B] => cross1[e] for free cycle edges e.
+  /// 3-clauses: cross2[A,B] AND cross2[A',B'] => cross1[e] (gated on
+  /// MAX_CROSSINGS_FOR_3EQ).
+  void encodeEqualFlowConstraints(SATModel& model) {
+    size_t numEqFlowClauses2 = 0;
+    size_t numEqFlowClauses3 = 0;
+
+    for (size_t cross0 = 0; cross0 < possibleCrossings.size(); cross0++) {
+      const auto [e1_cross0, e2_cross0] = possibleCrossings[cross0];
+      const auto [u, v] = graph.edges[e1_cross0];
+      const auto [x, y] = graph.edges[e2_cross0];
+      CHECK(all_unique({x, y, u, v}));
+
+      markCrossing(u, v, x, y);
+
+      std::vector<std::pair<int, int>> takenCrossings0 = {{e1_cross0, e2_cross0}};
+      numEqFlowClauses2 += findEqualFlow(x, y, u, v, takenCrossings0, &model);
+      numEqFlowClauses2 += findEqualFlow(u, v, x, y, takenCrossings0, &model);

+      if (possibleCrossings.size() < MAX_CROSSINGS_FOR_3EQ) {
+        for (size_t cross1 = 0; cross1 < possibleCrossings.size(); cross1++) {
+          const auto [e1_cross1, e2_cross1] = possibleCrossings[cross1];
+          if (cross0 == cross1)
+            continue;
+          const auto [u2, v2] = graph.edges[e1_cross1];
+          const auto [x2, y2] = graph.edges[e2_cross1];
+
+          if (!canMarkCrossing(u2, v2, x2, y2))
+            continue;
+
+          markCrossing(u2, v2, x2, y2);
+          std::vector<std::pair<int, int>> takenCrossings = {{e1_cross0, e2_cross0}, {e1_cross1, e2_cross1}};
+          numEqFlowClauses3 += findEqualFlow(x, y, u, v, takenCrossings, &model);
+          numEqFlowClauses3 += findEqualFlow(u, v, x, y, takenCrossings, &model);
+          unmarkCrossing(u2, v2, x2, y2);
+        }
+      }
+
+      unmarkCrossing(u, v, x, y);
+    }
+
     LOG_IF(verbose, "  added %'9d equal-flow 2-clauses", numEqFlowClauses2);
     LOG_IF(verbose, "  added %'9d equal-flow 3-clauses", numEqFlowClauses3);
-
-    LOG_IF(verbose && clauseIndex, "num_extra_clauses: %d", clauseIndex);
-    return numClauses2;
   }

   /// Find triples of crossings that cannot happen in a 1-planar drawing
@@ -473,7 +498,16 @@
     return true;
   }

-  int findEqualFlow(int x, int y, int u, int v, const std::vector<std::pair<int, int>> &takenCrossings) {
+  // Callback signature: (cr1_e1, cr1_e2, edge_v1, edge_v2). Optional. When
+  // non-null, called for each free edge (v1,v2) on a found equal-flow cycle.
+  // findEqualFlow may run with model==nullptr if only the callback is needed
+  // (e.g., for the witness-rule's equal-flow witness recording pass).
+  using EqualFlowCallback = std::function<void(int, int, int, int)>;
+

+  int findEqualFlow(int x, int y, int u, int v,
+                    const std::vector<std::pair<int, int>> &takenCrossings,
+                    SATModel* model,
+                    const EqualFlowCallback* cb = nullptr) {
     CHECK(takenCrossings.size() == 1 || takenCrossings.size() == 2);
     int numAddedClauses = 0;
     const int degreeBound = std::min(graph.degree(u) - 1, graph.degree(v) - 1);
@@ -499,11 +533,17 @@
             continue;
           // the edge must cross smth
           numAddedClauses++;
-          if (takenCrossings.size() == 1) {
-            addEqFlowClause(takenCrossings[0], v1, v2);
-          } else {
-            CHECK(takenCrossings.size() == 2);
-            addEqFlowClause(takenCrossings[0], takenCrossings[1], v1, v2);
+          if (model) {
+            if (takenCrossings.size() == 1) {
+              addEqFlowClause(*model, takenCrossings[0], v1, v2);
+            } else {
+              CHECK(takenCrossings.size() == 2);
+              addEqFlowClause(*model, takenCrossings[0], takenCrossings[1], v1, v2);
+            }
+          }
+          if (cb) {
+            const auto& cr1 = takenCrossings[0];
+            (*cb)(cr1.first, cr1.second, v1, v2);
           }
         }
       }
@@ -586,28 +626,27 @@
     return crossCycleEdges;
   }

-  void addEqFlowClause(std::pair<int, int> cr1, int u, int v) {
+  void addEqFlowClause(SATModel& model, std::pair<int, int> cr1, int u, int v) {
     const auto [e1, e2] = cr1;
     const int divUV = graph.findDivIndex(u, v);
     model.addClause({
-        model.getCross2Var(e1 + n, e2 + n, false),
+        model.getCross2Var(e1 + n, e2 + n, false),
         model.getCross1Var(divUV, true)
     });
   }
-  void addEqFlowClause(std::pair<int, int> cr1, std::pair<int, int> cr2, int u, int v) {
+  void addEqFlowClause(SATModel& model, std::pair<int, int> cr1, std::pair<int, int> cr2, int u, int v) {
     const auto [e1, e2] = cr1;
     const auto [e3, e4] = cr2;
     const int divUV = graph.findDivIndex(u, v);
     model.addClause({
-        model.getCross2Var(e1 + n, e2 + n, false),
-        model.getCross2Var(e3 + n, e4 + n, false),
+        model.getCross2Var(e1 + n, e2 + n, false),
+        model.getCross2Var(e3 + n, e4 + n, false),
         model.getCross1Var(divUV, true)
     });
   }

 private:
-  SATModel& model;
   const InputGraph& graph;
   const Params& params;
   const int n;
@@ -629,20 +668,119 @@
   int clauseIndex = 0;
 };

+/// Witness-rule filter: forbid pair (A, B) if there exists an edge e (not
+/// adjacent to either edge of (A, B)) such that EVERY admissible crossing
+/// partner X of e has the pair-pair (A, B), (e, X) flagged as incompatible.
+int applyWitnessRuleFilter(const InputGraph& graph, const Params& params,
+                           std::vector<std::vector<bool>>& crossablePairs,
+                           ForbiddenTuples& cache) {
+  cache.clear();
+  if (graph.isDirected()) return 0;
+
+  const int n = graph.n;
+  const int m = (int)graph.edges.size();
+
+  SepCyclesStatic sct(graph, params, cache);
+  if (!sct.init(2)) return 0;
+
+  // Phase 1a: populate cache with K4 (cross2, cross2) pair-pair incompatibilities.
+  // Mirrors encodeK4Constraints — pair-pairs (divUV, divXY) and (divVB, divUA)
+  // where edge-pair shares K4-structured incident edges.
+  {
+    const auto& adjLocal = graph.adj;
+    for (int divUV = n; divUV < n + m; divUV++) {
+      for (int divXY = n; divXY < n + m; divXY++) {
+        if (divUV == divXY) continue;
+        if (!crossablePairs[divUV][divXY]) continue;
+        const int u = graph.edges[divUV - n].first;
+        const int v = graph.edges[divUV - n].second;
+        for (const int a : adjLocal[u]) {
+          if (a == v) continue;
+          const int divUA = graph.findDivIndex(u, a);
+          for (const int b : adjLocal[v]) {
+            if (b == u) continue;
+            const int divVB = graph.findDivIndex(v, b);
+            if (!crossablePairs[divUA][divVB]) continue;
+            if (divUV == divVB || divUV == divUA) continue;
+            const CrossingPair cp(m, divUV - n, divXY - n, divVB - n, divUA - n);
+            if (cache.contains(cp)) continue;
+            cache.insert(cp);
+          }
+        }
+      }
+    }
+  }
+
+  // Phase 1b: populate cache with sep-cycle (cross2, cross2) pair-pair
+  // incompatibilities using the current crossablePairs state.
+  sct.findForbiddenCrossingPairs();
+
+  // Phase 2: apply the witness rule. Mutate crossablePairs in place so later
+  // iterations see fewer admissible partners.
+  int numForbidden = 0;
+  for (int e1A = 0; e1A < m; e1A++) {
+    for (int e2A = e1A + 1; e2A < m; e2A++) {
+      if (!crossablePairs[e1A + n][e2A + n]) continue;
+
+      bool forbidden = false;
+      for (int eIdx = 0; eIdx < m && !forbidden; eIdx++) {
+        if (eIdx == e1A || eIdx == e2A) continue;
+        if (graph.adjacentEE(graph.edges[eIdx], graph.edges[e1A])) continue;
+        if (graph.adjacentEE(graph.edges[eIdx], graph.edges[e2A])) continue;
+
+        int admissibleCount = 0;
+        int blockedCount = 0;
+        for (int X = 0; X < m; X++) {
+          if (X == eIdx || X == e1A || X == e2A) continue;
+          if (!crossablePairs[eIdx + n][X + n]) continue;
+          admissibleCount++;
+          CrossingPair pp(m, e1A, e2A, eIdx, X);
+          if (cache.contains(pp)) blockedCount++;
+        }
+        if (admissibleCount > 0 && admissibleCount == blockedCount) {
+          forbidden = true;
+        }
+      }
+      if (forbidden) {
+        crossablePairs[e1A + n][e2A + n] = false;
+        crossablePairs[e2A + n][e1A + n] = false;
+        numForbidden++;
+      }
+    }
+  }
+

+  // Phase 3: prune cache so the encoder doesn't try to emit clauses on
+  // pruned pairs (their cross2 vars won't exist).
+  if (numForbidden > 0) {
+    ForbiddenTuples kept;
+    for (const auto& cp : cache.pairs()) {
+      if (!crossablePairs[cp.e1 + n][cp.e2 + n]) continue;
+      if (!crossablePairs[cp.p1 + n][cp.p2 + n]) continue;
+      kept.insert(cp);
+    }
+    cache = std::move(kept);
+  }
+
+  return numForbidden;
+}
+
 /// Add constraints based on separating cycles
 void encodeSepCyclesConstraints(SATModel& model, const InputGraph& graph, const Params& params) {
   const int verbose = params.verbose;
   const int numPairs = to_int(params.sepCycleConstraints);
   CHECK(2 <= numPairs && numPairs <= 3, "incorrect value of sepCycleConstraints");

-  SepCyclesStatic sct(model, graph, params, model.getForbiddenTuples());
+  SepCyclesStatic sct(graph, params, model.getForbiddenTuples());
   if (sct.init(numPairs)) {
-    // forbidden pairs (2-clauses)
-    const size_t numClauses2 = sct.build2Clauses();
+    // forbidden pairs (2-clauses) — may already be partially populated by the
+    // witness-rule filter inside initCrossablePairs; missing entries are added.
+    const size_t numClauses2 = sct.findForbiddenCrossingPairs();
     LOG_IF(verbose, "  found %'9d 2-clauses", numClauses2);

+    // equal-flow constraints
+    sct.encodeEqualFlowConstraints(model);
+
     if (numPairs == 3 || params.unsatLevel >= 2) {
-      // forbidden triples (3-clauses)
       const size_t numClauses3 = sct.build3Clauses();
       LOG_IF(verbose, "  found %'9d 3-clauses", numClauses3);
     }
