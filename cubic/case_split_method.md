# An automorphism-driven case-split for non-1-planarity of symmetric graphs

## 1. Setting

Let $G = (V, E)$ be a finite simple graph. A drawing of $G$ is **1-planar**
if every edge participates in at most one crossing; the **NIC-1-planar**
variant additionally requires that no two pairs of crossing edges share a
vertex. The decision problem is whether $G$ admits a (NIC-)1-planar drawing.

The method targets highly symmetric graphs on which plain CDCL does not
converge because the search re-explores symmetric branches; it uses
$\mathrm{Aut}(G)$ to remove that redundancy. The audit in Section 5 and the
inactivity claims are worked out on two exemplars — the **Coxeter graph**
$C$ (28 vertices, 42 edges, 3-regular, $|\mathrm{Aut}(C)|=336$, edge- and
vertex-transitive, girth 7) and the **Tutte-Coxeter graph** $T$ (30
vertices, 45 edges, 3-regular, $|\mathrm{Aut}(T)|=1440$, edge- and
vertex-transitive, girth 8) — because their girth-$\ge 7$, 3-regular
structure makes the audit clean. Section 8 reports the instances actually
closed with the method (cubic girth-7 benchmarks).

We use the SAT encoder `oops`. Section 5 audits the encoder configuration
required for the case-split to be sound, and Section 6 describes the
`-fix-cross1` option the driver relies on.

## 2. Notation

For $\sigma \in \mathrm{Aut}(G)$ and $e = \{u,v\} \in E$, write
$\sigma(e) := \{\sigma(u), \sigma(v)\}$, and
$\mathrm{Stab}(S) := \{\sigma \in \mathrm{Aut}(G) : \sigma(e) = e \text{ setwise} \;\forall e \in S\}$.
Automorphisms relabel vertices, mapping NIC-1-planar drawings to
NIC-1-planar drawings; for the encoder's variable $x_e$ ("edge $e$ is
crossed in the drawing"), $x_e(\sigma \cdot D) = x_{\sigma^{-1}(e)}(D)$.

The SAT contract we rely on (justified for $C$ and $T$ in Section 5):

> **Branch-safety.** For any disjoint $N, P \subseteq E$, the formula
> $F_G \wedge \bigwedge_{e \in N} \neg x_e \wedge \bigwedge_{e \in P} x_e$
> is satisfiable iff $G$ admits a NIC-1-planar drawing with $x_e=\mathtt{false}$
> for $e \in N$ and $x_e=\mathtt{true}$ for $e \in P$.

## 3. The proof method

A **node** is a pair $(N, P)$ with $N \cap P = \emptyset$. Let
$\mathcal{D}(N, P)$ be the set of NIC-1-planar drawings of $G$ with
$x_e = \mathtt{false}$ for $e \in N$ and $x_e = \mathtt{true}$ for $e \in P$.
Throughout the recursion we maintain:
$(\ast)$ $P$ is a union of full $\mathrm{Stab}(N)$-orbits on $E \setminus N$.

**Anchor.** Query the oracle on $F_G \wedge \bigwedge_{e \in E} x_e$ ("Case
B"). If UNSAT, every drawing has at least one uncrossed edge.

**Lex-leader case-split.** At a node $(N, P)$ satisfying $(\ast)$, let
$H := \mathrm{Stab}(N)$, and let $\mathcal{O}_0, \dots, \mathcal{O}_{k-1}$
be the $H$-orbits of $E \setminus (N \cup P)$ ordered by their lex-smallest
representatives $r_0 < r_1 < \cdots < r_{k-1}$. Generate $k+1$ children:

$$
\begin{array}{rl}
\textbf{B-child}: & (N,\; P \cup \bigcup_{i<k}\mathcal{O}_i), \\[2pt]
\textbf{A-child } i\;(0 \le i < k): & (N \cup \{r_i\},\; P \cup \bigcup_{j<i}\mathcal{O}_j).
\end{array}
$$

Each $\bigcup_{j<i}\mathcal{O}_j$ is a union of $H$-orbits, and
$\mathrm{Stab}(N \cup \{r_i\}) \le H$, so $(\ast)$ propagates.

**Recursion.** Query the oracle on each child. Propagate UNSAT or SAT;
recurse on TIMEOUT children. Bottom out at a verdict at every leaf or a
configurable depth cap.

## 4. Soundness

**Lemma 1 (orbit cover).** Fix a node $(N, P)$ satisfying $(\ast)$. For every
$D \in \mathcal{D}(N, P)$, there exists $\sigma \in \mathrm{Stab}(N)$ such
that $\sigma \cdot D$ lies in the constraint set of at least one child.

*Proof.* Either every edge of $E \setminus (N \cup P)$ is crossed in $D$
(then $\sigma = \mathrm{id}$ realizes the B-child), or there is
$e^\star \in E \setminus (N \cup P)$ with $x_{e^\star}(D) = \mathtt{false}$.
Let $i^\star$ be the smallest index with $e^\star \in \mathcal{O}_{i^\star}$;
choose $\sigma \in H$ with $\sigma(e^\star) = r_{i^\star}$ (possible by
$H$-transitivity within $\mathcal{O}_{i^\star}$). For $e \in N$, $\sigma$
fixes $e$ setwise, so $x_e(\sigma\cdot D) = \mathtt{false}$. For $e \in P$,
$(\ast)$ gives $\sigma^{-1}(e) \in P$, so $x_e(\sigma\cdot D) = \mathtt{true}$.
For $e \in \bigcup_{j<i^\star}\mathcal{O}_j$: $\mathcal{O}_j$ is
$H$-invariant, so $\sigma^{-1}(e) \in \mathcal{O}_j$; minimality of
$i^\star$ implies $x_{\sigma^{-1}(e)}(D) = \mathtt{true}$, hence
$x_e(\sigma\cdot D) = \mathtt{true}$. Finally
$x_{r_{i^\star}}(\sigma\cdot D) = x_{e^\star}(D) = \mathtt{false}$. $\square$

**Theorem (refutation).** If the anchor returns UNSAT and every leaf of the
recursion rooted at $(\emptyset, \emptyset)$ returns UNSAT, then $G$ admits
no NIC-1-planar drawing.

*Proof.* Induction on recursion depth. Branch-safety yields $\mathcal{D}=\emptyset$
at every UNSAT leaf. At an internal node, suppose $D \in \mathcal{D}(N,P)$;
Lemma 1 gives a child with $\sigma \cdot D$ in its constraint set, contradicting
the inductive hypothesis. $\square$

## 5. Branch-safety audit for `oops`

We audit the components active under
`oops -unsat=1 -nic -up-sepcycles -Cno-transitive` against the source.

**(i) Sound graph-theoretic constraints.** Each clause is valid in *every*
NIC-1-planar drawing of $G$, hence preserved under any $(N, P)$:
strict crossing semantics (`encodeStrictConstraints`); K4 forbidden
2-clauses; sep-cycle 3-clauses; overfull-separator 2-clauses; IC-1
2-clauses; the dynamic sep-cycle user propagator (`createSepCyclesDynamic`,
emits only conflict clauses over cross1/cross2 literals corresponding to
geometrically infeasible partial assignments).

**(ii) Inactive on $C$ and $T$.** Each of the following is guarded by a
condition that fails on these graphs:

- Degree-3 candidate-pair reduction (in `initCrossablePairs`): guarded by
  $\mathrm{adj}(v) \subseteq \{x,y,z\}$ for a triangle-like configuration;
  impossible at girth $\ge 7$.
- Almost-twin reduction (in `initCrossablePairs`): requires two edges
  swap-feasible with shared neighborhood; impossible at girth $\ge 7$.
- Sep-cycle filter on candidates (in `initCrossablePairs`): on these
  graphs no candidate pair is filtered.
- Density / size early-exit (in `main.cpp`): require $n < 7$ or $m < 18$;
  both fail.
- Skewness early-exit (in `main.cpp`): returns SAT only if computed
  skewness $\le 1$; both graphs have skewness $\ge 2$.
- Degree-2 lex SB on cross1 (in `encodeCross1Constraints`): guarded by
  $\mathrm{adj}(v) = 2$; vacuous on 3-regular graphs.
- Swap constraints (`encodeSwapConstraints`): emits forbidden 2/3-clauses
  derived from a "swap reduces crossings" argument, but the requisite
  short-cycle / shared-neighbor configurations do not occur on these graphs,
  so no clauses are produced.
- Satsuma symmetry detection (`applySatsuma`): runs on the encoded CNF, but
  with strict-mode encoding present no CNF-level group is detected.
- BreakID: not enabled.

**(iii) Encoder features whose branch-safety is not established.**
The following are sound at the unconstrained root but not proven sound
under arbitrary $(N, P)$:

- **Symmetry-breaking / WLOG constraints vs. the external cross1
  assignment.** The encoder adds constraints that are sound only because
  they pick a *canonical representative* among assignments the solver would
  otherwise treat as equivalent — label-dependent stack-order symmetry
  breaking (`encodeStackSymmetry`: "vertex $0$ first", twin-vertex
  tie-breaks), the degree-2 lex constraint, CNF-level symmetry tools
  (Satsuma, BreakID), and any reduction justified "WLOG, by an automorphism
  of $G$". A `-fix-cross1` assignment is external and can break the very
  invariance such a constraint assumes: fixing an edge incident to only one
  member of a symmetric pair makes the two no longer interchangeable, yet a
  canonicalizing clause added blind to the fixing would still eliminate one
  representative — discarding the branch's only feasible model and yielding
  a false UNSAT. The encoder therefore **disables** its label-dependent
  canonicalizations (twin tie-break, degree-2 lex) whenever `-fix-cross1` is
  non-empty. The remaining obligation is to confirm that every other
  $\mathrm{Aut}(G)$-based reduction the encoder emits is either invariant
  under arbitrary cross1 fixings or likewise guarded.
- **`-Cno-transitive`**: drops transitivity clauses on auxiliary `rel`
  variables, weakening $F_G$ to $F'_G$ with $\mathrm{Models}(F_G) \subseteq \mathrm{Models}(F'_G)$.
  Sound for UNSAT only ($\mathrm{UNSAT}(F') \Rightarrow \mathrm{UNSAT}(F)$).
  On a SAT branch, `oops` asserts internally on inconsistent rel-layer
  models. A reviewer should detect the assertion in any branch and either
  re-run that branch without `-Cno-transitive`, or use the flag only
  when targeting UNSAT (as in the proof of $C$, $T$).

**Conclusion.** Under (i)+(ii), every component of `oops` either holds in
every drawing or is inactive on $C$ and $T$. With the label-dependent
canonicalizations of (iii) disabled under `-fix-cross1`,
$F_G \wedge \bigwedge_{e\in N}\neg x_e \wedge \bigwedge_{e\in P} x_e$
is satisfiable iff $\mathcal{D}(N, P) \ne \emptyset$ for any reachable
$(N, P)$ in the recursion on $C$ or $T$, subject to the residual obligation
in (iii) — that every remaining $\mathrm{Aut}(G)$-based reduction be
invariant under arbitrary cross1 fixings.

## 6. Implementation: the `-fix-cross1` option

`oops` accepts a fixed cross1 assignment via `-fix-cross1`, a
`;`-separated list of clauses, each a `,`-separated list of literals
`eN+` (cross1 of edge `N` true) or `eN-` (false). The case-split driver
only ever passes unit literals, so each clause fixes one edge.

Each `-fix-cross1` clause is added to the solver as a clause over the
corresponding cross1 literals; the option is inert when unset, so default
`oops` behavior is unchanged. The fixings are injected alongside the
encoder's $\mathrm{Aut}(G)$-based symmetry-breaking, which is built without
knowledge of them and may not be invariant under them — see Section 5(iii).

The driver (`scripts/oops_split.py`) computes $\mathrm{Aut}(G)$, walks the
recursion of Section 3, and invokes `oops` per node with the appropriate
`-fix-cross1` argument.

## 7. Independent verification

To verify a refutation produced by the method, a reviewer:

1. Reproduces $\mathrm{Aut}(G)$ and, at every recursion node, the
   $\mathrm{Stab}(N)$-orbits and lex-leader children of Section 3.
2. Confirms invariant $(\ast)$ holds at every node.
3. Re-runs `oops` with `-fix-cross1` on the leaf $(N, P)$ and confirms
   UNSAT (subject to the Section 5(iii) symmetry-invariance risk).
4. Confirms the inactivity claims in Section 5(ii) on the chosen graph
   either by inspecting the guards in source or by examining the per-leaf
   log for zero clauses from each listed component.

Steps 1, 2, 4 are mechanical; step 3 is the SAT verification.

## 8. Experimental results

All runs use `-unsat` direction (`-Cno-transitive`), `-node-timeout=300`,
48 threads; wallclocks include $\mathrm{Aut}$ computation and the anchor.
Every graph is 3-regular of girth $\ge 7$ and closes within 8 h.

| graph          | $\lvert V\rvert,\lvert E\rvert$ | $\lvert\mathrm{Aut}\rvert$ | variant | verdict                  | wallclock |
|-----------------|-----------|------------|---------|--------------------------|-----------|
| Coxeter         | 28, 42    | 336        | NIC     | not NIC-1-planar (UNSAT) | 24 min    |
| Tutte-Coxeter   | 30, 45    | 1440       | NIC     | not NIC-1-planar (UNSAT) | 37 min    |
| Tutte-Coxeter   | 30, 45    | 1440       | plain   | not 1-planar (UNSAT)     | 2.6 h     |
| cubic girth-7   | 30, 45    | 16         | plain   | not 1-planar (UNSAT)     | 6.5 h     |

The two Tutte-Coxeter rows are the same graph (`wiki_nonplanar.cfg` part 27)
decided in two variants; Coxeter is `wiki_nonplanar.cfg` part 25; the last
row is a less-symmetric cubic girth-7 graph (`3reg-gir7.g6` part 256).

Coxeter and Tutte-Coxeter are edge-transitive, so each root split has a
single A-child (one edge orbit) and the proof reduces to that one WLOG case;
the less-symmetric cubic graph ($\lvert\mathrm{Aut}\rvert=16$) has several
edge orbits and thus several root A-children. A richer automorphism group
helps within otherwise-comparable instances: at the same $|V|,|E|$,
plain-1-planarity Tutte-Coxeter ($\lvert\mathrm{Aut}\rvert=1440$) closes in
2.6 h versus 6.5 h for the cubic graph ($\lvert\mathrm{Aut}\rvert=16$),
because the larger group collapses more of the search into orbit
representatives.

On satisfiable (1-planar) instances the `-unsat` configuration cannot
produce a witness: a satisfiable branch either makes `oops` assert under
`-Cno-transitive` (reported as CRASH) or times out. To decide such
instances, run the `-sat` configuration (drops `-Cno-transitive`), which
returns SAT on some branch; case-splitting cannot turn a satisfiable
formula unsatisfiable, so this is sound.

## 9. Limitations

The audit relies on inactivity claims (girth $\ge 7$, 3-regularity,
skewness $\ge 2$) that hold for $C$ and $T$. Applying the method to
another graph requires re-checking each item in Section 5(ii). The
case-split is efficient only when $\mathrm{Aut}(G)$ acts non-trivially on
edges. `-Cno-transitive` is sound for UNSAT only.

Most importantly, the injected `-fix-cross1` assignment can break the
invariance that the encoder's $\mathrm{Aut}(G)$-based symmetry-breaking
constraints assume, producing a false UNSAT (Section 5(iii)). This is
inactive on the twin-free, girth-$\ge 7$ targets but live on denser graphs,
where it must be controlled by disabling or guarding every such constraint
under `-fix-cross1` (as done for the twin tie-break and the degree-2 lex
constraint).
