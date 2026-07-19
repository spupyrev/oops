# Cubic graphs through 28 vertices: progress report

This file is the durable research log for the verification project described in
`verify_cubic28_plan.md`.  Entries distinguish proved statements, hypotheses,
implementation status, and measurements.  A reduction is not production-ready
until its mathematical statement and implementation predicate agree and both
have been tested.

Current selection (2026-07-15): the project uses the five-audit architecture
recorded in the final entry below and in `cubic/verification.md`. Earlier
entries are retained as a chronological account of experiments and may
describe proof stages that were subsequently rejected or removed.

## 2026-07-10: initial audit

### Target and proof structure

The enumerated class is simple, connected, unlabeled cubic graphs.  Since cubic
graphs have even order, the new orders are 26 and 28.  The intended proof is by
increasing order: the established computation through order 24 is the base;
order 26 is verified next; then order 28 may use order 26 as part of its base.

For a skipped graph `G`, it is not enough to recognize a visually plausible
configuration.  We need a smaller graph `H` in an already verified class and a
proof that every 1-planar drawing of `H` can be locally expanded into a
1-planar drawing of `G`.  If `H` is not simple, connected, and cubic, the proof
must separately explain why it is covered.

### Reproducibility and compute constraints

- Research/development budget: 48 hours; final production runs are excluded.
- Acceptable final verification runtime: approximately 2--3 days on a stronger
  32-core machine.
- Current development machine: 14 logical CPUs, Intel Core Ultra 5 125U.
- Initial repository state: clean, commit `00eb86d` (`initial planning for
  cub28`).
- No `geng`, `genreg`, `plantri`, or `planarg` executable was found on `PATH`.
- Existing graph data is under
  `/home/spupyrev/research/one_planar/data/cubic`.

The published table is consistent with 32-core wall time.  For example,
4,624,501 girth-at-least-6 graphs at 1.7 seconds each represent about 91 CPU
days, or 2.84 days on 32 cores, close to the reported 2 days 23 hours.

### Existing reduction status

`hasReducibleTriangle` currently skips:

1. a triangle whose three external neighbors are pairwise distinct; and
2. a pair of triangles sharing an external edge pattern, subject to one
   additional distinctness test.

The corresponding proof in `main_cubic.tex` is incomplete.  The cut-pair and
cut-triple experiments are disabled.  Therefore the current
`-skip-reducible` mode is experimental and must not yet be used as verification
evidence.

### Initial ranked directions

1. **Complete and validate the triangle classification.**  This is already
   implemented in part, and triangle contraction has a clean local expansion
   when the three external neighbors are distinct.  Enumerate all equality
   patterns forced by simplicity and cubicity, state the exact reductions, and
   make the code report rule-specific counts.
2. **Analyze small edge cuts and separations.**  For cubic graphs, replacing a
   side of a 2- or 3-edge cut by a small completion may produce smaller cubic
   graphs.  The key issue is whether drawings can be composed without assuming
   an unsupported common-face property.
3. **Classify 4- and 5-cycle attachment patterns.**  Do not aim initially for a
   blanket girth theorem.  Seek individual patterns admitting local expansion,
   leaving an explicitly generated exceptional family.
4. **Improve the verification pipeline.**  Current g6 processing scans each
   input once to count and again to solve; modulo shards each scan the whole
   file.  Production must instead use generator-side shards or streaming and
   must fail on unknown results or count mismatches.
5. **Optimize SAT only after measuring irreducible survivors.**  The project
   needs large combinatorial reduction before constant-factor solver work can
   dominate.

### Immediate next experiment

Derive the complete isomorphism classification of a triangle's three external
neighbors in a simple connected cubic graph.  Compare it to
`hasReducibleTriangle`, write exhaustive small graph tests for the predicate,
and measure rule frequencies separately on available order-18/20/22 data.

## 2026-07-10: triangle reduction audit

### Result

A complete proof that a minimal subcubic non-1-planar graph is triangle-free
has been drafted in `main_cubic.tex`.  The proof now includes:

- a precise minimal-subcubic formulation, which connects the cubic enumeration
  to reductions that temporarily produce lower-degree vertices;
- biconnectivity and cubicity lemmas;
- a planar two-terminal replacement lemma; and
- all equality patterns among a triangle's three external neighbors.

The difficult pattern is a diamond (two triangles sharing an edge).  If its
outside endpoints are nonadjacent, the diamond is replaced by an edge.  If
they are adjacent, that adjacent pair can be absorbed into the planar
two-terminal patch.  Repeating this operation must reach either a cut vertex
or a nonadjacent terminal pair, because cubicity saturates every previously
absorbed vertex.  This avoids assuming that a minimal counterexample has
girth five or six.

The TeX source passes a `pdflatex` syntax build when the locally available
LIPIcs class is added through `TEXINPUTS`.  A later proof-review pass is still
required; compilation is not a correctness check.

### Instrumentation

The experimental reducer now classifies triangle configurations separately:

- `triangle-contraction`;
- `diamond-to-edge`;
- `diamond-adjacent-terminals`;
- `diamond-cut-vertex`;
- `K4`; and
- `unclassified-triangle` (a defensive completeness category).

The new `-analyze-reducible` option classifies a graph collection without
running planarity or SAT.  Final summaries report every category.  This does
not yet change the production verification policy.

### Measurements

All counts below are over complete connected cubic collections unless noted.
The analyzer found zero `unclassified-triangle` instances.

| order/data | total | no detected triangle | triangle contraction | diamond to edge | adjacent terminals | cut vertex |
|---|---:|---:|---:|---:|---:|---:|
| 20, all | 510,489 | 97,546 | 329,731 | 67,491 | 8,160 | 7,561 |
| 22, all | 7,319,447 | 1,537,240 | 5,209,343 | 499,872 | 47,351 | 25,641 |

The order-22 scan took 27 seconds after the input's initial counting pass.
The triangle-free, biconnected, nonplanar file with 1,432,712 graphs contains
no recognized triangle configuration, as expected.

The available order-26 archives appear to contain all 100 residue shards of
the triangle-free, biconnected, nonplanar family.  The first 30 uncompressed
shards total about 7.6 GB, suggesting that triangle elimination alone leaves
hundreds of millions of order-26 instances.  Thus triangle elimination is a
sound first reduction but is far from sufficient.

### Next direction

Analyze 4-cycles as four-terminal patches.  The immediate goals are to:

1. classify equality and adjacency patterns of their external neighbors;
2. identify patterns that collapse to the proved two-terminal replacement;
3. test whether contraction to one degree-four vertex, with a prescribed
   rotation, gives a useful computational certificate for the remaining
   patterns; and
4. measure pattern frequencies on the existing triangle-free order-22 and
   sampled order-26 data before committing to a complex implementation.

## 2026-07-10: 4-cycle split certificate

### Structural data

The `geng` executable from nauty 2.9.3 is available at
`/home/spupyrev/research/book_embedding/solvers/nauty2_9_3/geng`.  Its relevant
options are `-t` (triangle-free), `-f` (4-cycle-free), `-p` (5-cycle-free),
`-C` (biconnected), and canonical `res/mod` sharding.  Therefore a proved
reduction for a cycle length can be applied during generation rather than by
scanning the larger family.

The archived order-26 triangle-free, biconnected, nonplanar corpus has exactly
432,416,412 graphs: its archived uncompressed size is 24,647,735,484 bytes and
each order-26 graph6 record is 57 bytes.  The existing girth-at-least-5 count is
31,478,584, so about 401 million graphs contain a 4-cycle.

On all 1,432,712 triangle-free, biconnected, nonplanar order-22 graphs, the
4-cycle analyzer counted the following incidences (categories can overlap
at graph level because one graph may contain several patterns):

| pattern | cycles | graphs containing pattern |
|---|---:|---:|
| four distinct external neighbors, no external edge | 1,897,878 | 1,136,190 |
| four distinct external neighbors, some external edge | 1,337,958 | 660,280 |
| one repeated opposite neighbor, no external edge | 433,230 | 148,168 |
| one repeated opposite neighbor, some external edge | 109,773 | 52,435 |
| both opposite pairs repeated | 18,000 | 17,763 |

No chord/invalid pattern was found, and a girth-at-least-5 sanity dataset
produced zero 4-cycles.  In a 38,258-graph triangle-free order-26 sample,
36,540 graphs (95.5%) contain a 4-cycle.

### Certificate and implementation

For a chordless 4-cycle, replace the four cycle vertices by two adjacent cubic
vertices, pairing two consecutive external incidences at each new vertex.
This produces a cubic graph with two fewer vertices, including when opposite
external neighbors coincide.  If the new joining edge has an uncrossed
1-planar drawing, the edge can be expanded inside a crossing-free disk to
recover the original 4-cycle.  The rotations may force two new opposite cycle
edges to cross each other once, which remains a valid 1-planar drawing.

The implementation:

- generates both consecutive pairings for the first suitable 4-cycle;
- forbids crossings on the new joining edge in OOPS;
- verifies the reduced SAT witness with the existing planarization checker;
- accepts the original graph by the proved local expansion; and
- falls back to the original graph if neither reduced certificate succeeds.

The fallback makes this a sound solver optimization without assuming that all
reduced instances succeed.

### Pilot results

- 96 order-26 graphs: 87 initially eligible generic reductions, all certified;
  mean time 129 ms versus 166 ms baseline.
- 479 order-26 graphs after generic support: 452 certified, zero failed
  reductions; mean 106 ms versus 160 ms baseline.
- After supporting opposite repeated neighbors, 460 of 479 graphs were
  certified, zero failed reductions; mean 100 ms.
- Disabling the extra SAT-oriented constraints (`-sat=0`) improved a matched
  96-graph run from 102 ms to 88 ms.
- None of the reduced graphs in a 479-graph pilot was planar, so a plain
  planarity shortcut does not explain the success.
- The brute-force and move encoders were substantially slower and were
  abandoned after short pilots.  Raising the preliminary crossing search from
  one to two crossings increased mean time to 636 ms.

Despite the consistent certificate rate, roughly 90--100 ms per 4-cycle graph
would still require around two weeks on 32 comparable cores for the complete
order-26 4-cycle-containing corpus.  This reduction is useful but does not by
itself meet the reproducibility target.

### Uncrossed-edge conjecture

Every edge of each of five sampled nonplanar order-26 graphs admitted a
constrained 1-planar drawing (195/195 edge tests).  Exhaustively, all edges of
all 376 nonplanar connected cubic graphs of order 14 also admitted such a
drawing; there were no timeouts.  This is evidence, not a proof, that a
subcubic 1-planar graph may admit an uncrossed drawing for any prescribed
edge.  A proof would turn the split certificate into a direct 4-cycle
reduction.  The obvious rerouting argument is insufficient: detouring a
crossed edge around an endpoint of its crossing partner can cross the other
two incident edges at a cubic vertex.

### Next direction

Search for either:

1. a rigorous uncrossed-edge lemma for subcubic 1-planar graphs (or a weaker
   lemma tailored to the split graphs); or
2. a counterexample using the new exhaustive `-analyze-uncrossed-edges` mode.

In parallel, derive the analogous 5-cycle replacement gadget and determine
which small forest must be uncrossed for its local expansion.  This will show
whether a common constrained-drawing lemma can eliminate both cycle lengths.

## 2026-07-10: 5-cycle certificate audit and encoder speedup

### 5-cycle split: corrected certificate

Replacing a chordless 5-cycle by a three-vertex path with attachment groups
2+1+2 gives a smaller cubic graph.  Requiring both path edges to be uncrossed
provides a crossing-free disk in which to attempt the reverse expansion.
However, the first version of this certificate was too broad.  An explicit
enumeration of the eight rotation systems of the path gives boundary orders

```text
04321 03421 02431 02341 01432 01342 01243 01234.
```

For orders `02431` and `01342`, the naïve chord drawing crosses the same
5-cycle edge twice.  The implementation no longer accepts a reduced SAT
answer merely because the two path edges are uncrossed.  It now:

1. requests the concrete two-page witness from OOPS;
2. reconstructs the cyclic rotation at each of the three path vertices;
3. traverses the boundary of the path's regular neighborhood;
4. computes the induced order of the five external incidences; and
5. accepts only if drawing the 5-cycle as chords crosses every new edge at
   most once.

This is a directly checkable sufficient condition, stated as the 5-cycle
split certificate in `main_cubic.tex`.  A reduced witness with a bad rotation
is discarded and the next split (or the original graph) is tried.  On the
matched 479-graph order-26 pilot, all 18 graphs left after the 4-cycle
certificate still received an accepted 5-cycle certificate.  Together the
two reductions certified all 479 graphs with no fallback; mean processing
time was 56 ms after the encoder optimization below.

On 91 sampled order-22 girth-at-least-5 graphs, 89 had a 5-cycle reduction
and all 89 passed the corrected rotation check (the other two appear to have
girth at least six).  With canonical transitivity and `-Cno-dedup`, mean time
was 32 ms, versus 67 ms for the earlier direct-solve baseline.

### Evidence for edge flexibility

The uncrossed-edge search was extended beyond the initial order-14 run:

- every edge was flexible in all 4,060 connected cubic graphs of order 16
  (681 planar and 3,379 nonplanar graphs);
- every edge was flexible in stratified samples of 104 order-18 graphs, 103
  order-20 graphs, and 103 triangle-free nonplanar order-22 graphs; and
- there were no forced-edge results or timeouts.

This is still experimental evidence, not a theorem.  In particular, the
obvious endpoint-detour argument does not prove the claim for subcubic
1-planar graphs.

### Canonical transitivity clauses

The stack encoder previously emitted each transitivity clause three times,
once for every cyclic rotation of its vertex triple, and relied on global
clause sorting and deduplication to remove the copies.  It now emits only the
rotation whose first vertex is the smallest.  Swapping the other two vertices
still emits the second logical orientation, so the formula is unchanged.

For a representative reduced order-24 instance, clause construction fell
from approximately 205,000 pre-deduplication clauses to approximately 75,000.
On a matched 96-graph order-26 4-cycle pilot with `-sat=0`:

| configuration | mean time |
|---|---:|
| before canonical transitivity | 88 ms |
| canonical generation, normal global deduplication | 65 ms |
| canonical generation, `-Cno-dedup` | 55 ms |

After canonical generation, only about 500 duplicate clauses remain in the
profiled instance, so skipping the global sort/dedup pass is beneficial.
The default remains deduplicated pending broader regression testing;
`-Cno-dedup` is the candidate production setting.  Setting `-skewness=0`
saved only about 1 ms in the same pilot.  Several short, matched experiments
were rejected: true polarity for relation variables, natural variable order,
phase saving zero, random activity initialization, disabling solver
simplification, the move encoder, and brute force were neutral or slower.

At 55 ms per 4-cycle-containing order-26 graph, the roughly 400.94 million
4-cycle corpus projects to about 8.0 days on 32 cores of this development
machine.  A 2.5x per-core speedup on the intended server would reduce that to
about 3.2 days.  This is now close enough that an actual server pilot and
generator-side sharding matter more than speculative micro-optimization, but
the complete order-28 projection still depends on the 5-cycle survivor
counts and timings.

### Immediate validation work

Before treating these changes as verification infrastructure:

1. independently review the stack-rotation extraction against the explicit
   two-page drawing convention;
2. measure certificate failures (including bad rotations) on a substantially
   larger girth-five sample; and
3. build the generator-side sharding and result-audit pipeline.

The first three validation tasks from the previous checkpoint are complete.
`data/cubic_reductions.g6` and `scripts/run_tests.sh` now cover 18 4-cycle
certificates, two 5-cycle certificates, and simultaneous one/two-edge
constraints.  The complete existing test suite passes.  The five documented
NIC non-1-planarity smoke instances also pass after canonical transitivity,
with zero unknown results.  `main_cubic.tex` passes a syntax build.

## 2026-07-10: constrained-forest certificates and composite reductions

### Corrected 5-cycle topology

The earlier rotation filter for the 5-cycle split was unnecessarily
restrictive.  For boundary orders `02431` and `01342`, drawing only the cycle
as chords crosses the closing edge twice.  The full patch also contains the
five attachment spokes: in both exceptional orders, crossing the closing edge
with the middle spoke gives a drawing with exactly one crossing.  Therefore
all eight path rotation systems expand 1-planarly once the two path edges are
uncrossed.  `main_cubic.tex` now states the unconditional certificate, and the
implementation no longer extracts or filters stack rotations.

### Incremental edge and path coverage

The first edge-flexibility analyzer rebuilt the entire CNF for every crossed
edge.  It now adds one selector `U_e` per edge with clauses saying that `U_e`
forces the edge's division vertex to be comparable with every other division
vertex.  Assuming `U_e` is therefore exactly an uncrossed-edge constraint,
without enabling the larger cross1/cross2 encoding.

One base witness covers every edge that it leaves uncrossed.  The remaining
edges are tested with assumptions on the same solver, retaining learned
clauses.  A group-testing step first assumes a whole group simultaneously and
recursively splits only an UNSAT group.  The same mechanism audits every pair
of adjacent edges: a base witness covers pairs of uncrossed edges, and each
remaining path query assumes both edge selectors.

Measurements:

- a systematic 1,033-graph order-18 sample took 156 seconds with rebuilt
  formulas, 98 seconds with the first incremental cross-variable version, and
  substantially less per representative graph with direct selectors;
- on 96 sampled order-26 triangle-free nonplanar graphs, all edges and all
  adjacent edge pairs were flexible, with no unknowns;
- adjacent-pair coverage averaged 88 ms per order-26 graph, essentially the
  same as the direct SAT solve, because almost every graph was covered by the
  base witness plus one grouped alternate witness.

This finite property is stronger than ordinary 1-planarity and would certify
all 4-cycle and 5-cycle expansions from order 28 to order 26.  Auditing all
roughly 432 million triangle-free order-26 graphs is plausible but currently
projects above the preferred three-day production budget on the development
machine.  It remains a fallback proof route, not yet the selected production
flow.

### Composite 4-cycle reduction

A 4-cycle split can create triangles in its smaller graph.  The reducer now
composes the split with sound local reductions:

1. a triangle with three distinct external neighbors is contracted;
2. if the distinguished split edge lies inside that triangle, its constraint
   disappears, because expanding the contracted vertex creates it inside a
   crossing-free disk; and
3. a diamond with nonadjacent outside terminals is replaced by an edge.  An
   internal distinguished edge again becomes automatically uncrossed; a
   boundary distinguished edge maps to the replacement edge.

These operations are repeated until no applicable triangle remains.  On the
38,258-graph order-26 sample, the original split was immediately triangle-free
for 36,068 of 36,540 4-cycle graphs.  Ordinary triangle contraction increased
this to 36,514, and diamond composition increased it to all 36,540.

On the complete 1,432,712-graph order-22 triangle-free biconnected nonplanar
corpus, 1,341,775 graphs contain a 4-cycle and 1,341,651 receive a triangle-free
composite reduction.  The remaining 124 contain nested false-twin / `K_{2,3}`
modules.  This isolates the outstanding 4-cycle theorem to a small structural
family, but no general proof for that family has yet been claimed.

### Exact girth counts and the 6-cycle certificate

The independent [House of Graphs cubic census](https://houseofgraphs.org/Cubic)
gives the exact order-28 counts:

- girth at least 4: 8,542,471,494;
- girth at least 5: 656,783,890; and
- girth at least 6: 4,624,501.

Thus a 4-cycle theorem is the decisive combinatorial reduction.  Even after
removing 4-cycles, direct girth-six solving remains costly: 21 sampled graphs
averaged 1.65 seconds, matching the historical 1.7-second figure.

For a chordless 6-cycle, replace the cycle by a four-vertex path with
attachment groups 2+1+1+2 and require all three path edges to be uncrossed.
There are sixteen local rotation systems.  An exhaustive anchored-planarity
check found a 1-planar patch for every one, using at most two pairwise
edge-disjoint crossings.  The full sixteen-case certificate is recorded in
`main_cubic.tex` and implemented as `-six-cycle-reduction`.

On the same 21 order-28 girth-six graphs, every reduced certificate succeeded
without fallback.  The best fixed split offsets averaged 0.56--0.58 seconds;
the default offset averaged about 0.58 seconds, a 2.8x improvement.  Sorting
the six splits by reduced girth was worse (1.01 seconds) and was reverted.
Although this is still slower than the 4-cycle and 5-cycle reductions, it
projects to hours rather than days on the intended server.

### Current proof architecture

The most promising complete flow is now:

1. use the proved minimal-subcubic argument to remove disconnected graphs,
   cut vertices, and triangles;
2. finish the nested false-twin case so every 4-cycle has a sound composite
   reduction to an order-26 graph with a one-edge obligation;
3. choose between an order-26 edge/path-flexibility audit and direct
   order-28 5-cycle certificates based on an actual 32-core server pilot;
4. use the 6-cycle certificate on the 4.6 million girth-six graphs; and
5. solve the 21 girth-at-least-seven graphs directly.

The open mathematical point is now narrow and explicit: a sound treatment of
the 124 order-22-style nested `K_{2,3}` 4-cycle exceptions, or an alternative
proof that the necessary split edge is flexible in that family.

## 2026-07-10: closing the order-22 4-cycle exceptions

The count of 124 above was an intermediate checkpoint.  Three additional
composable replacements were implemented and checked structurally:

1. a six-vertex false-twin ladder is a planar two-terminal patch and can be
   replaced by an edge;
2. a false-twin `K_{2,3}` with three distinct outside neighbors can be
   contracted to one cubic vertex; its reverse patch has a 1-planar drawing
   in any of the three-terminal rotations; and
3. any side of a two-edge cut is replaceable when adding an edge between its
   internal terminals makes that side planar.  This is an explicit
   cofaciality certificate for the planar two-terminal replacement lemma.

The reducers update or discharge the distinguished uncrossed edge while they
compose.  They reduced the 124 exceptions first to 25, then 16, and finally
four.  The four remaining graphs have archive indices 599299, 619422, 823396,
and 1226747.  In every case the composite reduction itself is planar; direct
checks with `-four-cycle-reduction` certify all four without invoking SAT on
the original graph (16--25 ms per graph on the development machine).

Consequently, the analyzer's earlier "triangle-free reduction" metric was too
narrow for measuring closure into the finite census: a triangle-free smaller
instance carries the remaining uncrossed-edge obligation, while a planar
smaller instance discharges it immediately.  It now reports this union as a
"census-ready" 4-cycle reduction.  The repeated complete order-22 scan
finished in 87 seconds:

- graphs processed: 1,432,712;
- graphs containing a 4-cycle: 1,341,775;
- graphs with a census-ready triangle-free or planar composite: 1,341,775; and
- 4-cycle graphs without such a composite: zero.

This closes the observed nested-`K_{2,3}` family at order 22.  It is strong
evidence for the 4-cycle decomposition, but the code paths and the local
three-terminal `K_{2,3}` drawing still require a clean formal lemma and an
independent implementation review before the result can be used as a proof.

## 2026-07-10: correction of coverage and local-cycle certificates

### Glucose assumptions were not enforced

The earlier incremental edge/path coverage numbers in this report are
invalid.  The vendored Glucose `solveLimited(assumptions)` copied assumptions
but its search loop never applied them.  Adding the standard assumption loop
exposed a second repeated-solve defect: after a solve ended under CHB or
distance branching, `cancelUntil(0)` repopulated only the active decision
heap; the next solve reset to VSIDS, saw a stale empty heap, and returned a
partial assignment as SAT.

Both defects are now fixed.  Every solve rebuilds all decision heaps, and
assumption solves disable the fork's nonstandard chronological-backtracking
shortcuts.  A three-variable regression in
`cubic/test_solver_assumptions.cpp` first solves a base formula and then
checks that contradictory assumptions are UNSAT.  For coverage, every SAT
assumption is checked in the model, a full drawing is reconstructed, and its
planarization is independently tested.  The solver also asserts that every
active original clause is satisfied before accepting an assumption model.

The coverage algorithm no longer assumes the union of a large group.  It
maintains the crossing sets of verified witnesses, finds a smallest set of at
most one/two/three edges that hits all current crossing sets, asks for a
drawing avoiding that set, and repeats.  Once no such hitting set exists, the
stored drawings themselves prove the requested flexibility property.

Corrected measurements (all with zero forced sets; stated samples also had
zero unknowns unless noted):

- 100 stratified order-24 graphs, all edge pairs: 0.571 s/graph;
- 20 stratified order-24 graphs, all edge triples: 1.437 s/graph;
- 111 stratified order-26 girth-five graphs, all individual edges:
  0.766 s/graph; and
- 18 order-26 girth-six graphs, all individual edges, using the static
  cross-check implementation: 3.15 s/graph.

The much faster 88--97 ms results recorded earlier were unconstrained repeat
solves and must not be used for production projections.

### Attachment-spoke flaw in the 5-cycle and 6-cycle tables

The earlier "all rotations" local certificates were also too broad.  Several
table entries crossed a new cycle edge with an attachment spoke.  That spoke
continues an old edge outside the replacement disk; if the old edge already
has a crossing, the expanded edge would be crossed twice.

The implementation now reconstructs the cyclic boundary order of the path's
regular neighborhood from the concrete two-page witness.  It accepts only
orders whose local drawing crosses cycle edges with cycle edges:

- six of the eight 5-cycle orders (excluding `02431` and `01342`); and
- eight of the sixteen 6-cycle orders.

The rotation reconstruction is explicit: it obtains the cyclic order of
incident segments from their spine positions and arc pages, then performs a
face walk of the ribbon path.  `main_cubic.tex` states only these conservative
certificates.  New order-26 5-cycle and order-28 6-cycle fixtures exercise
the witness-dependent checks.

Corrected pilots remain positive:

- all 111 order-26 girth-five graphs found an accepted 5-cycle rotation with
  no fallback, at 111 ms/graph (previously 88 ms without the filter); and
- all 21 order-28 girth-six graphs found an accepted 6-cycle rotation with no
  fallback, at 1.92 s/graph (previously about 0.58 s).

### Flexibility hierarchy under evaluation

Let `P_k(n)` mean that every set of at most `k` edges in every relevant graph
of order at most `n` can be simultaneously uncrossed.  A 4-cycle split maps a
prescribed `k`-set to at most `k+1` edges in a graph two vertices smaller.
Consequently, a complete `P_2(24)` audit would prove `P_1(26)` for the
4-cycle-containing order-26 family and would eliminate all order-28 4-cycles.

The remaining order-26 girth-at-least-five family has 31,478,584 graphs and
still needs a direct `P_1` audit or a stronger reduction.  On the corrected
development timings, the pair audit through order 24 and the direct
girth-five edge audit together are above the desired 2--3 day server budget;
the safe 6-cycle stage adds roughly another day under the conservative 2.5x
per-core server-speed assumption.  This is therefore a rigorous fallback
architecture, not yet the final selected production flow.

### Higher flexibility is not a blanket shortcut

The arbitrary-set analyzer was extended through sets of twelve edges.  On the
first triangle-free, biconnected, nonplanar order-22 graph, `P_4` and `P_5`
hold, but `P_6` already fails: after 56 independently verified drawings the
solver found a six-edge set that every drawing must meet.  Increasing the
requested bound through nine therefore fails on the same obstruction.  The
hoped-for one-line hierarchy based on `P_9(22)` is false and will not be used.

An isolated optimized-build check gave 0.803 s/graph for pair flexibility on
100 stratified order-24 graphs and 0.904 s/graph for edge flexibility on 111
stratified order-26 girth-five graphs.  These jobs ran simultaneously and are
not directly comparable with the isolated standard-build figures above, but
they show no release-build speedup large enough to change the architecture.
The project's `CHECK` validations remain active under `-DNDEBUG`; only C/C++
`assert` checks are removed.  Production evidence should therefore retain the
explicit witness reconstruction and planarization checks regardless of build.

### Reverse-census direction

The next candidate is to enumerate smaller graphs with marked paths and cover
all *inverse* 5-cycle/6-cycle expansions of each smaller graph.  A single
verified drawing can certify many marked paths, whereas processing the
order-28 census separately rebuilds and solves essentially the same smaller
graphs many times.  For a girth-five parent, every triangle or 4-cycle of the
smaller graph must intersect the distinguished replacement path; otherwise it
would survive in the parent.  Thus the smaller marked graphs are much more
restricted than the full order-26 census.  This route requires a rotation-aware
coverage analyzer and an exact census/generator argument before it can replace
the conservative direct flow.

## 2026-07-10: packed cycle contraction

### An initially unsound shortcut and its correction

Contracting a 5-cycle to one degree-five vertex, or a 6-cycle to one
degree-six vertex, makes a much smaller unconstrained SAT instance.  The first
prototype incorrectly treated every rotation at that vertex as expandable.
That repeats the topology error of the earlier path table: contracting the
cycle forgets the cyclic order of its attachments, and some boundary orders
force a cycle edge to cross twice unless an attachment spoke is used.  The
unrestricted 17 ms 5-cycle timing was therefore not a valid certificate and
is not used.

The corrected implementation reconstructs the rotation at every contracted
vertex from the concrete two-page witness.  An exhaustive local planarization
check gives all twelve 5-cycle boundary orders.  Six use cycle edges only;
five more have alternatives using one attachment spoke; the remaining order
has five alternatives using a neighboring pair of spokes.  A spoke-using
patch is accepted only when the corresponding old edge is uncrossed in the
outside witness.  When several cycles are contracted, a small backtracking
check also ensures that an edge between two contracted vertices is not used
by both endpoint patches.  Thus every accepted result has a direct composed
1-planar drawing.

The finite table and all of its alternatives are checked by
`cubic/test_cycle_expansions.cpp`.  For each crossing pattern the
test constructs its planarization, surrounds the labeled boundary vertices by
a wheel, and runs the independent planarity checker.  The wheel forces the
boundary rim to be a face, which verifies that the planarization is drawable
inside the replacement disk.  The table and the composition lemma are stated
in `main_cubic.tex`.

### Packing and measurements

The implementation chooses a deterministic low-conflict packing of
vertex-disjoint chordless cycles.  Cycles joined by more than one edge are not
packed together, avoiding parallel edges after contraction.  Contracting
several cycles composes because their vertex disks are disjoint and reduces
the SAT order by four vertices per 5-cycle.  Relabeling the same contracted
graph gives a cheap deterministic retry when the solver's first witness has
an unusable rotation; every retry is independently reconstructed and checked.

A minibaum shard supplied 9,483 order-28 girth-at-least-five graphs.  On a
stratified 1,054-graph subset, the packing distribution was:

- one 5-cycle: 4 graphs;
- two 5-cycles: 405 graphs;
- three 5-cycles: 636 graphs; and
- four 5-cycles: 9 graphs.

With at most three relabeled attempts, a complete 9,483-graph pilot certified
9,470 graphs by packed contraction and sent only 13 to the ordinary
original-graph fallback.  It used 10,805 contracted solves and finished in
195 seconds, or 20 ms per parent, with zero unknown or non-1-planar results.
The packing distribution was 31/3,578/5,838/36 graphs with respectively
one/two/three/four contracted 5-cycles.  These runs use
`-sat=0 -skewness=0 -Cno-dedup`; disabling clause deduplication is sound and
reduced construction overhead.

At 20 ms/graph, the 652,159,389 order-28 graphs of girth exactly five project
to about 151 CPU-days, or 4.7 days on 32 development-speed cores.  Under the
previously used conservative 2.5x per-core server factor this is about 45
hours.  This is now compatible with a two-to-three-day production budget,
subject to a broader server pilot and tail measurements.

The analogous packed 6-cycle contraction certified all 100 girth-six graphs
in a 101-record stratified sample (the other record has girth seven).  Before
adding spoke-conditioned 6-cycle orders it uses only thirteen cycle-edge-only
rotations, so its current 83 ms mean is conservative.  Even this projects to
only a few hours for the 4,624,501 graph girth-six census on the target server.

### Second correction: SAT witnesses do not determine adjacent-edge rotations

The rotation-dependent packed results immediately above are superseded.  The
stack encoder constrains crossings of nonadjacent segments but deliberately
does not encode a rotation among segments sharing a vertex.  Consequently,
sorting the displayed arcs after the solve does not recover a rotation that
is part of the SAT certificate.  An attempted expanded-witness reconstruction
caught this: several locally accepted patches produced a nonplanar combined
planarization.  The 9,483-graph timing remains useful as a performance
experiment, but its rotation-dependent acceptance rule is not a proof.

The implementation now avoids rotation inference entirely.  The complete
5-cycle table has the following stronger property: for every one of its
twelve boundary orders, some patch uses only cycle edges and a subset of the
two consecutive spokes `s0,s1`.  Thus requiring just those two old attachment
edges to be uncrossed makes the contraction expandable for *every* possible
rotation.  Packed cycles are now required to have no edge between them, so
their spoke obligations are disjoint.  No post-solve rotation is consulted.

Forcing all five spokes was first tested as a conservative version; it created
minute-long SAT tails and was rejected.  With only the two proved necessary
spokes, a 1,054-graph order-28 girth-five sample finished in 23 seconds
(21 ms/graph).  Every graph was certified on its first contracted solve, with
packing counts 4/787/263 for one/two/three mutually nonadjacent 5-cycles.
This restores the roughly two-day target-server projection for the girth-five
band on a rotation-independent certificate.

The 6-cycle contraction and path options are disabled in the binary pending
an analogous rotation-independent table.  Girth-six therefore falls back to
ordinary direct OOPS solving in the current rigorous flow; at the historical
1.7 seconds per graph this is about 2.8 days on 32 development-speed cores,
or roughly 27 hours under the conservative 2.5x target-server factor.

### Packed witnesses do not yet accelerate edge flexibility

To attack `P_1(26)`, packed contractions were varied adaptively.  After each
witness, the next packing avoids an uncovered edge and additionally requires
that preserved edge to be uncrossed.  A witness conservatively certifies only
uncrossed reduced edges outside all contracted cycle edges and the two local
spoke obligations.  This argument is rotation-independent and each ordinary
fallback drawing is checked.

The optimization is sound but not faster.  With ten packed attempts, only
four of six pilot graphs completed without fallback.  With six attempts all
six fell back, requiring 17 original constrained solves in total and averaging
0.91 s/graph, versus 0.766 s/graph for the existing incremental selector
audit.  The adaptive packed analyzer remains available for research, but is
not part of the selected production flow.  `P_1(26)` therefore remains the
dominant runtime stage and the principal target for further improvement.

## 2026-07-11: rotation-independent refinements and coverage search bias

### Proper-crossing regression strengthened

The local patch regression previously planarized a displayed crossing as an
ordinary degree-four vertex. Plain planarity does not by itself require the
two crossed curves to alternate at that vertex. The regression now expands
each crossing to a 4-cycle with the four endpoints attached in prescribed
alternating order. Contracting that 4-cycle recovers a proper crossing,
not a tangential four-way contact. All 5-cycle and conservative 6-cycle
table entries still pass. `main_cubic.tex` now describes this stronger finite
check.

Two small exhaustive search programs record further finite topology checks:

- `scripts/search_square_crossed_split.py` enumerates the eight boundary
  orders induced when the distinguished edge of a 4-cycle split is crossed.
  None admits a 4-cycle patch when old attachment edges are forbidden to gain
  another crossing. Thus the uncrossed distinguished-edge hypothesis in the
  4-cycle lemma cannot simply be dropped.
- `scripts/search_cycle_universal_spokes.py` enumerates all cycle boundary
  orders and proper crossing matchings for a prescribed set of usable spokes.
  It reproduces all twelve 5-cycle orders for `s0,s1`.

### Selector-biased incremental coverage

Profiling corrected `P_1(26)` showed that model construction takes only about
10--12 ms. Almost all runtime is in the six to fourteen incremental
assumption solves. A new optional heuristic keeps every proof obligation
unchanged but, after the initial witness, disables phase saving, assigns the
coverage selectors true-first polarity, and raises their decision activity.
The required target edge is still an explicit SAT assumption. Every returned
witness is reconstructed, checked by independent planarization, and checked
to leave the assumed edge uncrossed. Hence the change affects search order
only, not the certificate or soundness argument.

The option is `-coverage-selector-bias`; the activity magnitude can be varied
with `-coverage-selector-boost` (the selected default remains `1000000`). On
109 stratified order-26 girth-five graphs it proved complete edge flexibility
at 0.361 +/- 0.049 seconds per graph, versus the previous corrected 0.766
seconds. All 109 audits completed with zero forced or unknown edges. On 100
stratified order-24 graphs, complete edge-pair flexibility took 0.419 +/-
0.070 seconds per graph, versus 0.571 seconds without the bias.

Batching several uncovered edges into one assumption set was also tested.
Batch size two gave at most a few percent improvement on small samples and
was slightly slower on the matched 18-graph set; sizes four and eight created
bad tails (the size-eight sweep was stopped after 30 seconds). Production
therefore keeps singleton assumptions plus selector polarity bias.

### A possible lower-order hierarchy, not yet selected

A complete `P_3(22)` audit would handle most `P_2(24)` obligations produced
by another 4-cycle split. On 101 stratified order-22 graphs, selector-biased
triple coverage completed in 0.654 +/- 0.125 seconds per graph. Auditing the
whole 1,432,712-graph order-22 corpus would therefore be cheap.

This is not yet a replacement for full `P_2(24)`: if the prescribed pair in
the order-24 graph consists of two adjacent edges of the reduced 4-cycle, no
choice of the two ordinary 4-cycle splits makes both edges internal to the
replacement disks. A separate finite certificate or a direct audit of those
marked exceptions is needed. Earlier informal projections that simply
replace all of `P_2(24)` by `P_3(22)` would be incomplete; the selected flow
continues to use direct `P_2(24)`.

### Adjacent packed 5-cycles are safe with relabeled obligations

Forbidding every edge between two packed 5-cycles was stronger than needed.
The packer now permits at most one joining edge per pair, which preserves a
simple contracted graph. It accepts a packing only if every selected
5-cycle has two consecutive attachment positions not leading to another
selected cycle, and rotates the local labels so those positions are `s0,s1`.
The runtime assertion verifies that every forced spoke ends outside the
packing. Therefore every inter-5-cycle edge is used by neither local patch,
and it cannot acquire a local crossing at both endpoints.

On the established 1,054-graph order-28 girth-five sample, the packing counts
changed from 4/787/263 to 4/405/645 graphs with one/two/three 5-cycles. All
graphs were certified on the first reduced solve. Runtime fell from 21 ms to
15 ms per graph. Trying eight greedy starting choices found the same packing
distribution and added about 1 ms, so the single deterministic
minimum-conflict packing was retained.

### Stronger 5-cycle and 6-cycle tables rejected on runtime

The finite search found two additional rotation-independent facts:

- With four consecutive 5-cycle spokes available, every cycle edge and
  every one of the four spokes can individually be avoided by a valid patch
  for all twelve boundary orders. This permits one reduced `P_1` audit to
  certify all original edges. In practice it averaged 0.745 seconds on six
  graphs, versus 0.514 seconds for direct biased coverage, and was rejected.
- For 6-cycles, two, three, and four fixed consecutive spokes cover 32, 43,
  and 52 of the 60 boundary orders. Five fixed spokes cover all 60, with
  explicit proper-crossing matchings found by the search script. The
  resulting SAT instances have severe tails: with a two-second cap, five of
  the first thirteen stratified girth-six instances timed out. The run was
  stopped, and both 6-cycle contraction flags remain disabled. Direct OOPS
  solving remains the rigorous selected girth-six route.

### Updated conservative runtime projection

Using development-machine means and the same deliberately conservative 2.5x
per-core advantage for the 32-core target server gives approximately:

| stage | instances | development mean | projected target wall time |
|---|---:|---:|---:|
| `P_2(24)`, direct | 23,751,969 | 0.419 s | 1.44 days |
| `P_1(26)`, girth at least five | 31,478,584 | 0.361 s | 1.64 days |
| order-28 girth exactly five, packed contraction | 652,159,389 | 0.015 s | 1.42 days |
| order-28 girth at least six, direct | 4,624,501 | 1.7 s historical | 1.14 days |

The conservative sum is about 5.6 target days, before small generation and
filtering costs. This is a substantial reduction from the preceding flow but
still above the desired two-to-three-day reproduction window. The estimate
uses historical per-core scaling rather than a measurement on the actual
server; a short matched server pilot is required before treating it as a
schedule. The principal remaining mathematical target is a sound treatment
of the adjacent marked-pair exceptions that would allow the cheap `P_3(22)`
hierarchy to replace most of direct `P_2(24)`.

## 2026-07-11: closing the adjacent-pair 4-cycle exception

The adjacent marked-pair gap in the `P_3(22)` hierarchy has a small
rotation-independent contraction certificate. Contract one 4-cycle to a
degree-four vertex. Up to reversal, its boundary order is one of `0123`,
`0132`, and `0213`. To preserve the adjacent target pair `c0,c1`, the local
patches are respectively empty, `c3 x s1`, and `c2 x s0`. Thus a reduced
drawing in which all four spokes are uncrossed supports, by rotating labels,
a separate expansion preserving each of the four adjacent 4-cycle-edge pairs.
The proper-crossing regression now checks these three 4-cycle patches as well.

This gives a complete two-level architecture for `P_2(24)`:

1. audit `P_3` on all 556,471 connected cubic graphs through order 20 and on
   the 1,432,712 biconnected nonplanar triangle-free order-22 cores, for
   1,989,183 records total;
2. for each of the 23,751,969 order-24 core graphs, contract one 4-cycle and
   solve the 21-vertex graph with all four spokes forced uncrossed; and
3. on the 1,620,479 4-cycle-free (girth-at-least-five) graphs, use direct
   selector-biased `P_2` instead.

For every prescribed pair not consisting of adjacent edges on the selected
4-cycle, an ordinary 4-cycle split makes every prescribed 4-cycle edge internal
to an endpoint disk and maps the outside prescribed edges plus the middle
edge to a set of at most three edges at order 22. The contraction certificate
handles exactly the excluded adjacent pair. This proof is now stated in
`main_cubic.tex`, and the new analyzer is
`-analyze-square-pair-exceptions`.

The broadened lower-order base avoids assuming that composite reductions
preserve arbitrary triples. The new flexibility-core lemma proves directly
that blocks, degree-two suppression, triangle contraction, and planar
two-terminal replacement map a prescribed set without increasing its size.
Thus every non-core order-22 graph reduces to the fully audited order-at-most-
20 base, while the triangle-free biconnected nonplanar order-22 graphs are
audited explicitly.

On 100 stratified order-24 graphs, 96 4-cycle contractions and four direct
fallbacks averaged 65 ms/graph. A broader 476-graph sample contained 435
4-cycle contractions and 41 fallbacks. All 476 completed without a timeout
and averaged 0.178 seconds, but one 4-cycle-free fallback alone took about 25
seconds; a ten-second research cap correctly exposed that tail as an unknown.
The broader mean, not the attractive 100-graph mean, is used for scheduling.

An additional finite search showed that if all five 5-cycle spokes are
available, every pair of the ten local cycle/spoke edges can be avoided for
all twelve boundary orders. This did not help computationally: none of four
girth-five graphs in the 100-sample admitted the all-five-spoke reduced
condition, and the known hard fallback also rejected it. The shortcut was
removed rather than added to the production flow.

The stage-specific audit driver now supports `p3-22` and
`p2-24-exceptions`. It selects `-sat=1` only for direct order-28 girth-six
solving; on eleven matched instances this averaged 1.87 seconds, versus 2.15
seconds for `-sat=0`, while VSIDS and phase-saving overrides were worse.
Using the broader 0.178-second order-24 mean, the conservative five-stage sum
is approximately 5.1 target days. A target-server pilot is still needed, and
the flow remains above the requested two-to-three-day schedule under the
deliberately conservative 2.5x per-core scaling assumption.

A targeted solver portfolio improved the 4-cycle-free tail. The contracted
4-cycle instances retain the lean `-sat=0` formula, while direct `P_2`
fallbacks enable the sound SAT-direction strengthening constraints. On the
identified hard graph this reduced 25 seconds to 3.67 seconds. Repeating the
476-graph sample completed all 435 4-cycle contractions and 41 fallbacks under
a ten-second cap at 0.102 +/- 0.026 seconds per graph. The selected
conservative five-stage projection is therefore about 4.8 target days.

### Reverse 5-cycle census rejected

Deleting a 5-cycle from an order-26 girth-five cubic graph leaves a
21-vertex graph of girth at least five with degree sequence `2^5 3^16`.
Adding a hub adjacent to the five degree-two ports gives the degree-five
contracted graph. A complete `P_3` audit of these hub graphs would imply
`P_1` for every cyclic 5-cycle expansion.

The proposed smaller census is not actually small. Nauty command
`geng -c -t -f -d2 -D3 21 29:29 0/1000` generated 22,312 cores in one
residue, projecting to about 22.3 million. On 101 stratified hub graphs,
triple coverage processed only 49 in about 90 seconds before the experiment
was stopped, roughly 1.8 seconds per graph. This is much worse than the
0.361-second direct parent `P_1` audit, so the reverse census was rejected.

A narrower follow-up tested only the triples actually needed by 5-cycle
expansion: every pair of the five hub spokes together with every target edge,
about 340 candidate sets instead of all triples. It processed 60 of 101 hub
graphs in about 31 seconds (roughly 0.5 seconds per graph) before being
stopped. Multiplying by the 22.3-million-core census is still comparable to
or worse than direct parent `P_1` once generation is included. The specialized
prototype was removed from the production code.

## 2026-07-11: three-5-cycle contraction-core census

The most common order-28 girth-five reduction packs three 5-cycles and
contracts them to a 16-vertex graph with degree sequence `5^3 3^13`. Instead
of solving the same reduced structural class once per parent, a separate
unlabeled core census can certify every possible 5-cycle expansion.

There is a cheap exact realizability filter. Delete the three degree-five hub
vertices. The ordinary induced graph must have girth at least five. At each
hub, a cyclic order of its five neighbors is possible precisely when no two
consecutive ordinary neighbors are already adjacent; such an adjacency would
form a 4-cycle after expansion. The three order choices are independent. The
selected packing additionally requires a consecutive pair whose two edges do
not lead to other hubs. Exhaustive explicit expansion on the regression cores
agrees with this characterization.

The new `-analyze-hub-pair-sets` audit checks every triple of realizable cyclic
orders. For an uncovered triple it may choose any clean consecutive pair at
each hub and requires all six spokes in one SAT witness. The obligation is
existential in those pair choices, which is exactly what expansion needs; an
earlier prototype unnecessarily required every pair combination and had bad
tails. Ordering alternatives by their crossing distance from existing
witnesses reduced a pathological core from more than 30 seconds to 113 ms.
On 205 broader cores the corrected audit averaged 37 ms with no failure under
a two-second cap. Across 1,223 cores from five independent generator residues,
all completed after the same correction; parallel contention raised observed
means to 67--126 ms.

`scripts/generate_three_pentagon_cores.py` reproducibly generates connected
simple degree-sequence cores with nauty, filters realizable orders, hashes
every shard, and writes a generation manifest. The verification driver has a
new `n28-g5-core3` stage. The parent `n28-g5` shortcut is conditional: it is
enabled only when the driver receives a passed same-binary core manifest, and
the dependency hash is copied into the parent manifest. OOPS then limits the
packing to three 5-cycles and skips the reduced solve exactly when all three
are found.

On the established 1,054-parent girth-five sample, the conditional core stage
covered 645 graphs and the original contraction solved the remaining 409.
Mean parent time fell from 15 ms to 8 ms, with zero fallback. The attempted
complete 16-vertex count-only generation was stopped after its 12-way launch
oversubscribed the development machine; the production generator manifest,
not an extrapolated count, is authoritative.

The later complete generation finished in 1,530 seconds using 256 fixed nauty
residues. It enumerated exactly 9,366,340 connected `5^3 3^13`
degree-sequence cores and retained 2,938,761 realizable cores (63 MiB in 256
hashed graph6 shards). An independent line-count and SHA-256 pass matched the
manifest. The generation-manifest SHA-256 is
`372d7fa6e985c9011d419db70244ee714e650218dec12c482aa0854062846d45`;
the files currently reside in `/tmp/cubic28-core3-full` pending transfer to
the production server/evidence archive.

### The analogous `P_1(26)` core shortcut was rejected

Order-26 parents frequently also pack three 5-cycles: 126 of 212 matched
girth-five graphs packed three and another six packed four. Contracting three
gives a 14-vertex `5^3 3^11` core. This universe is small: nauty generated
339,516 degree-sequence cores in 45 seconds, and the realizability filter kept
exactly 135,774 (16 deterministic shards, all hashed).

The finite patch search found that a prescribed 5-cycle cycle edge `c_i` can
be protected, uniformly over all boundary rotations, only by allowing the
adjacent spoke pair `(s_i,s_{i+1})`. A prescribed spoke `s_i` permits the three
adjacent spoke pairs not containing `s_i`. A prototype therefore audited each
realizable order triple, every surviving core edge, and all fifteen local
cycle-edge targets, coordinating hub-to-hub spokes so no local patch used one
at both endpoints.

The universal claim does not currently verify. On a 90-core stratified sample,
one realizable core left the obligation `orders=(1,1,1)`, target cycle edge 2
at hub 2 unresolved. All eight legal pair combinations exhausted repeated
5,000-conflict caps; SAT-direction strengthening did not help, and a bounded
deeper run still did not produce a witness. Only about 3% of the 135,774 cores
are planar, so restricting to the trivial planar subcase would not materially
change the schedule. The prototype solver mode was removed. Direct biased
`P_1(26)` remains the rigorous selected stage.

### Crossing-pair batching did not reproduce

A final P1 search heuristic batched the two edges of an actual witness
crossing when both remained uncovered, falling back soundly to one edge when
the pair was infeasible. On the matched 212-graph sample its first run was
339 ms/graph versus an earlier 376 ms singleton baseline. An immediate paired
singleton repeat was faster at 332 ms/graph. Witness counts changed only
slightly and individual regressions were mixed, so the apparent improvement
was machine noise. The batching prototype was removed.

### Evidence hardening and current gate

The production driver now separates worker concurrency from residue parts per
input, enforces optional exact processed and reduction counts, records host
and driver hashes, and rejects duplicate input basenames. The three-5-cycle
parent shortcut requires a passed same-binary core audit whose processed count
equals its linked complete generation manifest; the linked manifest must
still exist with the recorded hash. A tiny smoke core audit therefore cannot
authorize any skipped parent solve.

`scripts/benchmark_cubic28.py` runs bounded deterministic residue pilots with
the exact production flags and log validators. It records per-task tails,
CPU-seconds per graph, parallel throughput, binary/input/driver hashes, and
host metadata. Suggested target-server moduli are in
`cubic/verification.md`.

After cleanup, `scripts/run_tests.sh` passes in full, including the Glucose
assumption regression, proper-crossing patch regression, core-filter explicit
expansion regression, core analyzer smoke, manifest rejection test, and both
verification and pilot drivers. `git diff --check` is clean, all Python tools
compile, and `main_cubic.tex` builds to a nine-page PDF (only the pre-existing
table/bibliography warnings remain).

The selected conservative target projection is now about 4.1 days: the core
shortcut halves the sampled order-28 girth-five parent cost from 15 to 8 ms,
but direct `P_1(26)` and girth-six solving remain dominant. This extrapolation
cannot settle the requested two-to-three-day reproducibility window. The next
required evidence is the matched pilot on the actual 32-core server; no host
or access method is present in the repository.

### Release/LTO build rejected for production

The default Makefile build already uses `-O3` with explicit OOPS checks
enabled. A fresh `make r` adds LTO and `-march=native`; OOPS checks remain, but
the changed binary layout changes CDCL search. On the matched 212-graph P1
sample, `oops_release` averaged 384 ms versus the latest 332 ms standard run
and produced one ten-second tail. On the 1,054-parent girth-five contraction
sample both builds rounded to 8 ms/graph. The release build therefore offers
no reproducible gain and the standard checked binary remains selected.

## 2026-07-13: simplifying the 4-cycle proof

The 4-cycle work was reconsidered from first principles, with the goal of
removing special cases even at the cost of a modest amount of additional
finite computation.

### A crossed split edge cannot be repaired by the obvious local move

The ordinary 4-cycle split replaces a 4-cycle by two adjacent cubic
vertices. Its reverse expansion is immediate when the new middle edge is
uncrossed. It was tempting to hope that a crossed middle edge could always be
repaired inside a neighborhood of that edge and its crossing.

`scripts/search_square_crossed_split.py` now checks this possibility for all
eight boundary orders of that neighborhood. It allows arbitrary matchings of
proper crossings among the four new cycle edges and the old edge that crossed
the middle edge. Repeating the check for all 16 subsets of attachment spokes
that may also receive a crossing still gives eight failures out of eight
boundary orders in every case. Replacing the first split by the other
consecutive split, while requiring its middle edge to be uncrossed, also
fails in all eight orders, even when all four attachment spokes are available.

This is a finite topological obstruction to a *local* proof. It does not
disprove the global conjecture that another drawing, or the other reduced
graph, always has an uncrossed middle edge. Complete and bounded experiments
continue to support that conjecture:

- all 7,279 order-18 core graphs containing a 4-cycle passed the two-split
  reduction (93 reduced graphs were planar), in 56 seconds;
- all 9,129 4-cycle graphs in residue `0/10` of the order-20 core passed (37
  planar), in 117 seconds; and
- all 13,377 4-cycle graphs in residue `0/100` of the order-22 core passed (27
  planar), in 277 seconds.

These results are evidence only. Because the same-drawing local statement is
false, the global conjecture is not selected as a proof obligation.

### Larger neighborhoods and separation pairs do not remove the generic case

On residue `0/1000` of the available order-26 core shard, 5,423 of 6,104
graphs contain a 4-cycle. The attachment-pattern census counted 4,882 graphs
with at least one 4-cycle having four distinct external neighbors and no edge
among those neighbors. Thus special repeated-neighbor and adjacent-neighbor
rules cannot cover about 90% of the 4-cycle graphs in this sample.

`scripts/analyze_square_neighborhoods.py` measures disjoint 4-cycles. In the
same sample, 3,761 of the 5,423 4-cycle graphs (69.4%) have two
vertex-disjoint 4-cycles, and 1,821 have three. Of the 1,662 graphs with
packing number one, 1,490 have exactly one 4-cycle. A larger-neighborhood
rule can therefore accelerate many instances, but it cannot replace the
single generic 4-cycle case.

Nauty's connectivity classification gives 2,605 connectivity-two and 94,536
connectivity-three graphs in the 97,141-graph order-20 core, and 28,274 versus
1,404,438 in the 1,432,712-graph order-22 core. Only 2.68% and 1.97% of these
cores, respectively, are not 3-connected. General decomposition at
separation pairs would introduce boundary-state bookkeeping for a class of
only about 2%; the existing planar two-terminal replacement remains useful,
but a full decomposition is rejected as a central proof device.

### Deleting all four vertices asks for a harder drawing state

Deleting the 4-cycle leaves four degree-two attachment vertices. Completing
them with the opposite matching gives an order-24 cubic graph. If the two new
matching edges cross each other, replacing their crossing by the deleted
4-cycle is a clean certificate. This saves four vertices at once, but it asks
for a *specified crossing*, rather than for specified edges to be uncrossed.

An order-24 cubic graph has 36 edges and 558 independent edge pairs. A
1-planar drawing contains at most 18 crossing pairs, while one drawing can
leave hundreds of prescribed pairs simultaneously uncrossed. Consequently
crossing-pair coverage has intrinsically poor amortization, and OOPS does not
currently expose a reviewed force-this-pair-to-cross audit. Completing the
ports with two merely uncrossed edges is insufficient: their curves need a
common-face/boundary-order condition. The delete-all-four direction is
therefore rejected for the production proof.

### A stronger single order-22 property is too slow

A five-edge flexibility statement at order 22 could absorb all local
4-cycle obligations into one finite property. A bounded run on 15 stratified
order-22 core graphs found that all 15 do satisfy the property, with zero
forced or unknown sets. It averaged 20.913 seconds per graph (311 CPU seconds
total), compared with 0.378 seconds in the current matched direct two-edge
pilot at order 24. This is approximately 55 times slower per graph and has much worse
tails, so stronger order-22 flexibility is rejected.

### Selected simple 4-cycle architecture

Only two finite flexibility statements are needed:

1. every relevant cubic graph through order 24 is two-edge flexible; and
2. every biconnected triangle-free cubic graph of order 26 with no 4-cycle is
   one-edge flexible.

For an order-28 graph containing a 4-cycle, perform one ordinary split and
call its new middle edge `e`. If the order-26 graph has no 4-cycle, statement
2 supplies a drawing with `e` uncrossed. Non-core outputs reduce to smaller
graphs by the flexibility-core lemma and are already covered by statement 1.
If the order-26 graph has another 4-cycle, split it. In the order-24 graph,
require its new middle edge and the image of `e` to be uncrossed. If `e` lies
on that 4-cycle, choose the split that makes `e` an internal cycle edge, so
only the new middle edge is required. Statement 1 therefore expands the
second 4-cycle and then the first.

This proof may use either first split but does not depend on finding a lucky
one. It removes the `P_3(22)` hierarchy, the degree-four adjacent-pair
certificate, and the 4-cycle-free order-24 exception stage from the selected
protocol. A fresh deterministic 96-graph single-core pilot with the exact
production flags averaged 0.378 seconds per graph. Repeats under varying host
load ranged up to 0.436 seconds; the prior independent 100-graph mean of 0.419
seconds remains the conservative scheduling value. Direct `P_2(24)` therefore
projects to about 115 CPU-days, versus about 43 CPU-days for the more elaborate
hierarchy. The simplicity premium is roughly 72 CPU-days, or 0.9 days under
the deliberately conservative 80-development-core-equivalent server model.
The actual 32-core server pilot remains the scheduling authority.

The selector activity boost was retested rather than specialized by stage.
On the 96-graph pair sample, no boost and the selected one-million boost were
within run-to-run noise (0.429 versus 0.436 seconds in the paired repeat). On
129 matched order-26 girth-five graphs, however, no boost averaged 0.665
seconds while the selected boost averaged 0.436 seconds. The common
one-million setting is therefore retained for both flexibility stages.

### Girth-census generation is now manifest-bound

The local Minibaum source was inspected rather than relying on an informal
command description. Its `s` option restricts output to the requested order,
`g` emits graph6, and `m RESIDUE MODULUS` partitions the canonical generation
tree. The command `minibaum5 10 3 s g` reproduced the known 19 connected cubic
graphs of order 10; the two parts contained 11 and 8 records.

`cubic/generate_minibaum_census.py` now runs these deterministic parts in
parallel, validates the order byte of every graph6 record, writes shards and
sidecars atomically, supports hash-checked resume, and records the binary,
source, command, count, and SHA-256 for every shard. The production driver
accepts `--census-manifest` for the five Minibaum-backed stages and rejects a
wrong graph class, incomplete count, altered shard, or incomplete shard set.
The exact generation commands and the distinction between the
656,783,890-record minimum-girth-five input and the 652,159,389-record
exact-girth-five audit are documented in `cubic/verification.md`.

### Contracting the whole 4-cycle: topologically clean, census rejected

The possibility of eliminating all four vertices of a 4-cycle was revisited
in its simplest form: contract the 4-cycle to one degree-four vertex.  Its
four incident edges have only three cyclic orders up to reflection.  The
finite proper-crossing search

```
python3 scripts/search_cycle_universal_spokes.py 4 ''
```

finds a local expansion for all three orders.  Order `0123` is planar;
orders `0132` and `0213` use one crossing between opposite 4-cycle edges.
No attachment spoke needs to be uncrossed.  Thus every 1-planar drawing of
the contracted graph expands, with no boundary state at all.

For an order-28 cubic parent the contracted graph has 25 vertices and degree
sequence `4^1 3^24`.  A direct finite audit would therefore need the
biconnected triangle-free nonplanar part of

```
geng -C -t -d3 -D4 25 38:38
```

This is not a reviewer-friendly replacement census.  On nauty 2.9.3, even
residue `1/1000000` with the lowest legal early split
`-X2 -x3000000` did not reach its first output or count in nine minutes on
one development core and was stopped.  The large fixed canonical-generation
cost is already much worse than generating the selected cubic girth
censuses.  Without a theorem reducing these degree-four graphs to a much
smaller already-audited class, the clean local patch merely trades the
two-edge property for a new and expensive graph class.  It is therefore not
selected.

### The order-at-most-24 census is now manifest-bound

The small flexibility audit was the last production input described only by
shell pipelines and a total line count.  The new
`scripts/generate_p2_census.py` runs the exact nauty classes used by the
proof: all connected cubic graphs through order 20, followed by the
biconnected triangle-free nonplanar cubic cores at orders 22 and 24.  The two
large orders use fixed canonical residues.  The wrapper validates every
graph6 order byte, checks every per-order count and the 25,741,152 total,
supports hash-checked restart, and writes a relocatable shard manifest with
the `geng`, `planarg`, and driver hashes.

`cubic/verify_cubic28.py` now accepts this manifest for `p2-24`, requires
the exact per-order inventory and exact shard set, and binds its hash and
generator hashes into the verification result.  The old unsharded nauty
commands remain in the documentation as a readable description of the graph
classes; the manifest-producing wrapper is the production command.

### Crossing minimization in coverage witnesses is worth keeping

The incremental coverage loop greedily removes avoidable crossings from each
checked SAT witness before using its crossed-edge set to find the next
obligation.  This is not needed for soundness, so its cost was measured on the
same deterministic 96-graph order-24 pair sample.  With minimization the run
used 45 seconds of OOPS time (470 ms/graph, 46.21 seconds wall); without it,
the less useful witness sets caused more assumption work and used 62 seconds
(655 ms/graph, 63.62 seconds wall).  Both runs completed all pair obligations.
The roughly 1.39x regression without minimization is large enough that the
existing greedy step remains selected.

### Three-5-cycle core shortcut rejected after a broader tail audit

The completed 2,938,761-core generation made it possible to test beyond the
earlier thousand-core pilots.  A deterministic quarter of its first shard
contains 1,967 realizable cores.  With the original clean-pair enumeration,
the run processed only six cores in roughly 4.5 minutes before it was stopped;
one graph/seed combination caused the tail.  A ten-second solver cap allowed
the quarter-shard to finish in 259.53 seconds (131 ms/record), but correctly
reported one incomplete hub-order obligation, so that run could not be
accepted as evidence.

Two attempts to remove the procedural tail were also bounded and rejected.
Reducing the per-alternative conflict budget fixed the first hard core but
soon found another.  Encoding the existential clean-pair choice directly as
one guarded SAT condition per order tuple made the first 900 cores fast, then
again developed a multi-minute tuple tail; caps from 5,000 down to 50
conflicts did not make the whole quarter-shard predictable.  No failed
mathematical obligation was observed, but a complete special census would
need additional portfolio and retry machinery.

The shortcut saves about 7 ms on the 61% of sampled girth-five parents that
pack three 5-cycles, roughly 49 development CPU-days over the whole parent
census.  That is only about 10--12% of the conservative end-to-end schedule,
while it adds a fifth audit, a separate non-cubic generator, a hub-order
theorem, and the least predictable SAT tails in the flow.  It is therefore
removed from the selected protocol.  Ordinary packed-5-cycle contraction
certified all 949 parents in a fresh matched sample at 15 ms/graph with zero
fallback; contracting exactly one or at most two 5-cycles was simpler but
slower at 55 and 26 ms/graph, respectively.  The selected parent reduction
continues to use the deterministic maximal packing and an ordinary checked
reduced solve.

## 2026-07-13: selected-surface audit and 6-cycle reduction

### Rejected modes removed from the reviewer-facing binary

The many exploratory reduction, retry, packed-coverage, and hub-order modes
were useful for choosing the proof but made the resulting program difficult
to audit.  They have now been removed.  The new reviewer-facing surface is
limited to two flexibility checks (`-analyze-uncrossed-edges` and
`-analyze-uncrossed-pairs`) and the selected 5-cycle and 6-cycle
contractions.  The pre-existing biconnected/cubic/triangle reductions were
restored verbatim.  The selected modes are mutually exclusive and reject
directed, IC/NIC, move-solver, UNSAT-only, transitivity-dropping, Satsuma,
BreakID, and unrelated fixed-edge combinations.

The constrained-edge path was also re-audited against the graph-specific
crossing filters.  The almost-twin crossing rule and twin-vertex ordering use
input automorphisms and can move a named edge, so they are disabled whenever
an uncrossed-edge requirement is present.  Density, adjacency, and
separating-cycle filters rule out topologically impossible crossings.  The
remaining degree-three rule has an explicit local redraw that only removes
crossings and therefore preserves every already-uncrossed edge; that proof is
now in `main_cubic.tex`.

On a 97-graph order-24 pair sample after this hardening, mean time was 291 ms
with zero forced or unknown pairs.  A fresh exact-Minibaum order-26
minimum-girth-five sample contained 130 records and averaged 283 ms, with
zero forced or unknown edges.  The rarer order-26 girth-six band averaged
697 ms over 65 sampled records with the safe degree-three filter enabled.
The documentation keeps the older conservative scheduling values of 419 ms
and 436 ms rather than extrapolating from these favorable samples.

### Repeated-assumption SAT calls now have a differential regression

The vendored Glucose fixes for repeated assumption calls and decision-heap
rebuilding are no longer covered only by a toy formula.
`cubic/test_solver_assumptions.cpp` deterministically generates 250 random
CNFs on seven variables and makes 80 alternating assumption queries on every
satisfiable formula.  Every SAT/UNSAT result is compared with exhaustive
truth-table evaluation, and every returned model is checked.  The complete
test takes about two seconds.  The solver also scans the original clauses
before accepting a SAT result under assumptions.

### 5-cycle certificate reduced to its actual twelve cases

The packed-5-cycle theorem now requires exactly two consecutive protected
attachments per selected 5-cycle.  Its local certificate was reduced from
26 alternative patches to one patch for each of the twelve boundary orders
up to reversal.  The C++ regression checks that all twelve orders occur
exactly once, that no spoke other than the protected pair is crossed, that
the listed crossings form a matching, and that the planarized patch embeds
in a disk.  A 949-parent order-28 girth-five sample still completed at about
15 ms/graph with zero fallback.

### Census manifests strengthened

Both census generators now reject malformed graph6 records by checking the
exact record length, order byte, and graph6 byte range.  Supplying a census
manifest to `cubic/verify_cubic28.py` automatically enforces the exact
processed total and rejects a conflicting manual total before any solver
work.  Verification manifests bind the complete generation-manifest hash,
the exact shard inventory, generator hashes, and (for Minibaum) the archived
source record.  The Minibaum source used in development has SHA-256
`eaed6d0a9c963a40503893f68378102169ed8e281f78d626a4e67f9b684820de`.

### Four protected spokes make a 6-cycle contraction practical

Direct order-28 minimum-girth-six solving remained the slowest per-instance
stage.  The local patch search first established the sharp behavior for the
restricted patch family: protecting zero through three spokes cannot cover
all 60 boundary orders.  Four protected spokes suffice when the two
unprotected spokes are opposite; representative missing pairs at distance
one or two still leave eight or five boundary orders uncovered.  The selected
certificate therefore contracts one chordless 6-cycle and protects
`s1,s2,s4,s5`, leaving opposite spokes `s0,s3` unprotected.

`cubic/test_cycle_expansions.cpp` now contains one explicit patch
for all 60 boundary orders.  For every row it verifies that crossings form a
matching, only the four protected spokes may cross, the alternating crossing
gadgets planarize, and the boundary wheel gives a disk embedding.  In a
cubic graph of girth at least six, a 6-cycle automatically has no chord and
six distinct outside neighbors, so no structural subcases are needed.

On the same 101-record residue sample, direct solving took 176 seconds
(1,747 ms/graph), whereas the contraction took 49 seconds (490 ms/graph), a
3.6-fold speedup with zero fallback.  A broader 1,001-record residue sample
with a ten-second cap exposed three constrained tails with the complete
encoding (`-sat=0`); two parents were then solved directly and one timed out.
Without the artificial cap, the exact same sample completed all 1,001
contractions with zero fallback or unknown in 458 ms/graph.  `-sat=1`
averaged 438 ms on the capped sample but still had two capped fallbacks.  Its
small mean advantage did not justify documenting two additional
redraw-based crossing constraints for four named protected edges, so the
selected production stage uses the simpler complete encoding, `-sat=0`.

The minimum-girth-six census has 4,624,501 graphs.  Exactly 4,624,480 contain
a 6-cycle and are subject to the contraction; the remaining 21 have girth
at least seven and are solved directly.  The verification driver now requires
those exact aggregate counters and zero contraction fallback.  The 21-record
subcount is independently reproducible by running Minibaum at minimum girth
seven.  With the conservative 32-core-times-2.5 scheduling model, this lowers
the current total estimate from about 6.2 to about 5.3 target days; the real
target-server pilot remains authoritative.

## 2026-07-14: flexibility hierarchy and a smaller order-22 audit

### The 4-cycle patch already protects three consecutive cycle edges

The four endpoint-rotation patches introduced for two adjacent edges of a
4-cycle use either no crossing or the single crossing `c3 x s1`. They
therefore protect not just `c0,c1`, but all three consecutive cycle edges
`c0,c1,c2`. The patch regression now states and checks this stronger fact;
no new local drawing is required.

This yields a useful additional reduction. To prove `P3(22)` for a graph
with a 4-cycle, split the cycle to an order-20 graph. If all three target
edges are outside the cycle, require those three images and the new middle
edge, a `P4(20)` obligation. With one or more target cycle edges the number
of requirements decreases; adjacent or three-cycle-edge cases use the
three-edge patch and one protected attachment. Consequently it suffices to
audit `P4` through order 20 and audit `P3` directly only at order 22 and
girth at least five.

The exact direct order-22 class has 90,938 records, versus 1,432,712 records
in the previous biconnected nonplanar triangle-free order-22 core. For the
`P4` base, all connected cubic graphs through order 18 plus the 97,141
biconnected nonplanar triangle-free order-20 cores give 143,123 records.
The flexibility-core lemma proves that this base covers every subcubic graph
through order 20.

Matched bounded pilots all completed with zero uncovered or unknown edge
sets:

- `P4(20)` on 501 order-20 cores: 901 +/- 70 ms/graph;
- `P4(18)` on 498 connected graphs: 170 +/- 16 ms/graph; and
- `P3(22)` on 304 minimum-girth-five graphs: 1,616 +/- 144 ms/graph.

These measurements project to about 1.1 development CPU-days for the `P4`
base and 1.7 days for the direct `P3(22)` class. Including the measured
16.8-day 4-cycle-free `P2(24)` audit, the lower-order flexibility work drops
from about 35 to about 20 CPU-days. Under the conservative target scaling,
the end-to-end projection falls from about 4.2 to about 4.0 days.

### One flexibility mode replaces four special-purpose analyzers

The implementation now exposes one reviewer-facing option,
`-edge-flexibility=k`, for `k=1,2,3,4`. It replaces separate edge, pair,
triple, and quadruple modes while using the same incremental selector and
finite-witness coverage algorithm. The solver reports one generic coverage
summary, and the verification driver checks the stage's expected `k`. Smoke
tests for all four values pass on the same 22-record reduction file.

The proof text also no longer claims that a selector is an exact encoding of
an uncrossed edge. Only the sufficient direction is needed: every accepted
SAT model is reconstructed and the required edges are checked directly. An
UNSAT or unknown selector call rejects the audit, so selector completeness is
irrelevant to the finite-witness certificate.

### Two SAT search-order experiments rejected

Combining two uncovered obligations into one larger assumption reduced the
number of witnesses but did not improve elapsed time. On triples, a matched
20-graph diagnostic reduced roughly 44 witnesses to 31 with essentially
unchanged runtime; the broader triple result was only about 3% faster, while
pair batching regressed from 1,014 to 1,206 ms/graph. The prototype was
removed.

A second prototype preferred edges that appeared least often in stored
crossing sets when constructing the next hitting set. On the same rebuilt
binary it regressed `P3(22)` from 705 to 774 ms/graph (10%) and `P2(24)` from
859 to 1,122 ms/graph (31%). This option and all supporting code were also
removed. The selected coverage loop keeps the simpler deterministic edge
order.

### Final bounded P1 search experiments

Several additional changes were tested on the dominant `P1(26)` stage and
rejected:

- Re-boosting every selector for an edge still crossed in all witnesses
  regressed the matched 101-graph girth-six sample from 1,216 to 1,400
  ms/graph and increased ten-second coverage tails from two to three.
- The full optional `-sat=1` encoding regressed the same sample to 1,496
  ms/graph with four tails. The selected production stage therefore keeps
  the smaller complete stack encoding, `-sat=0`.
- Trying four deterministic greedy orders when removing redundant witness
  crossings improved one hard sample by about 10%, but was exactly neutral
  on the representative 98-graph girth-five sample (399 versus 398
  ms/graph). Since girth five dominates the production census, the original
  one-pass cleanup remains selected.
- Choosing the still-uncovered edge with the fewest possible crossing
  partners was neutral on both matched samples (397 versus 398 ms/graph on
  girth five, and 1,076 versus 1,078 ms/graph on girth six). Inspection of
  the progress trace indicates that it usually selected the same edge as the
  simpler deterministic order. The specialization was removed.

### 5-cycle-to-path reduction does not eliminate P1 locally

A 5-cycle can be replaced combinatorially by a three-vertex path carrying
two, one, and two attachment edges. If both internal path edges are
uncrossed, the established `P2(24)` property is enough to restore the
5-cycle for an unconstrained drawing. To preserve an arbitrary outside edge
with only `P2`, however, one internal path edge would have to be allowed to
cross.

`scripts/search_pentagon_path_split.py` enumerates the eight boundary orders
of a neighborhood in which one path edge is uncrossed and the other crosses
an old edge. It allows arbitrary proper-crossing matchings among every new
5-cycle edge and the old crossing edge, while keeping attachment edges
unchanged. All eight orders fail. Thus this natural same-drawing local move
cannot replace the direct `P1(26)` audit.

### Production generator exercised on the real order-22 class

The documented Minibaum sharding command was first checked at order 10,
where its four shards contain exactly the Petersen graph. It was then run on
the complete order-22 minimum-girth-five class with 14 canonical shards. The
wrapper regenerated exactly 90,938 records in 0.65 seconds, validated every
graph6 record, and recorded the expected source hash
`eaed6d0a9c963a40503893f68378102169ed8e281f78d626a4e67f9b684820de`.
The verification driver now re-hashes every input and the census manifest at
the end of a run as well as at the beginning.

### Complete small-stage generation and sealed evidence files

The production `P4(20)` generator was run to completion with nauty 2.9.3.
It produced 143,123 records: 1, 2, 5, 19, 85, 509, 4,060, and 41,301
connected cubic graphs at orders 4 through 18, plus 97,141 biconnected,
nonplanar, triangle-free cubic cores at order 20. The fourteen order-20
residue classes and all smaller files were independently counted as graph6
records. The first generation took 22.6 seconds; a resume pass revalidated
the existing shards and manifest without regenerating them.

The complete order-22 minimum-girth-five Minibaum generation was also
repeated. Its fourteen canonical residue classes contain exactly 90,938
records. Thus both of the two smallest production input classes have now
been exercised end to end with the actual generators, rather than only with
test fixtures.

The three evidence scripts were hardened against files changing while a run
is in progress. The nauty and Minibaum wrappers hash their generator binary,
source or auxiliary binary, driver, and every output shard at both the
appropriate start and the end of generation. The verification driver hashes
the OOPS binary, all graph inputs, its own source, and each part log; a
resumed part is accepted only if the stored log hash still matches. The final
manifest includes the complete part-result inventory. Manifests are written
atomically, so an interrupted write cannot look like completed evidence.

### Local crossing checker strengthened again

Every displayed crossing in the local patch tests is now represented by a
4-cycle whose four vertices receive the edge ends in alternating order,
plus an apex adjacent to all four cycle vertices. The surrounding disk is
similarly fixed by a boundary cycle and an outer apex. The two wheels prevent
a planarity test from satisfying the gadget by reversing or folding an
unintended face. The 4-cycle, 5-cycle, and 6-cycle patch tables still pass.
The negative 5-cycle-to-path search uses the same stronger crossing gadget.

An adversarial review made the wheel argument fully executable rather than
implicit. For each displayed crossing, the checker now removes that wheel's
apex and rim edges while retaining all other planarized crossings, and
verifies that the entire exterior is one connected bridge incident with all
four rim vertices. The exterior bridge and the apex must therefore be on
opposite sides of the rim, leaving the apex side empty for the intended
proper crossing. It verifies the same connected-bridge condition between the
local patch and the outer boundary wheel. All 4-cycle, 5-cycle, and 6-cycle
rows still pass these stronger conditions.

The same review found and repaired one genuine omission in the written
order-22 reduction. With two adjacent prescribed 4-cycle edges, the original
patch used the attachment at their common vertex as its possible local
crossing. That is insufficient when this attachment is itself the third
prescribed edge. A second four-row table instead uses the next attachment:
the two nontrivial rotations use `c3 x s2`, preserving `c0`, `c1`, and `s1`.
The order-20 obligation is then `{middle edge, s1, s2}`, still within
`P4(20)`. The strengthened crossing-wheel checker verifies all four new rows.
Thus the six-stage architecture is unchanged, but the proof now covers the
previously missing prescribed-edge placement explicitly.

The ordinary 4-cycle split was also promoted from a prose drawing argument
to a third checked four-row 4-cycle table. Its two nontrivial rotations use
only `c1 x c3`; the two cycle edges created in the endpoint disks remain
uncrossed. The same proper-crossing, connected-bridge, and outer-boundary
checks pass, so every 4-cycle patch invoked in the proof now has an executable
local certificate.

The final proof overview was also simplified by removing the false-twin and
ladder reductions from the selected argument. They were inherited from
earlier filtering experiments, but the production census does not omit graphs
on either basis: every biconnected cubic triangle-free core is audited. The
only core reductions now cited are blocks, low-degree suppression, and the
triangle/diamond two-terminal argument actually used by the theorem.

The verification driver now independently enforces the exact vertex count,
edge count, and both minimum and maximum degree three for every fixed-order
stage. A small general `-min-degree` input filter was added to OOPS so the
mixed-order `P4(20)` stage can make the same direct cubicity check. Any
malformed or wrong-order record is filtered out and causes the exact
processed-count assertion to fail. New order-22,
order-26, and order-28 girth-five fixtures ensure the smoke tests exercise the
same structural filters as their production stages.

Minibaum provenance is now sealed at shard granularity. The source file is a
required generator argument; every shard sidecar records both its source hash
and executable hash, and resume accepts the shard only when both match. The
verification driver rejects a Minibaum census manifest without a valid source
hash. This prevents a changed source file from being associated with old
resumed shards under a newly written manifest.

### Separation pairs are too rare to simplify the main audits

Exact nauty connectivity counts show that decomposition at cut vertices or
separation pairs would affect almost none of the direct girth-five classes.
At order 22 the counts by vertex connectivity one, two, and three are 1, 10,
and 90,927. At order 24 they are 8, 152, and 1,620,319. In an order-26 sample
of 15,581 records, every graph was 3-connected. A special theorem and code
path for the exceptional fraction would add proof cases without materially
reducing runtime, so this direction was rejected.

### Updated order-26 timing and assumption-solver audit

The repeated-assumption repair in the vendored Glucose solver was kept: it
processes assumptions at their own decision levels, disables two incompatible
chronological-backtracking shortcuts while assumptions are active, and
rebuilds all branching heaps before a new solve. A proposed diagnostic scan
of every clause after a SAT result changed a matched 392-graph `P1(26)` run
from 0.28631 to 0.28570 seconds per graph, a neutral 0.2 percent difference,
so the scan was removed and the smaller solver patch retained.

A broader independent `P1(26)` pilot processed 974 exact-Minibaum records in
299 seconds, or 0.307 seconds per graph, with zero uncovered or unknown edge
requirements. The production estimate now budgets 0.350 seconds per graph.
That is about 128 development CPU-days for all 31,478,584 order-26 records.
Together with the selected lower-order and order-28 stages, the current
32-core projection is about 3.6 target days under the deliberately
conservative 2.5-fold per-core speed ratio. A matched target-server pilot is
still the authoritative estimate.

### A stronger path reduction is mathematically promising but slower

For each 5-cycle, `scripts/analyze_pentagon_path_reductions.py` replaces its
five vertices by a three-vertex path with attachment multiplicities two,
one, and two. Every one of 15,581 records in one order-26 sample had an
orientation whose reduced graph still had girth at least five. In a separate
1/1000 residue of 35,522 records there was one strict exception; in a
137/1000 residue of 1,445 records there were 31. All of those exceptions were
still reducible after allowing a new path edge to lie on a 4-cycle and
then applying the already-proved 4-cycle split. A 151,522-record 999/1000
residue had no strict exception.

The local patch search found the exact extra flexibility that would make this
route rigorous. It is enough for the reduced order-24 graph to have a drawing
in which all three edges incident with the middle path vertex and any one
additional specified edge are uncrossed. All four boundary orders (eight
before reflection) can then be expanded back to the 5-cycle; either outside
cycle edge adjacent to an endpoint can also remain uncrossed uniformly.

This stronger property was prototyped as a constrained SAT coverage mode and
measured on the same 203-record order-24 residue used for direct triple
coverage. It certified only 123 records in 8.5 minutes and did not finish
within the ten-minute experimental budget. The four-edge selector queries
create the same long tails as unrestricted `P3(24)`, so the new theorem does
not pay for its extra proof and implementation case. The prototype solver
mode was removed. The two analysis scripts remain as a documented rejected
direction; the selected production flow continues to use direct `P1(26)`.

## 2026-07-14: target-preserving 5-cycle contraction for P1

The rejected 5-cycle-to-path direction suggested a different use of a
5-cycle: contract the whole cycle and tailor the local expansion to the one
parent edge that must remain uncrossed. Two new finite tables close all twelve
boundary orders:

- if the target is 5-cycle edge `c0`, protecting the endpoint attachments
  `s0,s1` permits a patch that never crosses `c0`; and
- if the target is attachment `s0`, protecting `s0,s1,s2` permits a patch
  that never crosses `s0`.

The ordinary two-spoke 5-cycle table already handles a target edge outside
the 5-cycle. Therefore, after contracting one 5-cycle, every one of the 39
parent-edge obligations becomes one explicit requirement of size two or
three on the same order-22 reduced graph. The new
`-contract-pentagon-flexibility` mode shares one SAT encoding and its checked
witnesses among all 39 requirements. The strengthened crossing-wheel and
connected-bridge regression accepts both new twelve-row tables and explicitly
checks that the target edge is absent from every listed crossing.

On the same 312 exact-Minibaum records, direct `P1(26)` took 0.294 seconds per
graph and the contraction took 0.162 seconds, a 1.8-fold speedup. Two
independent 974-record residues completed with zero uncovered or unknown
requirements at 0.181 and 0.178 seconds per graph. The production budget is
0.220 seconds for the 31,297,357 order-26 records that contain a 5-cycle.

The minimum-girth-five census also contains 181,227 graphs of minimum girth
at least six, which have no 5-cycle. The mode fails closed to the ordinary
one-edge audit on exactly this direct tail. A 182-record sample completed with
zero failures at 1.605 seconds per graph, about 3.4 development CPU-days when
extrapolated to the complete tail. The entire order-26 stage now budgets about
83 CPU-days rather than 128, lowering the conservative 32-core target
projection from about 3.6 to about 3.0 days.

The production driver requires zero specialized fallback and, when linked to
the complete Minibaum manifest, exactly 181,227 no-5-cycle records. The P4
and Minibaum manifest checks were also strengthened to recompute the complete
canonical residue inventory from the shard entries. For P4 the verifier now
recomputes the per-order totals from the individual shards instead of trusting
only the top-level `records_by_order` field.

An offline evidence auditor now closes the archival loop. Given the six final
stage manifests, `cubic/audit_cubic28_evidence.py` re-hashes the common OOPS
binary, every graph input, census manifest, and part log; reparses every final
counter; reconstructs the complete input-by-residue task inventory; and
enforces all production totals and the two small structural subcounts. A
partial-mode regression applies the same checks to the six smoke manifests.

The complete project regression suite passed after this integration,
including both order-26 branches, all 4-cycle/5-cycle/6-cycle local tables,
the repeated-assumption differential test, all six stage drivers, wrong-order
input rejection, and the new offline six-manifest audit. All Python tools
compile, `git diff --check` is clean, and `main_cubic.tex` builds to a ten-page
PDF. The remaining operational gate is the matched pilot on the actual
32-core target server, followed by the production runs themselves.

## 2026-07-14: rejected multi-5-cycle P1 contraction

Contracting several pairwise separated 5-cycles at once is mathematically
valid: at each additional 5-cycle, protect the ordinary consecutive spoke
pair, and use the target-specific pair or triple at the 5-cycle containing
the requested parent edge. This can reduce an order-26 graph as far as order
18, but it increases every simultaneous requirement from two or three edges
to roughly twice the number of contracted 5-cycles.

A prototype reduced the one-record fixture to order 18 and completed in 31
ms, versus 61 ms for the selected one-5-cycle mode. The larger matched sample
showed the opposite behavior decisively: after about six minutes it had
processed only 246 of the 312 records, with an individual record running for
more than two minutes, so the experiment was stopped. The selected
single-5-cycle mode then completed all 312 records with no failures in 58.8
seconds (0.188 seconds per graph), still within the 0.220-second production
budget. The multi-5-cycle prototype was removed. Keeping one contraction
also leaves the proof and implementation simpler and guarantees that every
reduced requirement has size only two or three.

A final proof-neutral selection heuristic was also rejected. Instead of the
first clean 5-cycle in graph order, it chose the 5-cycle whose contraction
created the most triangles in the reduced graph. On the same 312 records it
took 66.2 seconds (0.211 seconds per graph), slower than the restored
first-5-cycle run at 58.8 seconds (0.188 seconds per graph). The selected
source and binary were restored exactly after the trial.

## 2026-07-15: 90-percent simplification experiments

Two bounded experiments measured how much of the proof hierarchy can be
deleted while retaining nearly all end-to-end performance.

First, direct three-edge flexibility was run on residue `137/5000` of the
complete 1,432,712-record order-22 biconnected nonplanar triangle-free core
census. All 287 records completed with zero uncovered or unknown sets in
181.88 wall seconds, averaging 0.632 seconds per graph. Extrapolation gives
10.5 CPU-days for order 22; the smaller core orders add well under one day.
Replacing the current 2.8-day `P4(20)` plus girth-five `P3(22)` base by direct
`P3` on every relevant core through order 22 therefore costs only about
8--9 CPU-days, around 3.5 percent of the full 241-CPU-day schedule. This
supports deleting the `P4` audit and the detailed order-22 4-cycle reduction.

Second, a temporary order-26 hybrid used the ordinary 5-cycle table for
nonlocal targets and audited either five cycle edges or all ten local edges
directly on the parent graph. On the same 312 records, the selected two-table
contraction took 48.82 seconds (0.155 seconds per graph). Five singleton
local requirements took 100.77 seconds (0.322 seconds per graph), and ten
singleton requirements took 102.72 seconds (0.328 seconds per graph). A
single simultaneous 5-cycle-edge requirement certified a 156-record
sample without failures but still averaged 0.313 seconds; requiring all ten
local edges simultaneously was already UNSAT on the one-record fixture after
about five seconds. The viable hybrids would add roughly 57--61 CPU-days at
order 26 and retain only about 80 percent of end-to-end throughput, so they
are rejected. The experimental mode was removed, and rebuilding reproduced
the previously validated binary hash exactly.

The hybrid analysis nevertheless exposes one proof-only simplification. For
an attachment target `s_i`, rotate the ordinary 5-cycle table so its
protected pair is `s_{i+1},s_{i+2}` and also require `s_i` uncrossed in the
contracted drawing. The ordinary table crosses no spoke outside its protected
pair, so it automatically preserves `s_i`. This is exactly the existing
three-edge SAT requirement and needs no new solve. Thus the dedicated
`pentagonsPreserveS0` twelve-row table is redundant and can be removed with
zero performance cost; only the cycle-edge-preserving table remains special.

## 2026-07-15: five-audit proof adopted throughout the project

The simplification measured above is now the selected verification flow, not
just a proposed alternative. The proof, implementation, production commands,
evidence auditor, smoke tests, and reviewer documentation all use the same
five audits:

1. three-edge flexibility on every connected cubic graph through order 18
   and every biconnected, nonplanar, triangle-free cubic core at orders 20
   and 22 (1,575,835 records total);
2. two-edge flexibility on the 1,620,479 order-24 minimum-girth-five graphs;
3. one-edge flexibility on the 31,478,584 order-26 minimum-girth-five graphs,
   using the one-5-cycle reduction and its direct no-5-cycle tail;
4. packed-5-cycle contraction on the 652,159,389 order-28 graphs of girth
   exactly five; and
5. 6-cycle contraction on the 4,624,501 order-28 graphs of minimum girth six.

The old `P4(20)` stage and the separate order-22 girth-five stage have been
deleted from the production driver and evidence protocol. Their census
wrapper was replaced by `cubic/generate_p3_census.py`, which generates all
connected cubic graphs through order 18 plus the exact order-20 and order-22
flexibility cores. `cubic/verify_cubic28.py` now exposes `p3-22` as the
single first stage, and `cubic/audit_cubic28_evidence.py` requires exactly
the resulting five manifests. The expected first-stage count and every
per-order shard count are checked from the generation manifest.

The theorem draft now derives two-edge flexibility through order 24 directly
from three-edge flexibility through order 22. It uses one ordinary checked
4-cycle table and one checked table for the sole exceptional placement,
two adjacent prescribed cycle edges. The entire former order-22 4-cycle
case analysis and its second adjacent-pair table are gone.

The dedicated attachment-target 5-cycle table has also been removed from
the C++ certificate regression. For target attachment `s_i`, the proof
rotates the ordinary 5-cycle table so that its protected pair is
`s_{i+1},s_{i+2}` and the target is the unprotected spoke `s_4`. The
regression now explicitly forbids `s_4` in every crossing of all twelve
ordinary rows. The reduced SAT obligations are unchanged, so this deletion
has no runtime cost. Only the cycle-edge-preserving 5-cycle table remains in
addition to the ordinary table.

Validation after the rewrite succeeded: the standard C++ build completed,
all Python production tools compiled, the strengthened local patch checker
accepted the remaining 4-cycle, 5-cycle, and 6-cycle tables, the full
`scripts/run_tests.sh` suite ended with `ALL TESTS PASSED`, `git diff --check`
was clean, and `main_cubic.tex` built successfully to a ten-page PDF. The
tested OOPS binary has SHA-256
`e41b3b24efe82ebc1619ca3b84753fa3a861b73cf466c57d55ad3543194735d7`.
The full production censuses remain to be run; this entry records the tested
verification architecture and not a claim that the final evidence already
exists.

## 2026-07-15: final diff minimization

The complete working tree was reviewed path by path with the goal of keeping
the selected proof and at least 90 percent of its measured speedup while
removing research scaffolding.  The changed/untracked set fell from 34 paths
to 19.  Removed paths included an accidental newline-only change in
`src/brute_force.cpp`, six rejected search and neighborhood-analysis programs,
a redundant 5-cycle fixture, and the optional
`scripts/benchmark_cubic28.py` pilot wrapper.  The latter affected neither
verification nor runtime; production timing is already recorded by the
verification driver.  `.gitignore` retains focused rules for the C++ release
and debug builds, Python bytecode, and `main_cubic.tex` build products so
ordinary builds do not dirty the working tree.

The seven initially separate graph6 regression fixtures were also reduced to
one `data/cubic_verification_smoke.g6` file.  Its 27 records contain the same
22 small cubic coverage cases and the order-22, two order-26, and two order-28
boundary cases.  Stage filters select the relevant records.  The existing
`data/test3.s6` supplies the order-24 boundary case.  Thus the consolidation
changes no production input or proof and retains the same smoke-test branch
coverage.

Within retained files, the cleanup removed diagnostic-only edge-fixing CLI
options, encoder and solver phase timers, a packed-5-cycle histogram, an
unused witness vector, a redundant minimum-degree input option, the unsafe
unhashed-driver mode, redundant direct tests of stages already exercised by
the production driver, and verification flags whose early-exit behavior is
disabled by these certificate modes.  The user's article preamble and author
changes in `main_cubic.tex` are retained.  Historical entries above still
name some deleted exploratory scripts because they document how rejected
ideas were evaluated; they are not part of the selected flow.

All performance-bearing changes remain: canonical transitivity clauses,
optional clause-deduplication removal, incremental flexibility coverage with
SAT phase bias, correct repeated solving under assumptions, and the selected
5-cycle and 6-cycle contractions.  The cleanup itself does not alter any of
their hot paths or proof obligations.

After minimization, `make -j` succeeded and the full
`scripts/run_tests.sh` suite ended with `ALL TESTS PASSED`.  The four Python
production tools compile, all five stage-driver smoke tests and the combined
evidence audit pass, the local finite-table checker and repeated-assumption
solver regression pass, and `main_cubic.tex` builds to a ten-page PDF using
`../tex/cubic/cubic.bib`.  The tested OOPS binary has SHA-256
`c035fc868f2e45e992e28f58d2c496a5ec6d72b8288ce100ae2e9e0dce3d4a2d`.

## 2026-07-15: matched incremental-assumption measurement

A temporary comparison build measured the current incremental coverage loop
against fresh SAT calls. The fresh mode used the same selector clauses and
the same witness/hitting-set algorithm, but rebuilt the complete SAT model for
each uncovered edge requirement and asserted that query's selectors as
ordinary unit clauses. The comparison code was not retained. Both modes used
`-edge-flexibility=1 -sat=0 -Cno-dedup`, without external user propagation,
UNSAT strengthening, or external symmetry breaking.

On two deterministic minimum-girth-five residues of
`cub26_Ctd3D3_nonplanar_24_100.g6`, both modes certified every graph with zero
uncovered or unknown requirements:

| residue | graphs | assumptions | fresh SAT | speedup |
|---|---:|---:|---:|---:|
| `0/5000` | 133 | 53.72 s | 241.09 s | 4.49x |
| `1/10000` | 57 | 48.62 s | 116.37 s | 2.39x |
| combined | 190 | 102.34 s | 357.46 s | 3.49x |

The 22-record order-12 pair-flexibility smoke set independently took 0.17 s
with assumptions and 0.57 s with fresh solvers, a 3.4x speedup. Thus solver
reuse is a material performance optimization in addition to requiring the
assumption-handling correctness repair.

## 2026-07-15: Minibaum versus geng generation

The two local generators were compared on the identical order-24 class of
connected cubic graphs with minimum girth five. Minibaum used
`minibaum5 24 5 s g m RESIDUE 256`; nauty 2.9.3 used
`geng -c -t -f -d3 -D3 24 36:36 RESIDUE/256`. Output was sent to `/dev/null`
so the measurement covers generation but not file-system throughput.

Five spread-out 256-way Minibaum shards took 0.24--0.28 seconds, mean 0.26;
the corresponding five `geng` shards took 14.99--23.22 seconds, mean 18.82.
That is a 72-fold mean shard-time advantage and projects to approximately
1.1 versus 80 CPU-minutes for all 256 shards. A complete unsharded Minibaum
run produced the expected 1,620,479 records in 45.78 seconds. A 1000-way
sample made `geng` worse because repeated split-tree overhead projected to
111 CPU-minutes. Thus Minibaum is retained for the large girth-restricted
censuses; `geng` remains suitable for the much smaller lower-order cores and
as an independent count cross-check.

## 2026-07-15: girth filtering removed from OOPS

The workflow no longer adds `-min-girth` or `-max-girth` to the core solver.
The order-24 and order-26 class guarantees come from their hashed Minibaum
generation manifests. At order 28, one run now consumes the complete
656,783,890-record minimum-girth-five census: it tries packed-5-cycle
contraction first, then 6-cycle contraction, and directly solves the 21
graphs with neither cycle. This removes the overlapping order-28
minimum-girth-six input, one production stage, and all generic girth-filter
code from `src/main.cpp`.

The proof-specific orchestration was subsequently moved to
`src/verify_cubic.cpp`. OOPS now exposes one dispatcher option,
`-verify-cubic=p3|p2|p1|n28`; ordinary processing in `testOnePlanar()` contains
no flexibility or cycle-contraction branches.

## 2026-07-15: matched clause-deduplication benchmark

Clause deduplication was remeasured on one pinned core of the local Intel Core
Ultra 5 125U. Every configuration was run twice in the alternating order
dedup/on, dedup/off, dedup/off, dedup/on. All runs used the same debug-checked
optimized binary, `-sat=0`, and identical graph residues. Every graph was
certified with zero fallback, unknown, uncovered, or non-1-planar results.

| stage and sample | graphs | global dedup | no dedup | speedup |
|---|---:|---:|---:|---:|
| order-26 5-cycle flexibility, `cub26_sample.g6`, `0/100` | 383 | 56.69 s | 53.36 s | 1.062x |
| order-28 packed 5-cycles, `cub28_girth5_sample.g6`, `0/5` | 1,897 | 36.59 s | 30.67 s | 1.193x |
| order-28 6-cycles, `cub28-gir6.g6`, `0/50000` | 93 | 30.75 s | 30.49 s | 1.008x |

For order 26, disabling deduplication reduced summed encoding time from 5.71
to 1.98 seconds, increased the mean solver clause count from 63,554 to 64,048
(0.78%), and reduced peak RSS from 44.8 to 28.3 MB. For order-28 5-cycles,
encoding fell from 11.97 to 5.36 seconds, clauses rose from 36,183 to 36,410
(0.63%), and peak RSS fell from 41.7 to 19.2 MB. For order-28 6-cycles,
encoding fell from 1.53 to 0.57 seconds and clauses rose from 71,293 to 71,730
(0.61%), but the longer incremental SAT work absorbed nearly all of the
wall-time saving. Initial-solve conflict, decision, and propagation counts
were identical between configurations on every sample.

The existing deduplication implementation already uses `tsl::robin_set`. A
temporary build replacing it with whole-vector `sort_unique` was also run
twice on the 1,897-graph order-28 5-cycle sample. It averaged 36.06 seconds
and 11.25 seconds of encoding, slightly better than the robin-set results of
36.59 and 11.97 seconds, but still much slower than no deduplication at 30.67
and 5.36 seconds. Peak RSS was 18.8 MB for sort/unique versus 41.7 MB for the
robin set because the latter temporarily stores copied clause vectors.

Conclusion: `tsl` does not rescue this pass. Skipping global deduplication is
a repeatable 16% wall-time reduction on the dominant order-28 5-cycle stage,
a 6% reduction at order 26, and neutral on the smaller 6-cycle stage. The
production `-Cno-dedup` setting is therefore retained.

## 2026-07-15: decomposition of the optional SAT encoding

The four components enabled together by `-sat=1` were temporarily made
independently switchable: the exact `cross2` definitions, the clauses saying
that an edge crosses at most one other edge, the K4 crossing exclusions, and
the IC-2 crossing exclusions. The experimental switches were removed after
measurement. All runs used one pinned core, `-Cno-dedup`, and identical graph
residues. Every configuration produced the same successful verification
counters, with zero fallback, unknown, uncovered, or non-1-planar results.

| encoding | n=26 5-cycle flexibility, 383 graphs | n=28 5-cycle contraction, 1,897 graphs |
|---|---:|---:|
| `-sat=0` | 55.07 s | 30.63 s |
| `cross2` definitions only | 61.78 s | 34.11 s |
| `cross2` + at-most-one crossing | 55.47 s | 35.33 s |
| `cross2` + K4 exclusions | 61.31 s | 38.92 s |
| `cross2` + IC-2 exclusions | 72.04 s | 41.80 s |
| full `-sat=1` | 69.69 s | 39.21 s |

On order 26, the at-most-one-crossing clauses recover the approximately 11%
cost of introducing `cross2`, but do not improve on `-sat=0`. K4 is neutral
relative to `cross2` alone and IC-2 is harmful. On the dominant order-28
stage, `cross2` alone costs 11%, and every additional family makes the result
slower. No component provides a 10% speedup over the smaller encoding, so the
production verification retains `-sat=0`.

## 2026-07-15: higher threshold for OOPS-specific optimizations

The final review adopted a stricter simplicity rule: verification-specific
performance machinery in the core solver must provide roughly a 30--40%
speedup unless it is a broadly applicable formula simplification. Under this
rule the optional clause-deduplication bypass was removed from `Params`, the
verification driver, and the tests. Its measured gain ranged from neutral to
19%, so it did not justify another OOPS mode. The canonical transitivity
rewrite remains: it emits the same logical clauses without generating the
three cyclic copies in the first place and is not verification-specific.

## 2026-07-16: SAT-component matrix on general and cubic samples

The earlier SAT-component decomposition was repeated on broader fixed samples.
All runs used the normal clause-deduplication pass and one process pinned to one
core. Each graph had a 60-second timeout. The temporary component switches were
removed after the experiment. Times below are user CPU seconds; this also makes
the `sat/cross2` row comparable despite the laptop suspending during that run.

The columns mean: no optional SAT encoding; exact `cross2` definitions only;
`cross2` plus the at-most-one-crossing clauses; `cross2` plus K4 exclusions;
`cross2` plus IC-2 exclusions; and the normal full `-sat=1` encoding.

| sample | graphs | `sat=0` | `cross2` | +at-most-one | +K4 | +IC-2 | full `sat=1` |
|---|---:|---:|---:|---:|---:|---:|---:|
| complete `combo_sat.cfg` | 29,826 | 256.7 | 371.1 | 353.9 | 377.9 | 346.4 | 380.0 |
| fixed random `list1_sat.cfg` sample | 1,000 | 260.8 | 246.7 | 292.9 | 190.4 | 173.1 | 233.7 |
| fixed random `sat.cfg` sample | 75 | 354.6 | 412.1 | 306.0 | 382.8 | 354.1 | 362.9 |
| fixed random order-26 verification sample | 2,000 | 362.3 | 354.0 | 352.2 | 347.4 | 396.4 | 371.1 |
| mixed order-28 verification sample | 9,733 | 310.8 | 318.2 | 312.4 | 367.6 | 322.0 | 357.6 |

All cubic configurations certified every graph, with zero fallback, uncovered
requirement, or unknown result. The general samples contain hard instances
that reached the fixed timeout. The unknown counts by row were respectively:
`combo` 1,1,1,2,0,2; `list1` 0,0,1,0,0,1; and `sat` 1,2,0,2,2,2.

There is no universally best optional family. IC-2 is 34% faster than
`sat=0` on `list1`, but 9% slower on order 26. At-most-one is 14% faster on
`sat`, but 12% slower on `list1`. K4 is 27% faster on `list1`, but 18% slower
on order 28. Full `sat=1` ranges from 10% faster on `list1` to 48% slower on
`combo`. Consequently the experiment gives no basis for deleting any family
from general OOPS or changing its global default. For cubic verification,
`sat=0` remains the simplest choice: no SAT component improves order 26 by
more than 4%, and every nontrivial choice is neutral or slower on order 28.

## 2026-07-16: shared SAT setup for cubic verification

The cubic verifier no longer duplicates the core solver initialization and
first-solve sequence.  Generic `initSATSolver` and `solveSATModel` helpers are
now used by both `runSolver` and the incremental cubic path; only the repeated
coverage assumptions and witness checks remain in `verify_cubic.cpp`.

The two label-preservation arguments were replaced by one generic distinction:
`Params::hasCross1Restrictions()` is true for either command-line
`-fix-cross1` units or incremental cross1 assumptions.  It guards the
almost-twin crossing reduction, degree-two lexicographic rule, and vertex-twin
ordering.  Cubic verification enables incremental assumptions without adding
cubic-specific fields to the core solver.

The refactor passed normal stack and move solves, the intended
`-fix-cross1`/UNSAT mode, DIMACS generation, the exhaustive repeated-assumption
regression, the order-at-most-12 two-edge-flexibility smoke test, and the
order-28 5-cycle/6-cycle contraction smoke test.

## 2026-07-16: strict minimal cubic verifier

The four configurable `-verify-cubic=p3|p2|p1|n28` modes were replaced by one
strict `-verify-cubic` mode.  Graph order now selects the unique proof
obligation.  The verifier constructs its fixed plain-stack parameters directly
instead of parsing general OOPS SAT options, validates connected cubic input
(and girth at least five from order 24 onward), and aborts on the first missing
drawing or flexibility witness.

The following unused compatibility machinery was deleted: `-skip-planar` and
`-skip-reducible` handling, generic `isOnePlanar` fallback, retries after a
failed contraction certificate, non-1-planar/unknown/skipped counters, timeout
support, configurable SAT modes, and per-graph timing statistics.  The C++ file
fell from 523 to 333 lines.  The Python evidence driver retains the four census
jobs, but all invoke the same strict binary mode and validate theorem-oriented
success counters.

Witness crossing minimization was also removed after three runs of the
order-24 pair test changed mean time only from 6.19 to 6.34 seconds.  The
nine-line cross1 phase/activity bias was retained: removing it increased the
same mean from 6.19 to 6.95 seconds, a 12% slowdown.  All four smoke stages and
the independent evidence auditor pass with the reduced implementation.

## 2026-07-16: proof-specific contractions localized

The packed-5-cycle, target-preserving 5-cycle, and 6-cycle contraction
routines are used only by the cubic verification proof.  Their declarations
and implementations were moved into the anonymous namespace in
`src/verify_cubic.cpp`; `graph_algorithms.h` and `reducible_subgraphs.cpp` are
again unchanged from upstream.  This is a source-only relocation: the order-26
5-cycle-flexibility and order-28 5-cycle/6-cycle smoke cases still pass.

## 2026-07-17: rebased onto upstream micro-optimizations

The local planning commit and the complete uncommitted verification tree were
rebased onto upstream commit `3efec3c` (`micro-perf optimizations`).  The
overlapping SAT files merged without conflicts; both upstream's propagation,
clause-storage, and cross-variable improvements and the verification's repeated
assumption support remain present.  The build, five standard non-1-planarity
smoke instances, the order-26 5-cycle case, and both order-28 contraction cases
pass after the rebase.

## 2026-07-17: crossing variables localized in the core encoder

The cubic verifier now sets `Params::cubicVerification`; the mode is checked
to require `sat=0` and `unsat=0`.  `encodeStackPlanar` handles the mode by
introducing only the crossing variables and implication clauses previously
added by the verifier callback.  Protected-edge and flexibility queries use
negative assumptions on those variables.  The full `cross2` and standard
`cross1` constraint families remain disabled.

The verifier-specific partial crossing encoding, stack-encoding callback, and
`forbidEdgeCrossings` API were deleted.  All four proof-stage smoke tests and
the repeated-assumption regression pass with the standard encoding.

## 2026-07-17: review cleanup

Whitespace-only changes in tracked source and test files were restored to the
upstream formatting.  The theorem draft was also made implementation-neutral:
it now describes finite witness coverage and the local expansion certificates
mathematically, without discussing OOPS internals, options, source files, or
evidence manifests.  The C++ build and LaTeX compilation both pass.

## 2026-07-17: cross1 branching bias removed

The verifier's cross1 search heuristic was tested on the same fixed general
SAT samples used for the encoding review.  Each matched run used one pinned
core, `sat=1`, `unsat=0`, a 60-second per-graph timeout, and identical input;
the experimental run alone disabled phase saving, preferred `cross1=false`,
and boosted all cross1 activities by `1e6`.

| sample | baseline | biased | change |
|---|---:|---:|---:|
| complete `combo_sat.cfg` (29,826 graphs) | 350.58 s | 356.08 s | 1.6% slower |
| fixed `list1_sat.cfg` sample (1,000 graphs) | 221.38 s | 233.54 s | 5.5% slower |
| fixed `sat.cfg` sample (75 graphs) | 349.70 s | 355.77 s | 1.7% slower |

Verdict and timeout counts matched within every pair.  Thus the heuristic is
not a general OOPS optimization; its previously measured gain is specific to
the repeated cubic-flexibility queries.  To avoid proof-specific SAT search
tuning, the phase, polarity, and activity changes were removed from the cubic
verifier.

## 2026-07-17: existing crossed-first priority on cubic flexibility

The generic `cross-priority` SAT heuristic was temporarily enabled in the
cubic verifier and compared with both no priority and the removed
uncrossed-first heuristic.  Runs used one pinned core and identical input.

| sample | no priority | crossed-first | uncrossed-first |
|---|---:|---:|---:|
| order 26, 390 girth-five graphs (first pair) | 113.35 s | 95.59 s | 221.19 s |
| order 26, same graphs (reverse-order pair) | 151.75 s | 142.97 s | not repeated |
| order 24, 203 girth-five graphs | 484.54 s | 443.87 s | >600 s (151/203) |

The balanced order-26 pairs totalled 265.10 seconds without priority and
238.56 seconds crossed-first, a 10.0% speedup; the individual paired gains
were 15.7% and 5.8%.  Crossed-first was 8.4% faster in the single order-24
pair.  All completed runs produced every required certificate.
Uncrossed-first was strongly harmful on both samples, confirming that its
earlier gain was sample-dependent.  The temporary wiring was removed after
measurement; the verifier remains at the no-priority baseline pending a
decision about whether reusing the existing generic heuristic justifies a
workflow-specific parameter choice.

The remaining phase-saving choice was then measured on the same 390-graph
order-26 sample.  Crossed-first with normal phase saving took 123.08 seconds;
disabling phase saving took 284.01 seconds, 2.3 times as long.  The existing
generic `crossPriority` implementation—activity boost, `cross1=true` preferred,
and normal phase saving—is therefore enabled for cubic verification.  The
stale command-line and `Params` comments claiming that this option disabled
phase saving were corrected.

The dominant order-28 packed-5-cycle path was then measured directly on a
9,483-graph girth-five sample split into two halves.  To control for machine
drift, the first half ran baseline then crossed-first and the second half ran
crossed-first then baseline.

| sample | no priority | crossed-first | reduction |
|---|---:|---:|---:|
| half 0 (4,742 graphs) | 108.40 s | 78.28 s | 27.8% |
| half 1 (4,741 graphs) | 110.16 s | 88.92 s | 19.3% |
| combined | 218.56 s | 167.20 s | 23.5% |

Every graph received a packed-5-cycle certificate in all four runs.  Since
this is the largest projected production stage and both independently ordered
halves improved substantially, crossed-first priority with normal phase
saving remains enabled in the strict cubic verifier.

Replacing the two local recursive cycle searches in the order-26 5-cycle and
order-28 6-cycle contractions with the generic `forEachCycle` helper was also
tested before adoption.  Although both implementations found valid
certificates for every sampled graph, their different deterministic cycle
choices changed the downstream SAT cost.

| stage | original DFS | `forEachCycle` | change |
|---|---:|---:|---:|
| order-26 5-cycles, 390 graphs | 128.62 s | 140.88 s | 9.5% slower |
| order-28 6-cycles, 186 graphs | 82.65 s | 96.00 s | 16.1% slower |

Each total combines two matched residues run in opposite order.  The generic
helper refactor was rejected and the original local searches were restored.

The small `std::set<int>` objects used for cycle membership and distinct
external-neighbor checks were then replaced by 64-bit masks.  This preserves
the original DFS and cycle-selection order while removing heap allocation,
including inside the dominant packed-5-cycle enumeration.  Each contraction
now checks its `n <= 63` mask precondition.  The order-26/order-28 contraction
smoke cases, a 95-parent packed-5-cycle sample, and the finite patch tests all
pass with the representation-only change.

## 2026-07-17: paper-level statement of the upper bound

The `Upper Bound` section of the theorem draft was rewritten as a
self-contained research-paper argument.  It now states the order-at-most-28
goal, explains why direct recognition of the complete order-26 and order-28
censuses is infeasible, and presents the computation as five separate finite
claims with graph counts and positive-witness certification arguments.  The
summary theorem then derives 1-planarity of every simple cubic graph through
order 28 from those claims and the preceding reduction lemmas.  Project
workflow, source-code, command-line, and conditional-audit language was
removed.  The document compiles successfully.

## 2026-07-17: packed-5-cycle tail on the Minibaum `000/256` shard

Minibaum generated 1,221,746 order-28 graphs of minimum girth five in shard
`000/256`. Selecting every eighth record gave a deterministic 152,719-graph
sample. The current packed verifier on residue `0/48` stopped making
observable progress at extracted graph 88. This is zero-based sample record
4176 after accounting for the residue selection. Running that graph alone
used one CPU continuously for more than 180 seconds, with a maximum resident
set size of only 52 MiB; the process did not assert or exhaust memory.

The two simplification variants solve the same graph immediately:

| method | time |
|---|---:|
| packed 5-cycles | more than 180 s |
| one 5-cycle | 0.05 s |
| direct original graph | 0.03 s |

On the complete 3,182-graph residue `0/48`, one-5-cycle contraction verified
all graphs in 511.47 seconds: 3,150 used a 5-cycle and 32 used a 6-cycle.
Direct verification processed 1,598 of the 3,182 graphs before a 600-second
experimental cap, projecting to approximately 20 minutes if the remaining
records have similar cost. Thus packing several 5-cycles can create severe
SAT tails even though it is faster on the earlier 9,483-graph sample. The
one-5-cycle method is both substantially simpler and more robust on this
larger deterministic sample. All temporary tracing and variant changes were
removed after the experiment, and the baseline binary was rebuilt.

## 2026-07-17: runtime table recalibrated from server residues

The earlier paper table assumed that each target-server core would be 2.5
times faster than the development machine. Two completed target-server
residues do not support that assumption. The order-26 minimum-girth-five
residue `000/256` contains 60,221 graphs and took 645 seconds on 32 cores with
the selected 5-cycle contraction. Scaling its measured rate to the complete
31,478,584-graph census gives 124.9 equivalent core-days, or 3.9 days under
ideal 32-core distribution. Direct verification of the same residue took
1,440 seconds.

The order-28 minimum-girth-six residue `000/256` contains 20,442 graphs and
took 340 seconds on 32 cores with the selected 6-cycle contraction. Scaling
to 4,624,501 graphs gives 28.5 equivalent core-days, or 21.4 hours on 32
cores. Direct verification of that residue took 1,620 seconds.

The projected times for these two claims in `main_cubic.tex` were updated
accordingly. The order-28 girth-five projection and the total are now marked
as not established: the deterministic packed-5-cycle experiment above
invalidates the earlier tail-free mean, but does not yet provide a completed
replacement measurement.

## 2026-07-18: 5-cycle-to-path reduction replaces packed mode

All experiments in this entry used one OOPS process. Direct verification of
2,000 consecutive graphs from Minibaum order-28 minimum-girth-five residue
`000/256` took 209.96 seconds, or 104 milliseconds per graph. This establishes
the development baseline and a 5.2-millisecond target for a twenty-fold
speedup.

Two SAT-based alternatives were investigated first. Replacing simultaneous
fixed-edge obligations by existential selectors reduced a known packed-mode
tail from more than 180 seconds to 12 milliseconds. An exact reduced-graph
witness cache then reduced the 2,000-graph mean to 6.8 milliseconds. Although
sound, this remained short of the target and retained the complicated packed
reduction. Packed mode was therefore removed from the active verifier.

The replacement is structural. Delete a 5-cycle and introduce a path
`a-b-c`, attaching two consecutive outside neighbors to `a`, the next one to
`b`, and the remaining two to `c`. If the reduced order-26 graph has a drawing
in which the three edges incident with `b` are uncrossed, one of four checked
local patches restores the 5-cycle. The order-26 audit was correspondingly
strengthened: for every vertex it now verifies a drawing in which its complete
three-edge star is uncrossed. This property also implies the former one-edge
flexibility claim, so it replaces rather than supplements that job.

The C++ reduction scan took 95 milliseconds on the same 2,000 consecutive
parents. After the final connectivity checks and logging were in place, it
took 408 milliseconds on the first 10,000 parents and 245 milliseconds on a
disjoint consecutive block of 10,000 parents. All 20,000 graphs in the two
larger samples admitted a simple connected cubic order-26 reduction of girth
at least five. The observed 0.025--0.041 milliseconds per graph is over 2,500
times faster than direct recognition. At this point the intended production
evidence required every nonplanar order-28 graph of girth exactly five to use
this reduction. The later `102/256` residue experiment below disproved that
stronger expectation, so the final workflow records and validates direct
fallbacks instead.

As a broader tail check, the exact-girth-five portion of Minibaum residue
`000/256` was located at records 0 through 1,209,375; record 1,209,376 is the
first girth-six graph. One OOPS process scanned all 1,209,376 consecutive
girth-five records in 32.8 seconds. Every graph used the 5-cycle-to-path
reduction and none used the direct fallback. The resulting 0.027 milliseconds
per graph is about 3,800 times faster than the 104-millisecond direct baseline.

The stronger order-26 property is more expensive. It took 3.33 seconds on ten
consecutive order-26 girth-five graphs and 69.84 seconds on the first 100. The
latter rate is about twice the effective per-graph cost of the previous
server order-26 computation. This is the principal remaining scheduling
tradeoff of the much simpler order-28 proof and needs a matched server pilot.

## 2026-07-18: Minibaum residue `102/256` exposes path-reduction fallbacks

Minibaum generated 400,999 order-28 graphs of minimum girth five in residue
`102/256`. The first 400,471 records have girth exactly five; the remaining
528 have girth at least six. The boundary was checked directly on records
400,470 and 400,471.

One OOPS process verified the complete exact-girth-five prefix in 334.88
seconds. Of the 400,471 graphs, 398,671 admitted a simple connected cubic
5-cycle-to-path reduction of girth at least five. The remaining 1,800 graphs
(0.45 percent) did not admit such a reduction and were verified directly.
Every graph passed.

The overall mean was 0.836 milliseconds per graph, a 124-fold speedup over
the 104-millisecond direct baseline. Thus the twenty-fold performance target
still holds on this much less favorable residue. However, this experiment
disproves the draft assertion that every order-28 graph of girth exactly five
admits the path reduction. The final computational claim must explicitly
include verified direct drawings for the structural fallback, and the
production evidence must count those fallbacks rather than require zero.

## 2026-07-18: reuse of the tree replacement at other orders

The analogous simplification for Claim 5 was tested before changing the
verifier. A 6-cycle can be replaced by a cubic claw: its three leaves each
receive two of the six outside attachments, and the order-26 three-edge-star
claim keeps all three internal claw edges uncrossed. For each of the 15 ways
to pair the six cyclic attachments, all rotation systems of the claw were
enumerated. Local 6-cycle patches were then sought using crossings only among
the new cycle edges, since none of the six attachment edges is protected.

No pairing works for every rotation. The two best consecutive pairings each
leave four of their eight boundary orders uncovered; the other pairings leave
six or eight uncovered. Handling those rotations would require additional
protected edges or constraints on the rotation at the claw. This would add a
new SAT mode and a new case analysis, so the idea was rejected. The existing
6-cycle contraction remains worthwhile: it is more elaborate than direct
verification, but the matched server residue was about 4.8 times faster.

The 5-cycle-to-path replacement also does not directly simplify Claims 1--3.
Those claims require drawings preserving arbitrary sets of three edges, two
edges, or an entire three-edge star. Expanding a replaced 5-cycle already
requires the middle three-edge star to be uncrossed. Preserving the original
claim simultaneously would therefore require up to six protected edges in
the smaller graph, while the preceding claims guarantee only three or two.
Establishing the stronger properties would add computational jobs and local
lemmas rather than remove them. No such reduction is included in the final
workflow.

## 2026-07-18: final runtime table reconciliation

The Claim 3 row in `main_cubic.tex` still used the preliminary 100-graph local
measurement after the stronger three-edge-star property was adopted. The
completed Minibaum residue `000/256` contains 60,221 graphs and took 1,440
seconds on 32 server cores. Scaling this measured rate to all 31,478,584
graphs gives 278.8 equivalent core-days, or 8.7 days on 32 cores. Claims 1
and 2 retain their latest single-core pilot estimates; Claims 4 and 5 already
use the completed residues recorded above. The reconciled total is about 325
core-days, or 10.2 days under ideal 32-core distribution.

## 2026-07-18: equivalent 6-cycle labeling accelerates Claim 5

Several ways to reuse Claim 3 for the girth-six family were rejected before
changing production code. A four-vertex path has three internal edges, but no
single three-edge star protects all three, so it has no crossing-free regular
neighborhood. The cubic-claw replacement had already failed its rotation
test. An exhaustive local search then tested all 20 choices of three protected
6-cycle attachments against all 60 boundary orders; none covered every order.
The same search covered all 60 orders with the existing four protected
attachments, reproducing the checked lemma. Among 232 sampled graphs, no
alternative 6-cycle produced a planar contraction. An unconstrained-first
SAT solve took 129.5 seconds versus the 126.0-second constrained baseline.

The 6-cycle lemma leaves one opposite pair of attachments unprotected, and
there are three equivalent choices. Their times on the same 232-graph sample
were 126.0, 67.8, and 51.0 seconds. The best choice was then compared with the
original on a disjoint 232-graph sample: 97.2 versus 217.1 seconds. In the
final implementation the selected DFS cycle is rotated by one position and
the existing protected-position set `{1,2,4,5}` is unchanged. Thus the local
table, proof, contraction, and number of SAT assumptions are unchanged; only
the association between graph edges and SAT variable indices differs.

Across the two matched samples the rotation reduces 343.1 seconds to 148.5
seconds, a 2.31-fold speedup. Applying this factor to the completed server
residue revises the Claim 5 projection from 28.5 to 12.3 core-days, or from
21.4 to 9.2 hours on 32 cores. Claims 1, 2, and 5 now project to 17.7 hours
combined. A server rerun of one complete residue remains necessary to replace
this ratio-based projection with a direct measurement.

## 2026-07-18: final deletion-oriented verification audit

The final implementation was reviewed against a requirement to retain at
least 90 percent of end-to-end performance. No further core OOPS mechanism can
be removed within that budget. Rebuilding the SAT formula for each requirement
was previously 3.49 times slower than repeated assumptions. The full existing
SAT encoding was about 23--28 percent slower on the constrained verification
samples than the small crossing-variable encoding. Removing either order-28
cycle reduction also exceeds the ten-percent total-runtime allowance. The
remaining core surface is therefore the verification entry point, the lean
crossing-variable mode, public construction/solve helpers, and the correctness
repair for repeated assumptions.

The evidence layer admitted a substantial no-cost reduction. Claim 1 has only
1,575,835 records and can be generated directly by the ten documented nauty
commands, so the 237-line resumable Claim 1 generator was deleted. The
187-line offline evidence auditor duplicated hash, count, and log checks made
during generation and verification; it was also deleted. The verification
driver retains restartable execution, exact theorem counters, binary and input
hashes, and full binding to the Minibaum manifests needed for the much larger
orders 24, 26, and 28. This reduces the new Python implementation from 1,057
to 568 lines without affecting any SAT computation.

Finally, the separate local-6-cycle and 6-cycle-expansion lemmas in the
paper were merged. The same exhaustive 60-order certificate remains, but the
argument now has one statement and one computational proof rather than two
conceptual layers.

The standalone repeated-assumption regression was also removed after a final
one-time stress campaign. Five million generated seven-variable formulas were
tested with up to 200 alternating SAT and UNSAT assumption queries per
satisfiable formula on one reused solver. Every result was compared with
exhaustive truth-table evaluation and every SAT model was checked against its
assumptions. The campaign completed in 293.6 seconds without a discrepancy.
The production verifier still reconstructs every SAT drawing and checks every
required uncrossed edge, while a false UNSAT answer can only abort a run;
therefore the permanent 147-line solver regression is not part of the minimal
verification package.
## 2026-07-18: one paper and one verification checklist

The verification presentation was reduced to two authoritative documents.
`main_cubic.tex` is now exclusively a self-contained mathematical paper: it
defines the graph-theoretic terminology used in the argument and contains no
program names, commands, filenames, or implementation discussion.  The old
mixed proof-and-operations guide was replaced by a concise claim-by-claim
checklist in `cubic/verification.md` containing the exact generation,
verification, counting, and hashing commands.

The two Python orchestration programs, `verify_cubic28.py` and
`generate_minibaum_census.py`, were deleted.  Their mathematical checks were
already enforced by the strict OOPS verification mode, while deterministic
Minibaum residue generation, restart at file boundaries, record counts, and
hashing can be expressed directly in the checklist.  The test suite retains
the strict OOPS smoke test and the independent finite local-drawing test.

## 2026-07-18: first attempt to weaken the order-26 requirement

The 5-cycle expansion was checked to determine whether only the two
internal path edges, rather than all three edges incident with the middle
vertex, need be uncrossed.  Three of the four boundary orders already have
local drawings whose crossings use only new 5-cycle edges.  For the remaining
boundary order `01342`, all eleven matchings of nonincident 5-cycle-edge
pairs were exhaustively tested with the same planarization, alternating-order,
and connected-bridge checks used by the permanent local-drawing regression.
None gives a valid local drawing.  The existing patch crosses the middle
attachment in this order.  Thus the three-edge-star requirement cannot be
replaced by the two internal path edges without a different reduction or an
additional condition on the drawing.

## 2026-07-18: cycle terminology standardized

The paper, verification instructions, active C++ implementation, diagnostics,
and local-expansion test now use the graph-theoretic terms `4-cycle`,
`5-cycle`, and `6-cycle` consistently.  Reduction details such as replacing a
5-cycle by a path remain in the definitions rather than in lemma or certificate
names.  Historical command-line options, filenames, and source identifiers are
preserved verbatim in earlier log entries.

## 2026-07-18: final implementation and notation audit

The production verifier and the finite local-drawing test were reread against
the exact commands in `verification.md` and the notation in `main_cubic.tex`.
The verifier no longer stores the selected 6-cycle after its vertex mask and
outside neighbors have been recorded.  Claim 3 now constructs the incident
edge sets in one pass over the edges, in exactly the same order as before, and
does not revalidate sets that are implied by the preceding cubic-input check.
The input check uses the adjacency lists already supplied by the parser.

Implementation names now follow the paper: the replacement path has vertices
`a`, `b`, and `c`; outside neighbors and attachment edges are distinguished;
drawings are “verified”; and the order-28 return value records the verification
method rather than calling it a certificate.  In the local-drawing test,
“patch” and “spoke” were replaced by “local drawing” and “attachment edge.”
An unused omitted-vertex argument was removed from the connectivity check, and
the generic table checker used only for the 6-cycle was replaced by the exact
6-cycle check.

The local-drawing test, the 27-graph verification smoke test, and a full build
all pass.  On the same 19 graphs from residue class `0/10000` of
`cub26-gir6.g6`, the Claim 3 runtime changed from 468.11 seconds before the
cleanup to 462.61 seconds after it.  This 1.2-percent difference is within
run-to-run variation and shows no performance regression.
