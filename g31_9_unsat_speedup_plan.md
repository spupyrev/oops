# g.31.9 UNSAT Speedup

## Project Goal

```bash
./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=27 -up-sepcycles -timeout=600
```

must finish with UNSAT for `g.31.9.graphml` within 600 seconds of raw SAT solver time.

The project is successful only if the implementation produces an UNSAT result for this target instance, while existing SAT instances remain SAT.

You have at most 24h to work on the project, so be careful with running long tests.

## Instance And Baseline

- Graph: `g.31.9.graphml`
- Size: `|V| = 31`, `|E| = 58`
- SAT model: 1355 possible `cross2` variables and 58 `cross1` variables
- Encoding time is about 27 seconds, so the target is dominated by SAT solving.

## Build And Test Commands

Build from the repository root:

```bash
make -j
```

UNSAT benchmark cases:

```bash
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=13 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=18 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=22 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=23 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=30 -timeout=600
```
[optionally add `-up-sepcycles` as an option for an alternative solving strategy]
 

| Part | Graph     | SAT solver took | SAT solver took (with -up-sepcycles) |
|------|-----------|-----------------|--------------------------------------|
| 13   | Hoffman   |  13 sec         |  1 sec                               |
| 18   | Robertson |  96 sec         |  Time Limit Exceeded                 |
| 22   | Folkman   |  63 sec         | 17 sec                               |
| 23   | Brinkmann | 133 sec         |  7 sec                               |
| 30   | Holt      |  10 min         |  5 sec                               |

SAT benchmark cases:

```bash
  ./oops -i=data/combo_sat.cfg -verbose=1 -unsat=1 -up-sepcycles -timeout=120
  ./oops -i=data/sat.cfg -verbose=1 -unsat=1 -up-sepcycles -timeout=120
```

Be mindful of the runtimes: some tests take substantial amount of time, so constantly 
running tests might quickly exceed the allocated 24h of development time.

## Available Levers

This is a list of possible directions to explore. Do not limit yourself to the list, explore
all options. It is very likely that these ideas alone are not sufficient to reach the goal.

1. Improve UP strategy [-up-sepcycles] 
   SAT search control: Decide which variables Glucose branches on and how strongly `cross2` 
   participates in the decision heuristic. This is based on an earlier observation that the 
   solver spends too much search away from useful crossing choices.

   See Direction 1/2/3 below as potential next steps.

  Do not treat UP as the whole project. UP is one lever. The larger problem is making CDCL 
  search the right crossing space.

2. Explore alternative SAT solvers. 
   Check ../sat_solvers/ folder with possible alternatives. We can easily test alternatives
   via dumping CNF to dimacs files (-dimacs=...), and running a solver on the dump.

3. Improving existing SAT encoding.
   It is hard but not impossible to add more constraints/variables to the model. In 
   particular, adding 2-clauses (or 1-clauses!) would be very useful, as it directly reduces 
   the search space.
   Added constraints/variables should be sound, and based on a theoretical proof that the 
   resulting encoding is correct.

4. Incremental solving.
   One idea that has not been explored yet: Start with an encoding, run the solver for some 
   period of time, then stop and analyze learnt clauses. If there are new 1-clauses and 
   2-clauses (and perhaps, some 3-clauses), it might be worth re-starting the solver from 
   scratch by adding the BASE encoding with these extra small clauses.

## Process

Stricly follow the following development:

1. Come up with an idea. Write it down to g31_9_unsat_progress_report.md
2. Read the code, think about the best possible implementation, decide what logging is needed.
3. Implement the feature, run smallish SAT/UNSAT tests. Iterate over the implementation 
   until all tests pass as expected.
4. Run g.31.9.graphml, possibly with smallish timeout (e.g., -timeout=60), analyze logs
5. Iterate until the feature is fully implemented and thoroughly tested. If there are 
   promising follow-ups, continue to Step 2.
6. Record the results and the analysis to g31_9_unsat_progress_report.md
7. Go to Step 1


## Earlier (stale) directions for improving UP

## Direction 1: Cross2-Centered SAT Search

Hypothesis:

Glucose is spending too much search away from meaningful `cross2` choices. UP cannot help if the solver rarely reaches relevant crossing assignments.

Implementation direction:

- Add a generic solver mechanism for decision eligibility and/or initial activity.
- The model builder marks solver variables by family, but Glucose only receives generic weights or decision flags.
- First test: keep `cross2` variables as preferred decisions and prevent irrelevant auxiliary variables from being selected as decisions when possible.
- Keep `cross1` secondary. It may remain branchable, but it should not outrank `cross2`.

Required logging:

- decisions by variable family: `cross2`, `cross1`, auxiliary;
- decision level of first `cross2=true`;
- number of active `cross2=true` at UP checks and UP conflicts;
- conflicts/sec, because strict priority may slow propagation even if it improves search shape.

Keep criteria:

- `cross2` decisions increase substantially;
- UP pending calls increase for actual `cross2=true` assignments;
- 600s run shows clear progress toward UNSAT or solves the instance.

Stop criteria:

- `cross2` decisions increase but UP conflicts/search quality do not improve;
- conflicts/sec collapses without better proof progress;
- behavior becomes sensitive only to arbitrary priority constants.

## Direction 2: Targeted Static Cross2 Pruning

Hypothesis:

Some crossing exclusions should be available before search. Static clauses are stronger than dynamic rediscovery, but only if they can be generated with tight filters.

Implementation direction:

- Extend static sep-cycle/equal-flow preprocessing only for targeted short `cross2` clauses.
- Use longer cycles only with DFS pruning and degree/path bounds.
- Prefer pure `cross2` clauses or very short mixed clauses.
- Do not broadly enumerate all 3-clauses.

Metrics:

- encoding time;
- number of added clauses by size and type;
- SAT solver time;
- whether solving improves on g.31.9, not just small graphs.

Keep criteria:

- encoding time remains acceptable;
- added clauses substantially reduce SAT time;
- the clause family has a clear mathematical explanation.

Stop criteria:

- clause count grows by millions without clear SAT-time improvement;
- the rule depends on one-off observations from g.31.9 that do not have a sound general statement;
- the work becomes a maxflow/DFS micro-optimization project.

## Direction 3: Selective Dynamic Cross2=false

Hypothesis:

Conflict-only UP rejects bad crossing combinations too late. The useful dynamic implication target is `cross2=false`, but only if the candidate set is tiny and separator-local.

Implementation direction:

- Extend the generic user-propagator interface only enough to return one implied literal plus a reason clause.
- Trigger only from pending `cross2=true`.
- When a separator is found, reuse that separator to test nearby candidate crossings.
- For each candidate, verify soundness by temporarily assuming candidate `cross2=true` and running the separator proof.
- Return only `cross2=false` implications or conflict clauses. Do not reintroduce `cross1=true`.
- Keep reasons minimal enough to be useful: include only literals that support the separator proof.

Metrics:

- candidates tested per UP call;
- successful `cross2=false` implications;
- exclusions per successful separator;
- reason sizes;
- UP time;
- effect on decisions, conflicts, and result.

Keep criteria:

- candidate count stays small;
- more than one useful exclusion per separator on average;
- UP remains a small fraction of solver time;
- added `cross2=false` clauses materially reduce decisions or solve g.31.9.

Stop criteria:

- UP time again dominates solving;
- many `cross2=false` implications are added but the 600s result does not improve;
- candidate selection becomes broad/global.
