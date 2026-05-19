# g.31.9 UNSAT Speedup

## Project Goal

```bash
./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=27 -up-sepcycles -timeout=600
```

must finish with UNSAT for `g.31.9.graphml` within 600 seconds of raw SAT solver time.

The project is successful only if the implementation produces an UNSAT result for this
target instance, while existing SAT instances remain SAT.

You have at most 24h to work on the project, so be careful with running long tests.

## Instance And Baseline

- Graph: `g.31.9.graphml`
- Size: `|V| = 31`, `|E| = 58`
- SAT model: 1355 possible `cross2` variables and 58 `cross1` variables
- Encoding time is about 27 seconds, so the target is dominated by SAT solving.

Current result:
```bash
  SAT encoding took 29 seconds
  SAT solver took 18 minutes
  graph 'g.31.9.graphml' (index 0) with |V| = 31 and |E| = 58 is not 1-planar
```

## Build And Test Commands

Build from the repository root:

```bash
make -j
```

### UNSAT benchmark cases (named_paper.cfg):

```bash
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=13 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=18 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=22 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=23 -timeout=600
  ./oops -i=data/named_paper.cfg -verbose=1 -unsat=1 -part=30 -timeout=600
```
[optionally add `-up-sepcycles` as an option for an alternative solving strategy]


| Part | Graph     | solver time (no -up-sepcycles) | solver time (with -up-sepcycles) |
|------|-----------|--------------------------------|----------------------------------|
| 13   | Hoffman   |  13 sec                        | 576 ms                           |
| 18   | Robertson |  98 sec                        |   7 sec                          |
| 22   | Folkman   |  62 sec                        |  25 sec                          |
| 23   | Brinkmann | 134 sec                        |   6 sec                          |
| 30   | Holt      |  Timeout (540 sec)             |   6 sec                          |

### UNSAT benchmark cases (list1_unsat.cfg):

```bash
  ./oops -i=data/list1_unsat.cfg -verbose=1 -unsat=1 -part=62 -timeout=60
  ./oops -i=data/list1_unsat.cfg -verbose=1 -unsat=1 -part=74 -timeout=60
  ./oops -i=data/list1_unsat.cfg -verbose=1 -unsat=1 -part=90 -timeout=60
  ./oops -i=data/list1_unsat.cfg -verbose=1 -unsat=1 -part=267 -timeout=60
  ./oops -i=data/list1_unsat.cfg -verbose=1 -unsat=1 -part=272 -timeout=60
```

| Part | Graph     | solver time (no -up-sepcycles) | solver time (with -up-sepcycles) |
|------|-----------|--------------------------------|----------------------------------|
| 62   | g6_344    |  3 ms                          |  2 ms                            |
| 74   | g6_365    |  3 ms                          |  3 ms                            |
| 90   | g6_391    |  6 ms                          |  7 ms                            |
| 267  | g6_158    | 49 ms                          | 41 ms                            |
| 272  | g6_163    | 43 ms                          | 33 ms                            |

### UNSAT benchmark cases (rome50_unsat.cfg):

```bash
  ./oops -i=data/rome50_unsat.cfg -verbose=1 -unsat=1 -part=15 -timeout=180
  ./oops -i=data/rome50_unsat.cfg -verbose=1 -unsat=1 -part=15 -up-sepcycles -timeout=180
  ./oops -i=data/rome50_unsat.cfg -verbose=1 -unsat=1 -part=4 -timeout=180
  ./oops -i=data/rome50_unsat.cfg -verbose=1 -unsat=1 -part=4 -up-sepcycles -timeout=180
  ./oops -i=data/rome50_unsat.cfg -verbose=1 -unsat=1 -part=6 -timeout=180
  ./oops -i=data/rome50_unsat.cfg -verbose=1 -unsat=1 -part=6 -up-sepcycles -timeout=180
```

| Part | Graph                | result/time (no UP) | result/time (with UP) |
|------|----------------------|---------------------|-----------------------|
| 15   | grafo7753.40.graphml | Timeout (181 sec)   | Timeout (179 sec)     |
| 4    | grafo3365.43.graphml | UNSAT (111 sec)     | UNSAT (96 sec)        |
| 6    | grafo6547.38.graphml | UNSAT (71 sec)      | UNSAT (49 sec)        |

### UNSAT benchmark cases (unsat.cfg):

```bash
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=27 -timeout=180
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=27 -up-sepcycles -timeout=180
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=26 -timeout=180
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=26 -up-sepcycles -timeout=180
```

| Part | Graph           | result/time (no UP) | result/time (with UP) |
|------|-----------------|---------------------|-----------------------|
| 27   | g.31.9.graphml  | Timeout (299 sec)   | Timeout (178 sec)     |
| 26   | g.47.27.graphml | UNSAT (21 sec)      | UNSAT (45 sec)        |


### SAT benchmark cases:

```bash
  ./oops -i=../data/list1_sat.cfg -verbose=1 -unsat=1 -part=2195 -timeout=180
  ./oops -i=../data/list1_sat.cfg -verbose=1 -unsat=1 -part=2195 \
    -up-sepcycles -timeout=180

  ./oops -i=../data/sat.cfg -verbose=1 -unsat=1 -part=73 -timeout=180
  ./oops -i=../data/sat.cfg -verbose=1 -unsat=1 -part=73 \
    -up-sepcycles -timeout=180
  ./oops -i=../data/sat.cfg -verbose=1 -unsat=1 -part=78 -timeout=180
  ./oops -i=../data/sat.cfg -verbose=1 -unsat=1 -part=78 \
    -up-sepcycles -timeout=180
  ./oops -i=../data/sat.cfg -verbose=1 -unsat=1 -part=95 -timeout=180
  ./oops -i=../data/sat.cfg -verbose=1 -unsat=1 -part=95 \
    -up-sepcycles -timeout=180

  ./oops -i=../data/combo_sat.cfg -verbose=1 -unsat=1 -part=21511 \
    -timeout=180
  ./oops -i=../data/combo_sat.cfg -verbose=1 -unsat=1 -part=21511 \
    -up-sepcycles -timeout=180
```

| Config         | Part  | Graph                    | no UP        | with UP       |
|----------------|-------|--------------------------|--------------|---------------|
| list1_sat.cfg  | 2195  | list_2000_graphs.g6_552  | SAT (23 sec) | SAT (88 sec)  |
| sat.cfg        | 73    | grafo4911.57.graphml     | SAT (126 sec) | SAT (125 sec) |
| sat.cfg        | 78    | g.59.0.graphml           | SAT (32 sec) | Timeout (189 sec) |
| sat.cfg        | 95    | TC_minus_0_17            | SAT (58 sec) | SAT (85 sec)  |
| combo_sat.cfg  | 21511 | cub30-gir7.g6_544        | SAT (54 sec) | SAT (78 sec)  |


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
   particular, adding 2-clauses (or 1-clauses!) would be very useful, as it
   directly reduces
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
2. Read the code, think about the best possible implementation, decide what logging is
   needed.
3. Implement the feature, run smallish SAT/UNSAT tests. Iterate over the implementation
   until all tests pass as expected.
4. Run g.31.9.graphml, possibly with smallish timeout (e.g., -timeout=60), analyze logs
5. Iterate until the feature is fully implemented and thoroughly tested. If there are
   promising follow-ups, continue to Step 2.
6. Record the results and the analysis to g31_9_unsat_progress_report.md
7. Go to Step 1
