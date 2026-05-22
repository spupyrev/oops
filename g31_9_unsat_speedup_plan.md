# g.31.9 UNSAT Speedup

## Project Goal

```bash
./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=27 -up-sepcycles-timeout=600
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
[add `-up-sepcycles` as an option for an alternative solving strategy]


| Part | Graph     | solver time (no -up-sepcycles) | solver time (with -up-sepcycles) |
|------|-----------|--------------------------------|----------------------------------|
| 13   | Hoffman   |  19 sec                        | 467 ms                           |
| 18   | Robertson | 103 sec                        |   3 sec                          |
| 22   | Folkman   |  69 sec                        |  34 sec                          |
| 23   | Brinkmann | 231 sec                        |   4 sec                          |
| 30   | Holt      |   7 min                        |   5 sec                          |

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
| 74   | g6_365    |  5 ms                          |  3 ms                            |
| 90   | g6_391    |  4 ms                          |  3 ms                            |
| 267  | g6_158    | 59 ms                          | 14 ms                            |
| 272  | g6_163    | 23 ms                          | 20 ms                            |

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
| 15   | grafo7753.40.graphml | Timeout (179 sec)   | UNSAT (5 min)         |
| 4    | grafo3365.43.graphml | UNSAT (136 sec)     | UNSAT (128 sec)       |
| 6    | grafo6547.38.graphml | UNSAT (145 sec)     | UNSAT (73 sec)        |

### UNSAT benchmark cases (unsat.cfg):

```bash
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=27 -timeout=180
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=27 -up-sepcycles -timeout=180
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=26 -timeout=180
  ./oops -i=data/unsat.cfg -verbose=1 -unsat=1 -part=26 -up-sepcycles -timeout=180
```

| Part | Graph           | result/time (no UP) | result/time (with UP) |
|------|-----------------|---------------------|-----------------------|
| 27   | g.31.9.graphml  | Timeout (176 sec)   | UNSAT (8 min)         |
| 26   | g.47.27.graphml | UNSAT (34 sec)      | UNSAT (32 sec)        |


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
| list1_sat.cfg  | 2195  | list_2000_graphs.g6_552  | SAT (10 sec) | SAT (34 sec)  |
| sat.cfg        | 73    | grafo4911.57.graphml     | SAT (141 sec) | SAT (159 sec) |
| sat.cfg        | 78    | g.59.0.graphml           | Timeout (191 sec) | Timeout (182 sec) |
| sat.cfg        | 95    | TC_minus_0_17            | SAT (141 sec) | SAT (105 sec) |
| combo_sat.cfg  | 21511 | cub30-gir7.g6_544        | SAT (141 sec) | Timeout (179 sec) |


Be mindful of the runtimes: some tests take substantial amount of time, so constantly
running tests might quickly exceed the allocated 24h of development time.


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
