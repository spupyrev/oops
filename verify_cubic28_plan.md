# Verify 1-planarity of all cubic graphs with n <= 28 vertices

## Project Goal

The goal of the project is to (computationally + analytically) verify that ALL cubic
graphs with up to 28 vertices are 1-planar. The verification process should be
reproducible and take less than 24h of compute time on a 32-core machine.

There are 40,497,138,011 such graphs with n=28 and 2,094,480,864 with n=26. We do know 
that all smaller instances (n<=24) are 1-planar.

## Background

OOPS is a SAT-based solver for determining 1-planarity of a graph. It is pretty fast: 
takes ~100ms to test 1-planarity of a 24-vertex cubic graph and <1sec for a 30-vertex one.
Based on earlier experiments with oops, we know that ALL cubic graphs having n<=24 nodes are
1-planar. The utlimate goal is to find the smallest cubic non-1-planar instance.
We already know that there are 2 instances with n=30 nodes. What remains to do is to test
all smaller instances with n<=28, which is challenging due to the large number of them.

We have a powerful multi-core machine which we could use for running oops in parallel.
However, testing all 40*10^9 instances would take more than thousand of machine-days, which
is impractical. Hence, the goal is to come up with alternative strategy. At this moment, we
are open to all ideas; see below for specifics.

We work with simple, connected, unlabeled cubic graphs without loops and parallel edges.

## Data And Baselines

Build oops with `make -j`; see README.md for details.

To run oops on a single instance, use

```bash
./oops -i=benchmarks/3reg-gir7.g6 -verbose=1 -part=0
```

To run oops on a collection of instances, use

```bash
./oops -i=data/sat.cfg -verbose=0
```

An extensive benchmark of cubic graphs is available in `/home/spupyrev/research/one_planar/data/cubic`. For example, to test all instances of cubic graphs with n=18, use

```bash
./oops -i=/home/spupyrev/research/one_planar/data/cubic/cub18.g6 -verbose=0
```

It takes ~10 min (on a single thread). To run the same job on 6 threads, use

```bash
python3 ./scripts/proc.py
```

It completes in ~270 sec. 
One may need to adjust/extend the script to test other benchmarks or pass alternative options.

Finally, we should only look at bi-connected non-planar instance, since the smallest instance 
cannot be disconnected or one-connected:
```bash
./oops -i=/home/spupyrev/research/one_planar/data/cubic/cub18_Ctd3D3_nonplanar.g6
```
which ends in 3 minutes on one thread.

More data is available in `data/cubic/`; for example, `data/cubic/cub26_0.zip` contains
a subset (~1%) of all cubic graphs with n=26. `data/cubic/cub28-b.g6` are all bipartite
cubic graphs with n=28. Do use existing datasets for testing; when we have a solid solution and
ready for a production run, we can generate actual full benchmarks from plantri.

## Ideas

1. The primary idea so far was to filter out instances that are provably 1-planar or provably 
not smallest non-1-planar (since we're looking for a smallest non-1-planar cubic graph).
For example, we can show and provide a theoretical proof that such an instance contains no
triangle face. See main_cubic.tex for (incomplete) details.
The corresponding implementation is in reducible_subgraphs.cpp. To use, run
```bash
./oops -i=/home/spupyrev/research/one_planar/data/cubic/cub18.g6 -skip-reducible
```
which completes much faster, around 3 min.
The idea is to design (and show correctness!) more reducible subgraphs that can substantially
speed up the process for n=26 and n=28.

2. Currently we use plantri to generate cubic graphs. An idea would be intercept the generation
process, prove that certain branches/instances are not needed, and test with oops 
fewer instances. This is an open, exploratory idea.

3. Currently oops skips planar instances, since they're clearly 1-planar and planarity is a
very fast check. Can we develop any other quick filters? This is an open, exploratory idea.

## Process

Stricly follow the following development:

1. Come up with an idea. Write it down to cubic/cubic28_progress_report.md
2. Read the code, think about the best possible implementation, decide what logging is
   needed.
3. Implement the feature, run smallish SAT/UNSAT tests. Iterate over the implementation
   until all tests pass as expected.
5. Iterate until the feature is fully implemented and thoroughly tested. If there are
   promising follow-ups, continue to Step 2.
6. Record the results and the analysis to cubic/cubic28_progress_report.md
7. Go to Step 1

You have at most 48h for research/development on the project, so be careful with 
running long tests. The actual testing time of real instances is excluded from the timeline.

YOU CAN NEVER RUN MULTI-THREADED OOPS JOBS BY YOURSELF; A HUMAN WILL DO IT ON A DEDICATED SERVER

## Result

As a result of the work, I expect a reproducible verification flow, e.g., 
a new version of oops binary and supporting documentation with correctness proofs and
explanations. An independent reviewer should be able to re-run the verification process
on their machine and be convinced that all instances with n<=28 are 1-planar.