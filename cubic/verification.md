# Verification of the five computational claims

This document gives the commands needed to reproduce the five computational
claims in `main_cubic.tex`.  The mathematical definitions, reductions, and
proofs are given only in that paper.

Run all command blocks in the same Bash session.  The commands below are
serial for clarity.  On a parallel machine, distinct
input files may be assigned to independent processes.  A verification is
successful only if every process exits successfully and the total number of
processed records equals the stated family size.

## Preparation

Build OOPS and Minibaum, and identify the nauty executables:

```bash
set -euo pipefail

make -j
gcc -O3 /path/to/minibaum5.c -o minibaum5

NAUTY=/path/to/nauty2_9_3
MINIBAUM=./minibaum5
PARTS=256

mkdir -p data evidence
sha256sum ./oops "$MINIBAUM" /path/to/minibaum5.c \
  "$NAUTY/geng" "$NAUTY/planarg" > evidence/programs.sha256
```

The following shell function generates one complete Minibaum enumeration in 256
deterministic residue classes.  A completed file is not regenerated.

```bash
generate_minibaum() {
  n=$1
  girth=$2
  output=$3
  mkdir -p "$output"
  for part in $(seq 0 $((PARTS - 1))); do
    file=$(printf '%s/part-%03d.g6' "$output" "$part")
    if [ ! -f "$file" ]; then
      "$MINIBAUM" "$n" "$girth" s g m "$part" "$PARTS" > "$file.tmp"
      mv "$file.tmp" "$file"
    fi
  done
}
```

Before the graph families are processed, verify the finite local drawings
used in the 4-cycle, 5-cycle, and 6-cycle expansion lemmas:

```bash
g++ -std=c++17 -O2 -Isrc \
  cubic/test_cycle_expansions.cpp src/planarity_test.cpp \
  -o test_cycle_expansions
./test_cycle_expansions
```

## Verify Claim 1

Generate every connected cubic graph of orders 4 through 18, and every
biconnected, nonplanar, triangle-free cubic graph of orders 20 and 22:

```bash
mkdir -p data/claim1 evidence/claim1

for n in 4 6 8 10 12 14 16 18; do
  "$NAUTY/geng" -cq -d3 -D3 "$n" > "data/claim1/cub${n}.g6"
done

for n in 20 22; do
  "$NAUTY/geng" -Cq -t -d3 -D3 "$n" |
    "$NAUTY/planarg" -v > "data/claim1/cub${n}-core.g6"
done

wc -l data/claim1/*.g6

for input in data/claim1/*.g6; do
  name=$(basename "$input" .g6)
  ./oops -i="$input" -verify-cubic -colors=0 \
    > "evidence/claim1/${name}.log" 2>&1
done
```

The expected record counts, in increasing order, are

```text
1, 2, 5, 19, 85, 509, 4,060, 41,301, 97,141, 1,432,712.
```

Their sum is 1,575,835.  Every nonplanar record must be reported as
3-flexible; planar records satisfy the claim by a planar drawing.

## Verify Claim 2

Generate the connected cubic graphs with $n=24$ and girth at least 5, then
verify 2-flexibility:

```bash
generate_minibaum 24 5 data/claim2
wc -l data/claim2/*.g6

mkdir -p evidence/claim2
for input in data/claim2/*.g6; do
  name=$(basename "$input" .g6)
  ./oops -i="$input" -verify-cubic -colors=0 \
    > "evidence/claim2/${name}.log" 2>&1
done
```

The expected number of records is 1,620,479.  Every nonplanar record must be
reported as 2-flexible.

## Verify Claim 3

Generate the connected cubic graphs with $n=26$ and girth at least 5.  For
each graph and each vertex, verify the existence of a 1-planar drawing in
which the three incident edges are uncrossed:

```bash
generate_minibaum 26 5 data/claim3
wc -l data/claim3/*.g6

mkdir -p evidence/claim3
for input in data/claim3/*.g6; do
  name=$(basename "$input" .g6)
  ./oops -i="$input" -verify-cubic -colors=0 \
    > "evidence/claim3/${name}.log" 2>&1
done
```

The expected number of records is 31,478,584.  For every nonplanar graph and
every vertex, the computation must find a 1-planar drawing in which the three
edges incident with that vertex are uncrossed.

## Verify Claim 4

Generate every connected cubic graph with $n=28$ and girth at least 5.  This
same computation will also verify Claim 5:

```bash
generate_minibaum 28 5 data/claims4-5
wc -l data/claims4-5/*.g6

mkdir -p evidence/claims4-5
for input in data/claims4-5/*.g6; do
  name=$(basename "$input" .g6)
  ./oops -i="$input" -verify-cubic -colors=0 \
    > "evidence/claims4-5/${name}.log" 2>&1
done
```

The complete enumeration has 656,783,890 records.  Exactly 652,159,389 have girth
5.  Each of these must be planar, verified by a 5-cycle expansion,
or verified directly.

## Verify Claim 5

Use the input family and the computation from Claim 4; no second OOPS run is
required.  Summed over all logs, the required results are:

- 4,624,480 records verified by a 6-cycle expansion; and
- 21 records of girth at least 7 verified directly.

Thus the final two quantities sum to the 4,624,501 records of girth at least
6 required by Claim 5.

## Final checks

Every OOPS invocation must exit successfully.  Sum `#planar` and
`#1-planar` over the logs and compare the result with the appropriate family
size above.  Finally, record hashes of all inputs and logs:

```bash
find data evidence -type f ! -name files.sha256 -print0 | sort -z | xargs -0 sha256sum \
  > evidence/files.sha256
```
