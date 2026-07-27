# Verification of the five computational claims

This document gives the commands needed to reproduce the five computational
claims in `main_cubic.tex`.  The mathematical definitions, reductions, and
proofs are given only in that paper.

Run all command blocks in the same Bash session.  The commands are serial for
clarity, as a proof-of-concept reproduction rather than a production
scheduler.  A verification is successful only if every command exits
successfully and the total number of processed records equals the stated
family size.

## Prerequisites

The verification uses the following software:

- A C++17 compiler, a C compiler, GNU Make, Bash, and standard Unix utilities
  including `awk`, `sha256sum`, `find`, `sort`, `xargs`, and `wc`.  Building
  OOPS also requires the zlib development headers and library.
- **nauty and Traces 2.9.3**, which supplies the `geng`, `pickg`, and `planarg`
  executables used for Claims 1--5.  Obtain the sources from the
  [official nauty page](https://users.cecs.anu.edu.au/~bdm/nauty/) or download
  the [nauty 2.9.3 source archive](https://users.cecs.anu.edu.au/~bdm/nauty/nauty2_9_3.tar.gz).
  After extracting the archive, run `./configure` and `make`; `geng`, `pickg`,
  and `planarg` will be created in the resulting `nauty2_9_3`
  directory.
- **Minibaum 5**, which generates the connected cubic graphs of prescribed
  minimum girth for Claims 2--5.  The author's
  [Minibaum page](https://caagt.ugent.be/minibaum/) provides the
  [minibaum5.c source](https://caagt.ugent.be/minibaum/minibaum5.c) and the
  [Minibaum manual](https://caagt.ugent.be/minibaum/minibaummanual.pdf).

## Preparation

Build OOPS and Minibaum, and identify the nauty executables:

```bash
set -euo pipefail

make -j
gcc -O3 /path/to/minibaum5.c -o minibaum5

NAUTY=/path/to/nauty2_9_3
MINIBAUM=./minibaum5
PARTS=256
export LC_ALL=C

g++ -std=c++17 -O3 -Isrc -I"$NAUTY" \
  cubic/prepare_cubic.cpp "$NAUTY/nauty.a" \
  -o prepare_cubic

mkdir -p data evidence
sha256sum ./oops ./prepare_cubic "$MINIBAUM" /path/to/minibaum5.c \
  "$NAUTY/geng" "$NAUTY/pickg" "$NAUTY/planarg" \
  > evidence/programs.sha256
```

The following shell function generates only the biconnected, nonplanar graphs
from a Minibaum enumeration, in 256 deterministic residue classes.  The
unfiltered graphs are passed directly between the three programs and are not
stored.

```bash
generate_biconnected_nonplanar() {
  local n=$1
  local girth=$2
  local output=$3
  local girth_mode=${4:-minimum}
  local pickg_args=(-q -c2)
  if [ "$girth_mode" = exact ]; then
    pickg_args+=("-g$girth")
  fi
  mkdir -p "$output"
  for part in $(seq 0 $((PARTS - 1))); do
    local file
    file=$(printf '%s/part-%03d.g6' "$output" "$part")
    "$MINIBAUM" "$n" "$girth" s g m "$part" "$PARTS" |
      "$NAUTY/pickg" "${pickg_args[@]}" |
      "$NAUTY/planarg" -q -v > "$file"
  done
}
```

The program `pickg -c2` tests biconnectivity, while `planarg -v` tests
nonplanarity.  Neither program performs both tests, so both filters are
needed.  Passing `exact` as the fourth argument additionally applies
`pickg -gN`, retaining only graphs of girth exactly `N`.

Before the graph families are processed, verify the finite local drawings
used in the 4-cycle and 5-cycle expansion lemmas:

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

Generate the biconnected, nonplanar cubic graphs with $n=24$ and girth at
least 5, then verify 2-flexibility:

```bash
generate_biconnected_nonplanar 24 5 data/claim2-biconnected-nonplanar

wc -l data/claim2-biconnected-nonplanar/*.g6
claim2_count=$(wc -l data/claim2-biconnected-nonplanar/*.g6 | awk 'END {print $1}')
test "$claim2_count" -eq 1620470

mkdir -p evidence/claim2
for input in data/claim2-biconnected-nonplanar/*.g6; do
  name=$(basename "$input" .g6)
  ./oops -i="$input" -verify-cubic -colors=0 \
    > "evidence/claim2/${name}.log" 2>&1
done
```

The family contains 1,620,470 records.  Every record must be reported as
2-flexible.

## Verify Claim 3

Generate the order-26 family once.  Girth-five graphs are mapped to their
lexicographically least canonical 5-cycle contraction, while graphs of
larger girth pass through unchanged.  Sorting discards equal order-22 cores:

```bash
generate_biconnected_nonplanar 26 5 data/claim3-biconnected-nonplanar
wc -l data/claim3-biconnected-nonplanar/*.g6

for input in data/claim3-biconnected-nonplanar/*.g6; do
  ./prepare_cubic claim3 < "$input"
done | sort -u > data/claim3.g6
wc -l data/claim3.g6

mkdir -p evidence/claim3
./oops -i=data/claim3.g6 -verify-cubic -colors=0 \
  > evidence/claim3/claim3.log 2>&1
```

Every retained core must be reported as a 5-cycle core, and every direct
order-26 record must be reported as 1-flexible.  The complete preparation
records the exact counts; bounded samples project about 2.8 million cores.
The connected input family contains 31,478,584 graphs before the
biconnectivity and nonplanarity filters.

## Verify Claim 4

Generate every biconnected, nonplanar cubic graph with $n=28$ and girth
exactly 5.  Each parent emits one canonical pair `(H,b)`.  Sorting removes
duplicate pairs; the grouping pass combines all pairs with the same `H` into
one graph carrying one pendant marker at each required vertex:

```bash
generate_biconnected_nonplanar 28 5 data/claim4-biconnected-nonplanar exact
wc -l data/claim4-biconnected-nonplanar/*.g6

for input in data/claim4-biconnected-nonplanar/*.g6; do
  ./prepare_cubic claim4 < "$input"
done |
  sort -u |
  ./prepare_cubic group-claim4 > data/claim4-marked-groups.g6
wc -l data/claim4-marked-groups.g6

mkdir -p evidence/claim4
./oops -i=data/claim4-marked-groups.g6 -verify-cubic -colors=0 \
  > evidence/claim4/claim4.log 2>&1
```

Every record must be reported as a marked-star group.  The verifier removes
the pendant markers, builds one SAT instance for `H`, and checks the
three-edge star at every marked vertex.  The complete preparation records the
exact number of groups; bounded forward and reverse samples project about
29 million.  The connected exact-girth-five family has 652,159,389 records
before the biconnectivity and nonplanarity filters.

## Verify Claim 5

Generate every biconnected, nonplanar cubic graph with $n=28$ and girth at
least 6, then process the residue files:

```bash
generate_biconnected_nonplanar 28 6 data/claim5-biconnected-nonplanar
wc -l data/claim5-biconnected-nonplanar/*.g6

mkdir -p evidence/claim5
for input in data/claim5-biconnected-nonplanar/*.g6; do
  name=$(basename "$input" .g6)
  ./oops -i="$input" -sat=1 -unsat=0 -colors=0 \
    > "evidence/claim5/${name}.log" 2>&1
done
```

The family contains 4,624,501 records.  The completed verification reported
4,624,501 1-planar graphs and no other verdict; every invocation exited
successfully.  Its wall-clock time was approximately three days when the 256
independent residue files were distributed over 48 workers.

## Final checks

Every OOPS invocation must exit successfully.  Sum `#planar` and
`#1-planar` over the logs and compare the result with the appropriate family
size above.  Finally, record hashes of all inputs and logs:

```bash
find data evidence -type f ! -name files.sha256 -print0 | sort -z | xargs -0 sha256sum \
  > evidence/files.sha256
```
