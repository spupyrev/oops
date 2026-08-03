# Reproducing the five computational claims

This file contains the commands used for the five claims in
`main_cubic.tex`.  The paper contains the proofs.

Run the commands in one Bash session.  `JOBS` defaults to the number of
available cores and may be overridden.  A claim passes only if every command
succeeds and all generated graphs are checked.

## Prerequisites

The verification uses the following software:

- A C++17 compiler, a C compiler, Python 3, GNU Make, Bash, and standard Unix
  utilities including `awk`, `cmp`, `find`, `nproc`, `paste`, `sha256sum`,
  `sort`, `split`, `xargs`, and `wc`.
  Building OOPS also requires the zlib development headers and library.
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

make
gcc -std=gnu89 -O3 /path/to/minibaum5.c -o minibaum5

NAUTY=/path/to/nauty2_9_3
MINIBAUM=./minibaum5
PARTS=256
JOBS=${JOBS:-$(nproc)}
SHARDS=${SHARDS:-256}
export LC_ALL=C
export NAUTY MINIBAUM PARTS

g++ -std=c++17 -O3 -Isrc -I"$NAUTY" \
  cubic/prepare_cubic.cpp "$NAUTY/nauty.a" \
  -o prepare_cubic

mkdir -p data evidence
sha256sum ./oops ./prepare_cubic "$MINIBAUM" /path/to/minibaum5.c \
  "$NAUTY/geng" "$NAUTY/pickg" "$NAUTY/planarg" \
  cubic/test_cycle_expansions.cpp cubic/test_hub_expansions.py \
  > evidence/programs.sha256
```

This function divides Minibaum generation into 256 fixed parts.  It keeps
only biconnected, nonplanar graphs.

```bash
generate_biconnected_nonplanar() {
  local n=$1
  local girth=$2
  local output=$3
  local girth_mode=${4:-minimum}
  mkdir -p "$output"
  seq 0 $((PARTS - 1)) |
    xargs -n1 -P "$JOBS" bash -c '
      set -euo pipefail
      n=$1
      girth=$2
      output=$3
      girth_mode=$4
      part=$5
      file=$(printf "%s/part-%03d.g6" "$output" "$part")
      if [ "$girth_mode" = exact ]; then
        "$MINIBAUM" "$n" "$girth" s g m "$part" "$PARTS" |
          "$NAUTY/pickg" -q -c2 "-g$girth" |
          "$NAUTY/planarg" -q -v > "$file"
      else
        "$MINIBAUM" "$n" "$girth" s g m "$part" "$PARTS" |
          "$NAUTY/pickg" -q -c2 |
          "$NAUTY/planarg" -q -v > "$file"
      fi
    ' bash "$n" "$girth" "$output" "$girth_mode"
}

make_shards() {
  local input=$1
  local output=$2
  mkdir "$output"
  split -n "l/$SHARDS" -d -a 3 --additional-suffix=.g6 \
    "$input" "$output/shard-"
  cat "$output"/*.g6 | cmp - "$input"
}
```

`pickg -c2` keeps biconnected graphs.  `planarg -v` keeps nonplanar graphs.
The optional argument `exact` also removes graphs whose girth is larger than
the requested value.

Check the small drawing tables used in the paper:

```bash
g++ -std=c++17 -O2 -Isrc \
  cubic/test_cycle_expansions.cpp src/planarity_test.cpp \
  -o test_cycle_expansions
./test_cycle_expansions
python3 cubic/test_hub_expansions.py
```

The last command must report 60 expandable `BALOM` orders and 355 expandable
`BALOOM` orders.

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

Confirm that the counts, in the same order, are

```text
1, 2, 5, 19, 85, 509, 4,060, 41,301, 97,141, 1,432,712.
```

Their sum must be 1,575,835.  Accept the result only if every nonplanar graph
is reported as 3-flexible.  Planar graphs require no further check.

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

Confirm that the total is 1,620,470.  Accept the result only if every graph
is reported as 2-flexible.

## Prepare Claim 3

Generate the order-26 graphs.  For a graph of girth five,
`prepare_cubic` replaces one 5-cycle by a degree-five vertex.  Equal results
are removed.  Graphs of larger girth are kept unchanged.

```bash
generate_biconnected_nonplanar 26 5 data/claim3-biconnected-nonplanar
wc -l data/claim3-biconnected-nonplanar/*.g6

mkdir data/claim3-prepared
find data/claim3-biconnected-nonplanar -type f -name '*.g6' -print0 |
  sort -z | xargs -0 -n1 -P "$JOBS" bash -c '
    set -euo pipefail
    input=$1
    name=$(basename "$input" .g6)
    ./prepare_cubic claim3 < "$input" > "data/claim3-prepared/$name.g6"
  ' bash
cat data/claim3-prepared/*.g6 | sort -u > data/claim3.g6
wc -l data/claim3.g6

awk 'substr($0, 1, 1) == "U"' data/claim3.g6 \
  > data/claim3-cores.g6
awk 'substr($0, 1, 1) == "Y"' data/claim3.g6 \
  > data/claim3-girth6.g6
```

Confirm that the filtered input contains 31,478,465 graphs: 31,297,238 of
girth five and 181,227 of larger girth.  (The unfiltered connected family has
31,478,584 graphs.)  Record the number left after removing duplicates; this
is the output of `wc`.

## Prepare Claim 4

Generate the order-28 graphs of girth five.  For each graph,
`prepare_cubic` prints one of:

- `D`, when Claim 2 completes the proof;
- `B code`, for a smaller graph with one degree-6 vertex;
- `C code`, for a smaller graph with one degree-7 vertex; or
- `R code`, for a 5-cycle-to-path reduction whose marked three-edge star
  must remain uncrossed.

Every input graph must produce exactly one output line.

```bash
generate_biconnected_nonplanar 28 5 data/claim4-biconnected-nonplanar exact
mkdir -p evidence/claim4
mkdir data/claim4-classified evidence/claim4/preparation

find data/claim4-biconnected-nonplanar -type f -name '*.g6' -print0 |
  sort -z | xargs -0 -n1 -P "$JOBS" bash -c '
    set -euo pipefail
    input=$1
    name=$(basename "$input" .g6)
    ./prepare_cubic claim4-joint < "$input" \
      > "data/claim4-classified/$name.txt" \
      2> "evidence/claim4/preparation/$name.log"
  ' bash

claim4_parents=$(cat data/claim4-biconnected-nonplanar/*.g6 | wc -l)
test "$claim4_parents" -eq 652157758
for input in data/claim4-biconnected-nonplanar/*.g6; do
  name=$(basename "$input" .g6)
  test "$(wc -l < "$input")" -eq \
       "$(wc -l < "data/claim4-classified/$name.txt")"
done

awk '
  ($1 == "D" && NF == 1) ||
  (($1 == "B" || $1 == "C" || $1 == "R") && NF == 2) {
    count[$1]++
    next
  }
  { exit 1 }
  END {
    print "D", count["D"] + 0
    print "B", count["B"] + 0
    print "C", count["C"] + 0
    print "R", count["R"] + 0
  }
' data/claim4-classified/*.txt

awk 'NF == 2 {print $2}' data/claim4-classified/*.txt |
  sort -u > data/claim4.g6
wc -l data/claim4.g6
```

The class counts must be 254,354,420 `D`, 350,383,009 `B`, 45,970,275
`C`, and 1,450,054 `R`; they sum to `claim4_parents`.  There must be
3,788,793 distinct records in `data/claim4.g6`.  Keep the classification
files: they are the expensive artifact, whereas the parent files can be
regenerated quickly.

## Verify Claim 3

```bash
mkdir -p evidence/claim3
mkdir evidence/claim3/pass1-logs evidence/claim3/pass1-residues
make_shards data/claim3.g6 evidence/claim3/pass1-shards

find evidence/claim3/pass1-shards -type f -name '*.g6' -size +0c -print0 |
  sort -z | xargs -0 -n1 -P "$JOBS" bash -c '
    set -euo pipefail
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -timeout=120 \
      -Cverify-cubic-residue="evidence/claim3/pass1-residues/$name.g6" \
      -colors=0 > "evidence/claim3/pass1-logs/$name.log" 2>&1
  ' bash

cat evidence/claim3/pass1-residues/*.g6 |
  sort -u > evidence/claim3/claim3b-input.g6
test "$(wc -l < evidence/claim3/claim3b-input.g6)" -eq 705

mkdir evidence/claim3/claim3b-logs evidence/claim3/claim3b-residues
make_shards evidence/claim3/claim3b-input.g6 evidence/claim3/claim3b-shards
find evidence/claim3/claim3b-shards -type f -name '*.g6' -size +0c -print0 |
  sort -z | xargs -0 -n1 -P "$JOBS" bash -c '
    set -euo pipefail
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -Cclaim3b -timeout=120 \
      -Cverify-cubic-residue="evidence/claim3/claim3b-residues/$name.g6" \
      -colors=0 > "evidence/claim3/claim3b-logs/$name.log" 2>&1
  ' bash

cat evidence/claim3/claim3b-residues/*.g6 |
  sort -u > evidence/claim3/claim3b-residue.g6
test ! -s evidence/claim3/claim3b-residue.g6
```

The first-pass logs must account for all 4,881,876 records and report 705
unknown.  Claim 3b expands these records into 8,460 parents.  Its logs must
report 705 completed records, 8,460 1-flexible parents, and no unknown.

## Verify Claim 4

```bash
mkdir evidence/claim4/certificate-logs \
      evidence/claim4/certificate-residues
make_shards data/claim4.g6 evidence/claim4/certificate-shards

find evidence/claim4/certificate-shards -type f -name '*.g6' -print0 |
  sort -z | xargs -0 -n1 -P "$JOBS" bash -c '
    set -euo pipefail
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -timeout=60 \
      -Cverify-cubic-residue="evidence/claim4/certificate-residues/$name.g6" \
      -colors=0 > "evidence/claim4/certificate-logs/$name.log" 2>&1
  ' bash

cat evidence/claim4/certificate-residues/*.g6 |
  sort -u > evidence/claim4/residue.g6
```

A B record requires an uncrossed degree-6 hub, a C record requires a
permitted drawing of an uncrossed degree-7 hub, and an R record requires the
marked three-edge star to be uncrossed.  The D records need no solver input
because Claim 2 already proves them.

Recover and verify the parents represented by records that reached the time
limit:

```bash
mkdir evidence/claim4/residue-parent-records
for input in data/claim4-biconnected-nonplanar/*.g6; do
  name=$(basename "$input" .g6)
  paste "$input" "data/claim4-classified/$name.txt" |
    awk '
      FILENAME == ARGV[1] { residue[$1] = 1; next }
      NF == 3 && ($3 in residue) { print $1, $3 }
    ' evidence/claim4/residue.g6 - \
    > "evidence/claim4/residue-parent-records/$name.txt"
done

awk '
  FILENAME == ARGV[1] { residue[$1] = 1; next }
  { found[$2] = 1 }
  END { for (code in residue) if (!(code in found)) exit 1 }
' evidence/claim4/residue.g6 \
  evidence/claim4/residue-parent-records/*.txt

awk '{print $1}' evidence/claim4/residue-parent-records/*.txt |
  sort -u > data/claim4-residue-parents.g6

mkdir evidence/claim4/parent-logs
make_shards data/claim4-residue-parents.g6 evidence/claim4/parent-shards

find evidence/claim4/parent-shards -type f -name '*.g6' -print0 |
  sort -z | xargs -0 -n1 -P "$JOBS" bash -c '
    set -euo pipefail
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -sat=1 -unsat=0 -timeout=900 -colors=0 \
      > "evidence/claim4/parent-logs/$name.log" 2>&1
  ' bash
```

Accept Claim 4 only if:

- the counters in `certificate-logs` account for all 3,788,793 records, and
  their total `#unknown` equals the number of lines in `residue.g6`;
- the counters in `parent-logs` account for every line of
  `data/claim4-residue-parents.g6`, with `#non-1-planar = 0` and
  `#unknown = 0`.  These parents are nonplanar by construction, so their
  total `#planar` is also zero.

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

Confirm that the family contains 4,624,501 graphs.  Accept the result only if
all 4,624,501 graphs are reported as 1-planar and every command succeeds.

## Final checks

Every OOPS invocation must exit successfully.  Sum `#planar` and
`#1-planar` over the logs and compare the result with the appropriate family
size above.  Finally, record hashes of all inputs and logs:

```bash
find data evidence -type f ! -name files.sha256 -print0 | sort -z | xargs -0 sha256sum \
  > evidence/files.sha256
```
