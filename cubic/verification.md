# Reproducing the five computational claims

This file contains the commands used for the five claims in
`main_cubic.tex`.  The paper contains the proofs.

Run the commands in one Bash session.  They are written serially.  A claim
passes only if every command succeeds and all generated graphs are checked.

## Prerequisites

The verification uses the following software:

- A C++17 compiler, a C compiler, Python 3, GNU Make, Bash, and standard Unix
  utilities including `awk`, `sha256sum`, `find`, `sort`, `xargs`, and `wc`.
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
export LC_ALL=C

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

for input in data/claim3-biconnected-nonplanar/*.g6; do
  ./prepare_cubic claim3 < "$input"
done | sort -u > data/claim3.g6
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

- `D parent`, when Claim 2 completes the proof;
- `B parent code`, for a smaller graph with one degree-6 vertex;
- `C parent code`, for a smaller graph with one degree-7 vertex; or
- `R parent code`, for a 5-cycle-to-path reduction whose marked three-edge star
  must remain uncrossed.

Every input graph must produce exactly one output line.

```bash
generate_biconnected_nonplanar 28 5 data/claim4-biconnected-nonplanar exact
mkdir -p evidence/claim4
claim4_parents=$(
  wc -l data/claim4-biconnected-nonplanar/*.g6 |
    awk 'END {print $1}'
)

for input in data/claim4-biconnected-nonplanar/*.g6; do
  ./prepare_cubic claim4-joint < "$input"
done > data/claim4-classified.txt 2> evidence/claim4/preparation.log

test "$(wc -l < data/claim4-classified.txt)" -eq "$claim4_parents"
awk '
  ($1 == "D" && NF == 2) ||
  (($1 == "B" || $1 == "C" || $1 == "R") && NF == 3) { next }
  { exit 1 }
' data/claim4-classified.txt
awk '{print $2}' data/claim4-classified.txt |
  cmp - <(cat data/claim4-biconnected-nonplanar/*.g6)

awk 'NF == 3 {print $3}' data/claim4-classified.txt |
  sort -u > data/claim4.g6

for class in D B C R; do
  awk -v c="$class" '$1 == c {count++} END {print c, count + 0}' \
    data/claim4-classified.txt
done
wc -l data/claim4.g6
```

The four class counts must add up to `claim4_parents`.

## Verify Claim 3

```bash
mkdir -p evidence/claim3
./oops -i=data/claim3.g6 -verify-cubic -timeout=120 \
  -Cverify-cubic-residue=evidence/claim3/claim3b-input.g6 -colors=0 \
  > evidence/claim3/claim3.log 2>&1

if [ -s evidence/claim3/claim3b-input.g6 ]; then
  ./oops -i=evidence/claim3/claim3b-input.g6 \
    -verify-cubic -Cclaim3b -timeout=120 \
    -Cverify-cubic-residue=evidence/claim3/claim3b-residue.g6 \
    -colors=0 > evidence/claim3/claim3b.log 2>&1
else
  : > evidence/claim3/claim3b-residue.g6
fi
test ! -s evidence/claim3/claim3b-residue.g6
```

The first pass retains records that reach the 120-second limit.  Claim 3b
expands each retained order-22 Claim 3 core into at most 12 relevant
order-26 parents and verifies them directly; it retries a retained order-26
Claim 3 graph directly.  Accept Claim 3 only if the two logs account for
every input and parent and
`claim3b-residue.g6` is empty.  When several OOPS processes are used, give
each process a different residue filename.

In the first log, `#unknown` must equal the number of lines in
`claim3b-input.g6`.  In the Claim 3b log, `#claim3b-records` and
`#claim3b-completed` must both equal that number, and `#unknown` must be zero.

## Verify Claim 4

```bash
./oops -i=data/claim4.g6 -verify-cubic -timeout=60 \
  -Cverify-cubic-residue=evidence/claim4/residue.g6 -colors=0 \
  > evidence/claim4/verification.log 2>&1
```

A B record requires an uncrossed degree-6 hub, a C record requires a
permitted drawing of an uncrossed degree-7 hub, and an R record requires the
marked three-edge star to be uncrossed.  The D records need no solver input
because Claim 2 already proves them.

Verify the parents represented by records that reached the time limit:

```bash
awk '
  FILENAME == ARGV[1] { residue[$1] = 1; next }
  NF == 3 && ($3 in residue) { print $2; found[$3] = 1 }
  END { for (code in residue) if (!(code in found)) exit 1 }
' evidence/claim4/residue.g6 data/claim4-classified.txt |
  sort -u > data/claim4-residue-parents.g6

./oops -i=data/claim4-residue-parents.g6 -sat=1 -unsat=0 \
  -timeout=900 -colors=0 > evidence/claim4/residue-parents.log 2>&1
```

Accept Claim 4 only if:

- in `verification.log`, `#planar + #1-planar + #unknown` equals the number
  of lines in `data/claim4.g6`, and `#unknown` equals the number of lines in
  `residue.g6`;
- in `residue-parents.log`, `#planar + #1-planar` equals the number of lines
  in `data/claim4-residue-parents.g6`, while `#non-1-planar` and `#unknown`
  are both zero.

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
