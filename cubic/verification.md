# Reproducing the five computational claims

This script verifies that every cubic graph with $n\le 28$ is 1-planar.
It uses 48 parallel jobs; change `JOBS` below for a different machine.

## Prerequisites

- C and C++17 compilers, Python 3, GNU Make, Bash, zlib, and standard GNU utilities;
- [nauty and Traces 2.9.3](https://users.cecs.anu.edu.au/~bdm/nauty/nauty2_9_3.tar.gz)
  (`pickg` and `planarg`);
- [Minibaum 5](https://caagt.ugent.be/minibaum/minibaum5.c).

## Preparation

Set the source paths, build the executables, and choose the run parameters:

```bash
set -euo pipefail

NAUTY=/path/to/nauty2_9_3
MINIBAUM_SOURCE=/path/to/minibaum5.c
MINIBAUM=./minibaum5
SHARDS=288
JOBS=48
export LC_ALL=C
export NAUTY MINIBAUM SHARDS

# Build nauty.
(cd "$NAUTY" && ./configure && make)

# Build Minibaum.
gcc -std=gnu89 -O3 "$MINIBAUM_SOURCE" -o "$MINIBAUM"

# Build OOPS.
make

# Build the cubic-graph preparation tool.
g++ -std=c++17 -O3 -Isrc -I"$NAUTY" \
  cubic/prepare_cubic.cpp "$NAUTY/nauty.a" \
  -o prepare_cubic

# Build the drawing-table checker.
g++ -std=c++17 -O2 -Isrc \
  cubic/test_cycle_expansions.cpp src/planarity_test.cpp \
  -o test_cycle_expansions

mkdir verification
```

## Utilities

```bash
# Run a command once for each NUL-separated input.
run_in_parallel() {
  local command=$1
  xargs -0 -r -n1 -P "$JOBS" bash -c "
    set -euo pipefail
    $command
  " bash
}

# Generate biconnected, nonplanar cubic graphs of minimum or exact girth.
generate_minibaum_parts() {
  local n=$1
  local girth=$2
  local output=$3
  local girth_mode=${4:-minimum}
  export n girth output girth_mode
  mkdir -p "$output"
  seq 0 $((SHARDS - 1)) | tr '\n' '\0' | run_in_parallel '
    part=$1
    file=$(printf "%s/n%02d-part-%03d.g6" "$output" "$n" "$part")
    if [ "$girth_mode" = exact ]; then
      "$MINIBAUM" "$n" "$girth" s g m "$part" "$SHARDS" |
        "$NAUTY/pickg" -q -c2 "-g$girth" |
        "$NAUTY/planarg" -q -v > "$file"
    else
      "$MINIBAUM" "$n" "$girth" s g m "$part" "$SHARDS" |
        "$NAUTY/pickg" -q -c2 |
        "$NAUTY/planarg" -q -v > "$file"
    fi
  '
}

# Shuffle an input file and split it into balanced shards.
make_shards() {
  local input=$1
  local output=$2
  mkdir "$output"

  shuf "$input" > "$output/input"
  split -n "l/$SHARDS" -d -a 3 --additional-suffix=.g6 \
    "$output/input" "$output/shard-"
  cat "$output"/shard-*.g6 | cmp - "$output/input"
}

# Sum a counter over the final summary line of each log.
sum_counter() {
  local counter=$1
  shift
  for log in "$@"; do
    grep -F "$counter =" "$log" | tail -n 1
  done | tr -d ',' | awk -v counter="$counter" '
    { for (i = 1; i <= NF; i++) if ($i == counter) total += $(i + 2) }
    END { print total + 0 }
  '
}
```

## Verify the drawing tables

```bash
./test_cycle_expansions
python3 cubic/test_hub_expansions.py
```

## Verify Claim 1

Generate every biconnected, nonplanar, triangle-free cubic graph with
$n\le 22$:

```bash
mkdir verification/claim1

for n in 4 6 8 10 12 14 16 18 20 22; do
  generate_minibaum_parts "$n" 4 verification/claim1
done

cat verification/claim1/*.g6 > verification/claim1.g6
test "$(wc -l < verification/claim1.g6)" -eq 1538493

mkdir verification/claim1/logs
make_shards verification/claim1.g6 verification/claim1/shards
printf '%s\0' verification/claim1/shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -colors=0 \
      > "verification/claim1/logs/${name}.log" 2>&1
'

test "$(sum_counter '#3-flexible' verification/claim1/logs/*.log)" -eq 1538493
test "$(sum_counter '#unknown' verification/claim1/logs/*.log)" -eq 0
```

## Verify Claim 2

Generate the biconnected, nonplanar cubic graphs with $n=24$ and girth at
least 5, then verify 2-flexibility:

```bash
generate_minibaum_parts 24 5 verification/claim2-biconnected-nonplanar
cat verification/claim2-biconnected-nonplanar/*.g6 > verification/claim2.g6
test "$(wc -l < verification/claim2.g6)" -eq 1620470

mkdir verification/claim2
mkdir verification/claim2/logs
make_shards verification/claim2.g6 verification/claim2/shards
printf '%s\0' verification/claim2/shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -colors=0 \
      > "verification/claim2/logs/${name}.log" 2>&1
  '

test "$(sum_counter '#2-flexible' verification/claim2/logs/*.log)" -eq 1620470
test "$(sum_counter '#unknown' verification/claim2/logs/*.log)" -eq 0
```

## Verify Claim 3

### Prepare data

Generate and deduplicate the Claim 3 records.

```bash
generate_minibaum_parts 26 5 verification/claim3-biconnected-nonplanar
test "$(wc -l verification/claim3-biconnected-nonplanar/*.g6 | awk 'END {print $1}')" \
  -eq 31478465

mkdir verification/claim3-prepared
printf '%s\0' verification/claim3-biconnected-nonplanar/*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./prepare_cubic claim3 < "$input" > "verification/claim3-prepared/$name.g6"
  '
cat verification/claim3-prepared/*.g6 |
  sort --parallel="$JOBS" -u > verification/claim3.g6
test "$(wc -l < verification/claim3.g6)" -eq 4881876
```

### Verify records

```bash
mkdir verification/claim3
mkdir verification/claim3/pass1-logs verification/claim3/pass1-residues
make_shards verification/claim3.g6 verification/claim3/pass1-shards

printf '%s\0' verification/claim3/pass1-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -timeout=120 \
      -Cverify-cubic-residue="verification/claim3/pass1-residues/$name.g6" \
      -colors=0 > "verification/claim3/pass1-logs/$name.log" 2>&1
  '

cat verification/claim3/pass1-residues/*.g6 | sort -u > verification/claim3/claim3b-input.g6

mkdir verification/claim3/claim3b-logs verification/claim3/claim3b-residues
make_shards verification/claim3/claim3b-input.g6 verification/claim3/claim3b-shards
printf '%s\0' verification/claim3/claim3b-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -Cclaim3b -timeout=120 \
      -Cverify-cubic-residue="verification/claim3/claim3b-residues/$name.g6" \
      -colors=0 > "verification/claim3/claim3b-logs/$name.log" 2>&1
  '

cat verification/claim3/claim3b-residues/*.g6 | sort -u > verification/claim3/claim3b-residue.g6
test ! -s verification/claim3/claim3b-residue.g6

claim3_planar=$(sum_counter '#planar' verification/claim3/pass1-logs/*.log)
claim3_1planar=$(sum_counter '#1-planar' verification/claim3/pass1-logs/*.log)
claim3_unknown=$(sum_counter '#unknown' verification/claim3/pass1-logs/*.log)
test "$((claim3_planar + claim3_1planar + claim3_unknown))" -eq 4881876
test "$claim3_unknown" -eq "$(wc -l < verification/claim3/claim3b-input.g6)"

claim3b_records=$(wc -l < verification/claim3/claim3b-input.g6)
claim3b_parents=$(sum_counter '#claim3b-parents' verification/claim3/claim3b-logs/*.log)
test "$(sum_counter '#claim3b-records' verification/claim3/claim3b-logs/*.log)" \
  -eq "$claim3b_records"
test "$(sum_counter '#claim3b-completed' verification/claim3/claim3b-logs/*.log)" \
  -eq "$claim3b_records"
test "$(sum_counter '#1-flexible' verification/claim3/claim3b-logs/*.log)" \
  -eq "$claim3b_parents"
test "$(sum_counter '#unknown' verification/claim3/claim3b-logs/*.log)" -eq 0
```

## Verify Claim 4

### Prepare data

Generate and classify the Claim 4 graphs.  Each output line is `D` or a
`B`, `C`, or `R` record followed by its code.

```bash
generate_minibaum_parts 28 5 verification/claim4-biconnected-nonplanar exact
mkdir verification/claim4
mkdir verification/claim4-classified verification/claim4/preparation

printf '%s\0' verification/claim4-biconnected-nonplanar/*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    classified="verification/claim4-classified/$name.txt"
    ./prepare_cubic claim4-joint < "$input" \
      > "$classified" \
      2> "verification/claim4/preparation/$name.log"
    test "$(wc -l < "$input")" -eq "$(wc -l < "$classified")"
  '

test "$(wc -l verification/claim4-biconnected-nonplanar/*.g6 | awk 'END {print $1}')" \
  -eq 652157758
awk '
  { count[$1]++; if (NF == 2) print $2 }
  END { exit count["D"] != 254354420 || count["B"] != 350383009 ||
             count["C"] != 45970275 || count["R"] != 1450054 }
' verification/claim4-classified/*.txt |
  sort --parallel="$JOBS" -u > verification/claim4.g6
test "$(wc -l < verification/claim4.g6)" -eq 3788793
```

### Verify records

```bash
mkdir verification/claim4/certificate-logs \
      verification/claim4/certificate-residues
make_shards verification/claim4.g6 verification/claim4/certificate-shards

printf '%s\0' verification/claim4/certificate-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -timeout=60 \
      -Cverify-cubic-residue="verification/claim4/certificate-residues/$name.g6" \
      -colors=0 > "verification/claim4/certificate-logs/$name.log" 2>&1
  '

cat verification/claim4/certificate-residues/*.g6 | sort -u > verification/claim4/residue.g6

claim4_planar=$(sum_counter '#planar' verification/claim4/certificate-logs/*.log)
claim4_1planar=$(sum_counter '#1-planar' verification/claim4/certificate-logs/*.log)
claim4_unknown=$(sum_counter '#unknown' verification/claim4/certificate-logs/*.log)
test "$((claim4_planar + claim4_1planar + claim4_unknown))" -eq 3788793
test "$claim4_unknown" -eq "$(wc -l < verification/claim4/residue.g6)"
```

### Verify residue parents

```bash
paste <(cat verification/claim4-biconnected-nonplanar/*.g6) \
      <(cat verification/claim4-classified/*.txt) |
awk '
  FILENAME == ARGV[1] { residue[$1] = 1; next }
  NF == 3 && ($3 in residue) { print $1; found[$3] = 1 }
  END { for (code in residue) if (!(code in found)) exit 1 }
' verification/claim4/residue.g6 - |
sort --parallel="$JOBS" -u > verification/claim4-residue-parents.g6

mkdir verification/claim4/parent-logs
make_shards verification/claim4-residue-parents.g6 verification/claim4/parent-shards

printf '%s\0' verification/claim4/parent-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -colors=0 \
      > "verification/claim4/parent-logs/$name.log" 2>&1
  '

claim4_parent_count=$(wc -l < verification/claim4-residue-parents.g6)
test "$(sum_counter '#1-planar' verification/claim4/parent-logs/*.log)" \
  -eq "$claim4_parent_count"
test "$(sum_counter '#non-1-planar' verification/claim4/parent-logs/*.log)" -eq 0
test "$(sum_counter '#unknown' verification/claim4/parent-logs/*.log)" -eq 0
```

## Verify Claim 5

Generate every biconnected, nonplanar cubic graph with $n=28$ and girth at
least 6, then verify it directly:

```bash
generate_minibaum_parts 28 6 verification/claim5-biconnected-nonplanar
cat verification/claim5-biconnected-nonplanar/*.g6 > verification/claim5.g6
test "$(wc -l < verification/claim5.g6)" -eq 4624501

mkdir verification/claim5
mkdir verification/claim5/logs
make_shards verification/claim5.g6 verification/claim5/shards
printf '%s\0' verification/claim5/shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -colors=0 \
      > "verification/claim5/logs/${name}.log" 2>&1
  '

test "$(sum_counter '#1-planar' verification/claim5/logs/*.log)" -eq 4624501
test "$(sum_counter '#non-1-planar' verification/claim5/logs/*.log)" -eq 0
test "$(sum_counter '#unknown' verification/claim5/logs/*.log)" -eq 0
```
