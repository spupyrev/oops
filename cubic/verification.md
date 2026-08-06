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
PARTS=256
JOBS=48
export LC_ALL=C
export NAUTY MINIBAUM PARTS

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

mkdir -p data evidence
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
  seq 0 $((PARTS - 1)) | tr '\n' '\0' | run_in_parallel '
    part=$1
    file=$(printf "%s/n%02d-part-%03d.g6" "$output" "$n" "$part")
    if [ "$girth_mode" = exact ]; then
      "$MINIBAUM" "$n" "$girth" s g m "$part" "$PARTS" |
        "$NAUTY/pickg" -q -c2 "-g$girth" |
        "$NAUTY/planarg" -q -v > "$file"
    else
      "$MINIBAUM" "$n" "$girth" s g m "$part" "$PARTS" |
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
  split -n "l/$PARTS" -d -a 3 --additional-suffix=.g6 \
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
mkdir data/claim1 evidence/claim1

for n in 4 6 8 10 12 14 16 18 20 22; do
  generate_minibaum_parts "$n" 4 data/claim1
done

cat data/claim1/*.g6 > data/claim1.g6
test "$(wc -l < data/claim1.g6)" -eq 1538493

mkdir evidence/claim1/logs
make_shards data/claim1.g6 evidence/claim1/shards
printf '%s\0' evidence/claim1/shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -colors=0 \
      > "evidence/claim1/logs/${name}.log" 2>&1
'

test "$(sum_counter '#3-flexible' evidence/claim1/logs/*.log)" -eq 1538493
test "$(sum_counter '#unknown' evidence/claim1/logs/*.log)" -eq 0
```

## Verify Claim 2

Generate the biconnected, nonplanar cubic graphs with $n=24$ and girth at
least 5, then verify 2-flexibility:

```bash
generate_minibaum_parts 24 5 data/claim2-biconnected-nonplanar
cat data/claim2-biconnected-nonplanar/*.g6 > data/claim2.g6
test "$(wc -l < data/claim2.g6)" -eq 1620470

mkdir evidence/claim2
mkdir evidence/claim2/logs
make_shards data/claim2.g6 evidence/claim2/shards
printf '%s\0' evidence/claim2/shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -colors=0 \
      > "evidence/claim2/logs/${name}.log" 2>&1
  '

test "$(sum_counter '#2-flexible' evidence/claim2/logs/*.log)" -eq 1620470
test "$(sum_counter '#unknown' evidence/claim2/logs/*.log)" -eq 0
```

## Verify Claim 3

### Prepare data

Generate and deduplicate the Claim 3 records.

```bash
generate_minibaum_parts 26 5 data/claim3-biconnected-nonplanar
test "$(wc -l data/claim3-biconnected-nonplanar/*.g6 | awk 'END {print $1}')" \
  -eq 31478465

mkdir data/claim3-prepared
printf '%s\0' data/claim3-biconnected-nonplanar/*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./prepare_cubic claim3 < "$input" > "data/claim3-prepared/$name.g6"
  '
cat data/claim3-prepared/*.g6 |
  sort --parallel="$JOBS" -u > data/claim3.g6
test "$(wc -l < data/claim3.g6)" -eq 4881876
```

### Verify records

```bash
mkdir evidence/claim3
mkdir evidence/claim3/pass1-logs evidence/claim3/pass1-residues
make_shards data/claim3.g6 evidence/claim3/pass1-shards

printf '%s\0' evidence/claim3/pass1-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -timeout=120 \
      -Cverify-cubic-residue="evidence/claim3/pass1-residues/$name.g6" \
      -colors=0 > "evidence/claim3/pass1-logs/$name.log" 2>&1
  '

cat evidence/claim3/pass1-residues/*.g6 | sort -u > evidence/claim3/claim3b-input.g6

mkdir evidence/claim3/claim3b-logs evidence/claim3/claim3b-residues
make_shards evidence/claim3/claim3b-input.g6 evidence/claim3/claim3b-shards
printf '%s\0' evidence/claim3/claim3b-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -Cclaim3b -timeout=120 \
      -Cverify-cubic-residue="evidence/claim3/claim3b-residues/$name.g6" \
      -colors=0 > "evidence/claim3/claim3b-logs/$name.log" 2>&1
  '

cat evidence/claim3/claim3b-residues/*.g6 | sort -u > evidence/claim3/claim3b-residue.g6
test ! -s evidence/claim3/claim3b-residue.g6

claim3_planar=$(sum_counter '#planar' evidence/claim3/pass1-logs/*.log)
claim3_1planar=$(sum_counter '#1-planar' evidence/claim3/pass1-logs/*.log)
claim3_unknown=$(sum_counter '#unknown' evidence/claim3/pass1-logs/*.log)
test "$((claim3_planar + claim3_1planar + claim3_unknown))" -eq 4881876
test "$claim3_unknown" -eq "$(wc -l < evidence/claim3/claim3b-input.g6)"

claim3b_records=$(wc -l < evidence/claim3/claim3b-input.g6)
claim3b_parents=$(sum_counter '#claim3b-parents' evidence/claim3/claim3b-logs/*.log)
test "$(sum_counter '#claim3b-records' evidence/claim3/claim3b-logs/*.log)" \
  -eq "$claim3b_records"
test "$(sum_counter '#claim3b-completed' evidence/claim3/claim3b-logs/*.log)" \
  -eq "$claim3b_records"
test "$(sum_counter '#1-flexible' evidence/claim3/claim3b-logs/*.log)" \
  -eq "$claim3b_parents"
test "$(sum_counter '#unknown' evidence/claim3/claim3b-logs/*.log)" -eq 0
```

## Verify Claim 4

### Prepare data

Generate and classify the Claim 4 graphs.  Each output line is `D` or a
`B`, `C`, or `R` record followed by its code.

```bash
generate_minibaum_parts 28 5 data/claim4-biconnected-nonplanar exact
mkdir evidence/claim4
mkdir data/claim4-classified evidence/claim4/preparation

printf '%s\0' data/claim4-biconnected-nonplanar/*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    classified="data/claim4-classified/$name.txt"
    ./prepare_cubic claim4-joint < "$input" \
      > "$classified" \
      2> "evidence/claim4/preparation/$name.log"
    test "$(wc -l < "$input")" -eq "$(wc -l < "$classified")"
  '

test "$(wc -l data/claim4-biconnected-nonplanar/*.g6 | awk 'END {print $1}')" \
  -eq 652157758
awk '
  { count[$1]++; if (NF == 2) print $2 }
  END { exit count["D"] != 254354420 || count["B"] != 350383009 ||
             count["C"] != 45970275 || count["R"] != 1450054 }
' data/claim4-classified/*.txt |
  sort --parallel="$JOBS" -u > data/claim4.g6
test "$(wc -l < data/claim4.g6)" -eq 3788793
```

### Verify records

```bash
mkdir evidence/claim4/certificate-logs \
      evidence/claim4/certificate-residues
make_shards data/claim4.g6 evidence/claim4/certificate-shards

printf '%s\0' evidence/claim4/certificate-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -verify-cubic -timeout=60 \
      -Cverify-cubic-residue="evidence/claim4/certificate-residues/$name.g6" \
      -colors=0 > "evidence/claim4/certificate-logs/$name.log" 2>&1
  '

cat evidence/claim4/certificate-residues/*.g6 | sort -u > evidence/claim4/residue.g6

claim4_planar=$(sum_counter '#planar' evidence/claim4/certificate-logs/*.log)
claim4_1planar=$(sum_counter '#1-planar' evidence/claim4/certificate-logs/*.log)
claim4_unknown=$(sum_counter '#unknown' evidence/claim4/certificate-logs/*.log)
test "$((claim4_planar + claim4_1planar + claim4_unknown))" -eq 3788793
test "$claim4_unknown" -eq "$(wc -l < evidence/claim4/residue.g6)"
```

### Verify residue parents

```bash
paste <(cat data/claim4-biconnected-nonplanar/*.g6) \
      <(cat data/claim4-classified/*.txt) |
awk '
  FILENAME == ARGV[1] { residue[$1] = 1; next }
  NF == 3 && ($3 in residue) { print $1; found[$3] = 1 }
  END { for (code in residue) if (!(code in found)) exit 1 }
' evidence/claim4/residue.g6 - |
sort --parallel="$JOBS" -u > data/claim4-residue-parents.g6

mkdir evidence/claim4/parent-logs
make_shards data/claim4-residue-parents.g6 evidence/claim4/parent-shards

printf '%s\0' evidence/claim4/parent-shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -colors=0 \
      > "evidence/claim4/parent-logs/$name.log" 2>&1
  '

claim4_parent_count=$(wc -l < data/claim4-residue-parents.g6)
test "$(sum_counter '#1-planar' evidence/claim4/parent-logs/*.log)" \
  -eq "$claim4_parent_count"
test "$(sum_counter '#non-1-planar' evidence/claim4/parent-logs/*.log)" -eq 0
test "$(sum_counter '#unknown' evidence/claim4/parent-logs/*.log)" -eq 0
```

## Verify Claim 5

Generate every biconnected, nonplanar cubic graph with $n=28$ and girth at
least 6, then verify it directly:

```bash
generate_minibaum_parts 28 6 data/claim5-biconnected-nonplanar
cat data/claim5-biconnected-nonplanar/*.g6 > data/claim5.g6
test "$(wc -l < data/claim5.g6)" -eq 4624501

mkdir evidence/claim5
mkdir evidence/claim5/logs
make_shards data/claim5.g6 evidence/claim5/shards
printf '%s\0' evidence/claim5/shards/shard-*.g6 | run_in_parallel '
    input=$1
    name=$(basename "$input" .g6)
    ./oops -i="$input" -colors=0 \
      > "evidence/claim5/logs/${name}.log" 2>&1
  '

test "$(sum_counter '#1-planar' evidence/claim5/logs/*.log)" -eq 4624501
test "$(sum_counter '#non-1-planar' evidence/claim5/logs/*.log)" -eq 0
test "$(sum_counter '#unknown' evidence/claim5/logs/*.log)" -eq 0
```
