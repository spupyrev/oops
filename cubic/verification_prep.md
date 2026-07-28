Run one preparation stage first. Do not start full Claim 4 verification yet.

This uses 256 Minibaum residue tasks and automatically uses the machine’s core count—32, 48, or otherwise. There is no `JOBS=48` and no resume machinery.

```bash
set -euo pipefail

NAUTY=/path/to/nauty2_9_3
MINIBAUM=./minibaum5
PARTS=256
export NAUTY MINIBAUM PARTS LC_ALL=C

mkdir -p data/claim4-pairs evidence/claim4

seq 0 $((PARTS - 1)) |
  xargs -n1 -P "$(nproc)" bash -c '
    set -euo pipefail
    part=$1
    output=$(printf "data/claim4-pairs/part-%03d.pairs" "$part")

    "$MINIBAUM" 28 5 s g m "$part" "$PARTS" |
      "$NAUTY/pickg" -q -c2 -g5 |
      "$NAUTY/planarg" -q -v |
      ./prepare_cubic claim4 > "$output"
  ' bash
```

If you want fully serial execution, replace `-P "$(nproc)"` with `-P 1`. On the 48-core server, serial preparation would be impractically slow; the adaptive form should reproduce the measured 25–30× effective speedup.

Then perform the single global merge:

```bash
sort -u data/claim4-pairs/part-*.pairs \
  > data/claim4-marked-pairs.txt

./prepare_cubic group-claim4 \
  < data/claim4-marked-pairs.txt \
  > data/claim4-marked-groups.g6

pair_count=$(wc -l < data/claim4-marked-pairs.txt)
group_count=$(wc -l < data/claim4-marked-groups.g6)

printf 'distinct marked pairs: %s\n' "$pair_count"
printf 'distinct groups:       %s\n' "$group_count"

test "$group_count" -eq \
  "$(cut -d' ' -f1 data/claim4-marked-pairs.txt | uniq | wc -l)"
```

Record the exact marker distribution:

```bash
cut -d' ' -f1 data/claim4-marked-pairs.txt |
  uniq -c |
  awk '{ groups[$1]++ }
       END {
         for (markers in groups)
           print markers, groups[markers]
       }' |
  sort -n \
  > evidence/claim4/marker-histogram.txt

cat evidence/claim4/marker-histogram.txt

awk -v pairs="$pair_count" -v groups="$group_count" \
  'BEGIN { printf "mean markers/group: %.6f\n", pairs/groups }'
```

Finally, benchmark 100 uniformly spaced groups from the globally merged file:

```bash
awk -v count="$group_count" -v wanted=100 '
  BEGIN { step = count / wanted }
  NR >= int((selected + 0.5) * step) + 1 && selected < wanted {
    print
    selected++
  }
' data/claim4-marked-groups.g6 \
  > data/claim4-sample100.g6

/usr/bin/time -v -o evidence/claim4/sample100.time \
  ./oops -i=data/claim4-sample100.g6 -verify-cubic -colors=0 \
  > evidence/claim4/sample100.log 2>&1

tail -5 evidence/claim4/sample100.log
cat evidence/claim4/sample100.time
```

The required verifier counters are exactly:

```text
#planar = 0; #1-planar = 100
#3-flexible = 0; #2-flexible = 0
#1-flexible = 0; #5-cycle-cores = 0; #marked-star-groups = 100
```

Expected numbers

These are forecasts, not verification targets:

- Raw filtered parents/pair records before `sort -u`: probably about 640–648M; certainly no more than the 652,159,389 connected minimum-girth-five inputs.
- Distinct marked pairs after `sort -u`: roughly 400–500M.
- Distinct groups: current diagnostic estimate is about 30–31M, but this is the number the run must determine exactly.
- Mean markers per globally merged group: plausibly 10–17, with a hard maximum of 26.
- Preparation plus global sorting: approximately 12–15 hours on the 48-core server.
- The 100-group sample:
  - 27 seconds if the old 269 ms/group rate survives global merging.
  - More realistically, perhaps 3–8 minutes if the larger marker sets cost several seconds per group.
  - If it takes over 15 minutes, Claim 4 is already far outside the three-day budget.

Decision after the run

Let \(N\) be `group_count`, in millions, and \(t\) the sample seconds per group. The projected verification time is:

```text
N × 1,000,000 × t / effective_speedup
```

For hours:

```bash
awk -v groups="$group_count" -v seconds_per_group="REPLACE_ME" '
  BEGIN {
    printf "25x: %.1f hours\n", groups*seconds_per_group/25/3600
    printf "30x: %.1f hours\n", groups*seconds_per_group/30/3600
  }'
```

Even under the optimistic old rate of 0.269 s/group:

- At 19.5M groups, the joint Claims 3–4 pipeline is near the 72-hour limit at 25×.
- Above 23.4M groups, it exceeds 72 hours even at 30×.
- At the current rough expectation of 30M groups, verification alone is 75–90 hours, before preparation and Claim 3.

Therefore, if `group_count` is around 30M, or the final sample is materially slower than 269 ms/group, do not launch full Claim 4 verification.