#!/bin/bash

assert() {
  local condition="$1"
  local expect_failure="$2"
  local message="${3:-$condition}"

  if eval "$condition"; then
    if [[ "$expect_failure" -eq 1 ]]; then
      echo "❌ FAIL (unexpected success): $message"
      exit 1
    else
      echo "✅ PASS: $message"
    fi
  else
    if [[ "$expect_failure" -eq 1 ]]; then
      echo "✅ PASS (expected failure): $message"
    else
      echo "❌ FAIL: $message"
      exit 1
    fi
  fi
}

assert_pass() { assert "$1" 0 "$2"; }
assert_fail() { assert "$1" 1 "$2"; }


OP=${PWD}/oops
assert_pass "[ -f ${OP} ]" "Binary ${OP} exists"
DATA_DIR=${PWD}/data
assert_pass "[ -d ${DATA_DIR} ]" "Folder ${DATA_DIR} exists"

# no params
assert_fail "${OP}"
assert_pass "${OP} --help"

# read from different formats
assert_fail "${OP} -verbose=1 -i=${DATA_DIR}/test3.xxx"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test2.g6"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test3.s6"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test4.dot"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test5.gml"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test6.graphml"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test7.graphml"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=1"
# TODO: read directed

# gen random
assert_pass "${OP} -verbose=1 -i=gen-complete -graphs=5 -n=10"
assert_pass "${OP} -verbose=1 -i=gen-complete-bipartite -graphs=1 -n=10"
assert_pass "${OP} -verbose=1 -i=gen-grid -graphs=2 -n=20"

# write to files
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=0 -o=/tmp/test1.gml"

# test misc options
assert_fail "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -invalid-option"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -cross2"
assert_fail "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -cross1"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -cross2 -cross1"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -cross2 -cross1 -ic -timeout=2"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -cross2 -cross1 -nic -timeout=1"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -satsuma"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -breakID"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -cross1 -cross2 -satsuma"
assert_fail "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=1 -move-planar"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=1 -move-planar -cross2"

sparsify

# timeout
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test8.cfg -timeout=5"

# dimacs
assert_pass "rm -f /tmp/test8.dimacs"
assert_fail "[ -f /tmp/test8.dimacs ]"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test8.cfg -dimacs=/tmp/test8.dimacs"
assert_pass "[ -f /tmp/test8.dimacs ]"

# dimacs-result e2e
#assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test8.cfg -dimacs-result=/tmp/test8.dimacs"

echo "✅ ALL TESTS PASSED"