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
assert_pass "rm -f /tmp/test.gml"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=5 -o=/tmp/test.gml"
assert_pass "[ -f /tmp/test.gml ]"

assert_pass "rm -f /tmp/test.gml"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=5 -o=/tmp/test.gml"
assert_pass "[ -f /tmp/test.gml ]"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=5 -o=txt"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=5 -o=/tmp/test.txt"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=5 -o=/tmp/test.dot"

# test misc options
assert_fail "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -invalid-option"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -sat"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -sat=0"
# assert_fail "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -unsat"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -sat -unsat"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -sat -unsat -ic -timeout=2"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -sat -unsat -nic -timeout=1"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -satsuma"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -breakID"
assert_pass "${OP} -verbose=0 -i=${DATA_DIR}/test1.cfg -unsat -sat -satsuma"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test1.cfg -part=1 -move-planar"

# test sparsify

# timeout
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test8.cfg -timeout=5"

# dimacs
assert_pass "rm -f /tmp/test8.dimacs"
assert_fail "[ -f /tmp/test8.dimacs ]"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test8.cfg -dimacs=/tmp/test8.dimacs"
assert_pass "[ -f /tmp/test8.dimacs ]"

# dimacs-result e2e
assert_pass "rm -f /tmp/test9.dimacs"
assert_fail "[ -f /tmp/test9.dimacs ]"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test9.dot -dimacs=/tmp/test9.dimacs"
assert_pass "diff /tmp/test9.dimacs data/test9.dimacs"
assert_pass "${OP} -verbose=1 -i=${DATA_DIR}/test9.dot -dimacs-result=${DATA_DIR}/test9.dimacs_result"

echo "✅ ALL TESTS PASSED"