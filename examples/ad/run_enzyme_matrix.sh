#!/bin/zsh
# Run every Enzyme case × mode in its own process under a timeout; one summary line each.
# usage: examples/ad/run_enzyme_matrix.sh [julia-channel, e.g. +1.12] [project dir]  > log
cd "$(dirname "$0")/.."
JL=${1:-}
PROJ=${2:-.}
for c in truss_area truss_shape frame_lineload frame_joints springs variable_element diaphragm shells shells_tri cases spaceframe; do
  for m in reverse forward; do
    out=$(timeout 1500 julia $JL --project=$PROJ ad/enzyme_matrix.jl $c $m 2>&1)
    code=$?
    line=$(echo "$out" | grep '^RESULT' | tail -1)
    if [ -n "$line" ]; then echo "$line"
    elif [ $code -eq 124 ]; then echo "RESULT $c $m TIMEOUT(1500s)"
    else echo "RESULT $c $m FAIL: $(echo "$out" | grep -m1 -E 'ERROR|Error|Assertion|error' | cut -c1-200)"
    fi
  done
done
