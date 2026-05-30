#!/bin/sh
# Lint: every munit registration's display-name stem must equal its function
# pointer, so "X tests passed" can't lie about which routine was exercised.
# (See the swap_stride / idamax-as-isamax misregistration class.)
set -eu
f="tests/test_all.c"
bad=$(grep -oE '\{"/test_[a-z0-9_]+", *test_[a-z0-9_]+,' "$f" \
  | sed -E 's/\{"\/(test_[a-z0-9_]+)", *(test_[a-z0-9_]+),/\1 \2/' \
  | awk '$1 != $2 { print "  " $1 "  ->  " $2 }')
if [ -n "$bad" ]; then
  echo "FAIL: test registration name != function pointer in $f:"
  echo "$bad"
  exit 1
fi
echo "OK: all test registrations match their function pointers"
