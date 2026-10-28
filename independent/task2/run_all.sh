#!/bin/bash
cd "$(dirname "$0")/../.."   # the repository root
seq 1 20 | parallel -j 4 'timeout 18000 magma -b Idx:={} independent/task2/check_entry.m > independent/task2/entry_$(printf %02d {}).log 2>&1; echo "idx {} exit $?" >> independent/task2/run_all.status'
echo ALL_DONE >> independent/task2/run_all.status
