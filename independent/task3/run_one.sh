#!/bin/bash
cd "$(dirname "$0")"
mkdir -p logs
timeout 18000 magma -b Kind:=$1 N:=$2 rankcheck.m < /dev/null > logs/${1}_$2.log 2>&1
echo "$1 $2 exit $?" >> status.txt
