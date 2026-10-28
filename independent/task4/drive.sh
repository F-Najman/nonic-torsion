#!/bin/bash
# usage: drive.sh p e K k
cd "$(dirname "$0")"
p=$1; e=$2; K=$3; k=$4
mkdir -p hist
out=hist/h_${p}_${e}_${k}of${K}.txt; det=hist/d_${p}_${e}_${k}of${K}.txt
rm -f $out $det
echo "p=$p; e=$e; K=$K; k=$k; outfile=\"$out\"; detfile=\"$det\"; watch=[27,39,87,97,101,29,31,33,37,41,43,47,49,51,53,57,58,59,61,62,65,67,69,71,73,75,79,83,89,45,55,38]; read(\"hist.gp\");" | gp -q --default parisizemax=2000000000 > hist/log_${p}_${e}_${k}of${K}.txt 2>&1
echo "exit $?" >> hist/log_${p}_${e}_${k}of${K}.txt
