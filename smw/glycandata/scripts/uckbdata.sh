#!/bin/sh
set -x

PYTHON=apython
PYGLY=../../../pygly
DATA=../data

$PYTHON $PYGLY/GlycanResource/main.py UniCarbKB allgtc | \
  $PYTHON ./normalizebasecomp.py -f UCKBCOMP -p 1 -r '^comp_' -a | \
  awk '$1 ~ /comp_/ {print "-",$1,$2; next} {print $1,"-",$2}' | \
  sort -u -k1n,1 -k2,2 -k3,3 | \
  awk '$2 ~ /comp_/ {print $2,$3; next} {print $1"\t"$3}' \
  > $DATA/uc2gtc.txt 2>/dev/null
$PYTHON $PYGLY/GlycanResource/main.py UniCarbKB allpub  | \
  $PYTHON ./normalizebasecomp.py -f UCKBCOMP -p 1 -r '^comp_' -a | \
  awk '$1 ~ /comp_/ {print "-",$1,$2; next} {print $1,"-",$2}' | \
  sort -k1n,1 -k2,2 -k3n,3 | \
  awk '$2 ~ /comp_/ {print $2,$3; next} {print $1"\t"$3}' \
  > $DATA/uc2pubmed.txt 2>/dev/null
$PYTHON $PYGLY/GlycanResource/main.py UniCarbKB alltaxa | \
  $PYTHON ./normalizebasecomp.py -f UCKBCOMP -p 1 -r '^comp_' -a | \
  awk '$1 ~ /comp_/ {print "-",$1,$2; next} {print $1,"-",$2}' | \
  sort -k1n,1 -k2,2 -k3n,3 | \
  awk '$2 ~ /comp_/ {print $2,$3; next} {print $1"\t"$3}' \
  > $DATA/uc2taxa.txt   2>/dev/null
