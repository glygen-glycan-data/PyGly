#!/bin/sh
# set -x

PYGLY=../../../scripts
DATA=../data

$PYGLY/glyres.py UniCarbKBSourceFile allgtc | \
  ./normalizebasecomp.py -f UCKBCOMP -p 1 -r '^comp_' -a | \
  awk '$1 ~ /comp_/ {print "-",$1,$2; next} {print $1,"-",$2}' | \
  sort -u -k1n,1 -k2,2 -k3,3 | \
  awk '$2 ~ /comp_/ {print $2,$3; next} {print $1"\t"$3}' \
  | uniq > $DATA/uc2gtc.txt 2>/dev/null
$PYGLY/glyres.py UniCarbKBSourceFile allpubs  | \
  ./normalizebasecomp.py -f UCKBCOMP -p 1 -r '^comp_' -a | \
  awk '$1 ~ /comp_/ {print "-",$1,$4; next} {print $1,"-",$4}' | \
  sort -k1n,1 -k2,2 -k3n,3 | \
  awk '$2 ~ /comp_/ {print $2,$3; next} {print $1"\t"$3}' \
  | uniq > $DATA/uc2pubmed.txt 2>/dev/null
$PYGLY/glyres.py UniCarbKBSourceFile alltaxa | \
  ./normalizebasecomp.py -f UCKBCOMP -p 1 -r '^comp_' -a | \
  awk '$1 ~ /comp_/ {print "-",$1,$3; next} {print $1,"-",$3}' | \
  sort -k1n,1 -k2,2 -k3n,3 | \
  awk '$2 ~ /comp_/ {print $2,$3; next} {print $1"\t"$3}' \
  | uniq > $DATA/uc2taxa.txt   2>/dev/null
