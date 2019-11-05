#!/bin/sh
set -x

PYTHON=apython
PYGLY=../../../pygly
DATA=../data

$PYTHON $PYGLY/GlycanResource.py UniCarbKB allgtc  | sort -k1n,1 -k2,2 > $DATA/uc2gtc.txt    2>/dev/null
$PYTHON $PYGLY/GlycanResource.py UniCarbKB allpub  | sort -k1n,1 -k2,2 > $DATA/uc2pubmed.txt 2>/dev/null
$PYTHON $PYGLY/GlycanResource.py UniCarbKB alltaxa | sort -k1n,1 -k2,2 > $DATA/uc2taxa.txt   2>/dev/null

