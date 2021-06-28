#!/bin/sh

set -x

DIR=`pwd`
DIR=`readlink -f "$DIR"`
DATA=$DIR/../data
DATA=`readlink -f "$DATA"`

DDATRANS="--mintrans 4 --maxtrans 100"
DIATRANS="--mintrans 4 --maxtrans 8"

for sa in SA000001 SA000002 SA000004 SA000005; do
  $DIR/extract_transitions.py --sample $sa --inst Orbitrap --acqtype DDA $DDATRANS > $DATA/$sa.Orbitrap.DDA.5decoys.tsv
done
for sa in SA000001; do
  $DIR/extract_transitions.py --sample $sa --inst Orbitrap --acqtype DIA $DIATRANS > $DATA/$sa.Orbitrap.DIA.5decoys.tsv
  $DIR/extract_transitions.py --sample $sa --inst TripleTOF --acqtype DIA $DIATRANS > $DATA/$sa.TripleTOF.DIA.5decoys.tsv
done
