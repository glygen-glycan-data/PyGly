#!/bin/sh
set -x

DIA_DIRS=" \
  GP_SWATH_Orbitrap \
  Erlich_plasma_50acn/endogenous \
  Erlich_plasma_50acn/explicit \
  Erlich_plasma_water/endogenous \
  Erlich_plasma_water/explicit \
"

DIR1=`pwd`
DIR1=`readlink -f "$DIR1"`
DIR=`pwd`/../scripts
DIR=`readlink -f "$DIR"`

for workdir in $DIA_DIRS; do 
  ( cd $workdir; \
    $DIR/OpenSWATH.sh *.mzML.gz > process.log 2>&1 ; \
    $DIR/OpenSWATH2JSON.sh *.mzML.gz >> process.log 2>&1 ; \
    $DIR1/export.sh DEV *.mzML.gz > export.log 2>&1 \
  )
done

