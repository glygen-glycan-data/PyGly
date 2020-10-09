#!/bin/sh

#
#
#
set -x

LIB="SA000001.5decoy.tsv"
IRT="nrtgp.TraML"
WIN="ME000005.windows.tsv"

BIN="/home/nedwards/projects/OpenSWATH/docker/bin"
for f in "$*"; do
  BASE=`basename "$f" .mzML.gz`
  DIR=`dirname "$f"`
  TMP=`mktemp -d -p /scratch -t openswath.XXXXXX`
  $BIN/OpenSwathWorkflow.sh \
    -in "$f" \
    -tr "$DIR/$LIB" \
    -sort_swath_maps \
    -readOptions cache \
    -threads 4 \
    -tempDirectory "$TMP/" \
    -batchSize 1000 \
    -tr_irt "$DIR/$IRT" \
    -extra_rt_extraction_window 100000 \
    -rt_extraction_window 840 \
    -swath_windows_file "$DIR/$WIN" \
    -out_tsv "${BASE}_table.tsv" \
    -out_chrom "${BASE}_chrom.mzML" \
    > ${BASE}.log 2>&1
  rm -rf $TMP*
done

