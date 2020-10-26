#!/bin/sh

#
#
#

set -x

LIB="SA000001.5decoy.tsv"
TRAFOXML="/home/lzhang/OpenSWATH/Erlich/DataMS_UO1_SWATH_11202019_Erlich7.centroid.trafoXML"
WIN="ME000014.windows.tsv"

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
    -rt_norm "$TRAFOXML" \
    -extra_rt_extraction_window 100000 \
    -rt_extraction_window 840 \
    -swath_windows_file "$DIR/$WIN" \
    -out_tsv "${BASE}_erlich7_table.tsv" \
    -out_chrom "${BASE}_erlich7_chrom.mzML" \
    > ${BASE}.log 2>&1
  rm -rf $TMP*
done

