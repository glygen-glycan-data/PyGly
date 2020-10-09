#!/bin/sh

#
#
#
set -x

LIB="SA000001.5decoy.tsv"
WIN="ME000014.windows.tsv"

BIN="/home/nedwards/projects/OpenSWATH/docker/bin"
for f in "$*"; do
  BASE=`basename "$f" .mzML.gz`
  DIR=`dirname "$f"`
  TMP=`mktemp -d -p /scratch -t openswath.XXXXXX`
  $BIN/OpenSwathWorkflow.sh \
    -in "$f" \
    -tr "$DIR/$LIB" \
    -Scoring:Scores:use_rt_score false \
    -sort_swath_maps \
    -readOptions cache \
    -threads 4 \
    -tempDirectory "$TMP/" \
    -batchSize 1000 \
    -extra_rt_extraction_window 100000 \
    -rt_extraction_window -1 \
    -swath_windows_file "$DIR/$WIN" \
    -out_tsv "${BASE}_table.tsv" \
    -out_chrom "${BASE}_chrom.mzML" \
    > ${BASE}.log 2>&1
  rm -rf $TMP*
  
  apython fdr.py "${BASE}_table.tsv" 5 > "${BASE}_fdr.txt"
  apython getcoeff_fdr.py "${BASE}_table.tsv" "$DIR/$LIB" "${BASE}_fdr.txt" > "${BASE}.trafoXML"

  TRAFOXML="${BASE}.trafoXML"
  
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
    -out_tsv "${BASE}_final_table.tsv" \
    -out_chrom "${BASE}_final_chrom.mzML" \
    > ${BASE}.log 2>&1
  rm -rf $TMP*
done