#!/bin/sh

OUTDIR="/data/projects/CTGRC/GPTWiki/repository/DEV"
for SPEC in "$@"; do
  case "$SPEC" in                                                                                                            
    *.centroid.mzML.gz) BASE=`basename $SPEC .centroid.mzML.gz`;; 
    *.mzML.gz) BASE=`basename $SPEC .mzML.gz`;;                                                                              
    *.msp) BASE=`basename $SPEC .msp`;;                                                                                      
  esac
  if [ "$BASE" != "" ]; then
    echo "$BASE..."
    rsync --progress -a --delete $BASE/ trypsin:$OUTDIR/$BASE
  fi
done
