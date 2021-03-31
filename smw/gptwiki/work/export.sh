#!/bin/sh

if [ "$1" = "" ]; then
  echo "Usage: export.sh [PROD|DEV] *.mzML.gz"
  exit 1
fi
OUTPUT="$HOME/projects/GlyGen/PyGly/smw/gptwiki/repository/$1"
shift;
if [ ! -d "$OUTPUT" ]; then
  echo "Usage: export.sh [PROD|DEV] *.mzML.gz"
  exit 1
fi
for SPEC in "$@"; do
  case "$SPEC" in                                                                                                            
    *.mzML.gz) BASE=`basename $SPEC .mzML.gz`;;                                                                              
    *.msp) BASE=`basename $SPEC .msp`;;                                                                                      
  esac
  if [ -d $OUTPUT/$BASE ]; then
    rm -rf $OUTPUT/$BASE
  fi
  cp -r $BASE $OUTPUT
done
