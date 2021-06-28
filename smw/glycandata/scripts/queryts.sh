#!/bin/sh
if [ ! -f "$2" ]; then
  echo "Usage: $0 <database>.tdb <queryfile>"
  exit 1
fi
if [ ! -d "$1" ]; then
  echo "Usage: $0 <database>.tdb <queryfile>"
  exit 1
fi
tdb2.tdbquery --loc "$1" --file "$2" --results "TSV"
