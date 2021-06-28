#!/bin/sh
set -x
if [ ! -f "$1" ]; then
  echo "Usage: $0 <rdffile>.rdf.gz"
  exit 1
fi
DB=`basename "$1" .rdf.gz`
DB="$DB.tdb"
if [ -d "$DB" ]; then
  rm -rf "$DB"
fi
tdb2.tdbloader -loc "$DB" "$1"
