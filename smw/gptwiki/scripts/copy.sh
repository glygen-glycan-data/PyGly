#!/bin/sh
if [ $# -lt 2 ]; then
  echo "copy.sh FROM TO" 1>&2
  exit 1;
fi
set -x
./loadsite.py --smwenv "$2" ../wiki
./copysite.py "$1" "$2"
./refresh.py --smwenv "$2"
./refresh.py --smwenv "$2" - 
