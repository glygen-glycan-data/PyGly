#!/bin/sh
if [ "$1" = "" ]; then
  echo "Usage: gittagrelease.sh <release-tag>"
  exit 1;
fi
CURTAG="GlyGen-GlycanData-Export-Current"
set -x
git pull
git tag -a -m "Pointer to release $1" -f "$CURTAG" "$1"
git push -f --tags 
