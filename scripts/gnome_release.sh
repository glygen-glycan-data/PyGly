#!/bin/sh

set -x

VER="$1"
TITLE="$2"
shift; shift

if [ "$VER" = "" -o "$TITLE" = "" ]; then
  echo "Usage: gnome_release.sh VX.Y.Z <title> [ <comment> <comment> ... ]" 1>&2 
  exit 1
fi

TMP=`mktemp` || exit 1
touch $TMP
if [ "$1" ]; then 
  for comment in "$@"; do 
    echo "* $comment" >> $TMP
  done
fi

gh release create "$VER" -R glygen-glycan-data/GNOme -t "$TITLE" -F $TMP GNOme*.owl GNOme*.obo GNOme*.json

rm -f $TMP
