#!/bin/sh
set -x
GT="$2"
SMW="$1"
mkdir -p $GT/data/glygengct_nlinked
./dumpglycoct.py "$SMW" $GT/data/glygengct_nlinked N-linked
mkdir -p $GT/data/glygengct_olinked
./dumpglycoct.py "$SMW" $GT/data/glygengct_olinked O-linked
