#!/bin/sh
set -x
GT="$2"
SMW="$1"
( cd $GT; ./2_clear_data.sh )
( cd $GT; rm -rf data/glygensvg ; mkdir -p data/glygensvg )
./dumpglycoct.py "$SMW" ../glycoTree/data/gct
cat images-snfg-extended-svg.tbz.0* | tar -xjf - -C $GT/data/glygensvg --strip-components 3
( cd $GT; $GT/populate_N-tree.sh ./data/glygensvg/ )
