#!/bin/sh
set -x
GT="$2"
SMW="$1"
# ( cd $GT; ./2_clear_data.sh )
# ( cd $GT; rm -rf data/glygensvg ; mkdir -p data/glygensvg )
mkdir -p $GT/data/glygengct 
./dumpglycoct.py "$SMW" $GT/data/glygengct
# ( cd $GT/data/glygengct ; find . -type f ) | tr './' '  ' | awk '{print "snfg/extended/svg/"$1".svg"}' > $GT/data/glygensvg.lst
# cat images-snfg-extended-svg.tbz.0* | tar -xjf - -C $GT/data/glygensvg --strip-components 3 -T $GT/data/glygensvg.lst
# ( cd $GT; $GT/populate_N-tree.sh ./data/glygensvg/ )
