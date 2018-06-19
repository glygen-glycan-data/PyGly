#!/bin/sh
THEDIR=`dirname $0`
THEDIR=`readlink -f $THEDIR`
NAME=$1
if [ -z "$NAME" ]; then
  echo "makedb.sh directory"
  exit 1
fi
rm -f $NAME.gct*
rm -f $NAME/*.png
( cd $NAME; java -cp ../../pygly/GlycoCT2ImageBundle.jar GlycoCT2Image *.txt )
( cd $NAME; zip ../$NAME.gct *.txt *.png )
if [ -z "$PYTHONPATH" ]; then
  PYTHONPATH="$THEDIR/.."
else
  PYTHONPATH="$PYTHONPATH:$THEDIR/.."
fi
export PYTHONPATH
$THEDIR/GlyDbIndex.py force-build $NAME.gct
