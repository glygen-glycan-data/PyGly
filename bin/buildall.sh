#!/bin/sh
CMD="build"
if [ "$1" = "-f" ]; then
    CMD="force-build"
    shift
fi
THEDIR=`dirname $0`
THEDIR=`readlink -f $THEDIR`
if [ -z "$PYTHONPATH" ]; then
  PYTHONPATH="$THEDIR/.."
else
  PYTHONPATH="$PYTHONPATH:$THEDIR/.."
fi
export PYTHONPATH
for gdb in "$THEDIR"/../data/*.{gct,gdb}; do
    if [ ! -f $gdb ]; then
	continue
    fi
    gdb=`readlink -f $gdb`
    echo $gdb
    $THEDIR/GlyDbIndex.py $CMD $gdb
done
