#!/bin/sh
THEDIR=`dirname $0`
THEDIR=`readlink -f $THEDIR`
if [ -z "$PYTHONPATH" ]; then
  PYTHONPATH="$THEDIR/.."
else
  PYTHONPATH="$PYTHONPATH:$THEDIR/.."
fi
export PYTHONPATH
for gdb in "$THEDIR"/../data/*.{gdb,gct}; do
    if [ ! -f $gdb ]; then
        continue
    fi
    gdb=`readlink -f $gdb`
    $THEDIR/GlyDbIndex.py count $gdb "$@"
done
