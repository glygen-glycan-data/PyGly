#!/bin/sh
THEDIR=`dirname $0`
THEDIR=`readlink -f $THEDIR`
if [ -z "$PYTHONPATH" ]; then
  PYTHONPATH="$THEDIR/.."
else
  PYTHONPATH="$PYTHONPATH:$THEDIR/.."
fi
export PYTHONPATH
exec $THEDIR/GlyDbIndex.py "$@"
