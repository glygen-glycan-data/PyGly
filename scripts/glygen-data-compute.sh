#!/bin/sh

DOLOG=1
if [ "$1" != "DAEMON" -a "$1" != "--nolog" ]; then
   setsid "$0" DAEMON "$@" &
   exit 0
else
   if [ "$1" == "--nolog" ]; then
       DOLOG=0
   fi
   shift;
fi

BASE=`basename "$0" .sh`
DIR=`dirname "$0"`
if [ $DOLOG -eq 1 ]; then
  exec </dev/null >$DIR/$BASE.`date +"%Y-%m-%d_%H-%M-%S"`.log 2>&1 
fi

DATESTAMP=`date "+%Y%m%d"`
FORCE=0
NPROC=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d) DATESTAMP="$2"
            shift 2
            ;;
        -n) NPROC="$2"
            shift 2
            ;;
        -f) FORCE=1
            shift 1
            ;;   
        *) break ;;
    esac
done

if [ "$NPROC" != "" ]; then
  NPROC="-j $NPROC"
fi

if [ $FORCE -eq 1 ]; then
    echo make -f $DIR/$BASE.mk DATESTAMP="$DATESTAMP" clean
    make -f $DIR/$BASE.mk -j 1 DATESTAMP="$DATESTAMP" clean
fi

echo make -f $DIR/$BASE.mk DATESTAMP="$DATESTAMP" $@
make -f $DIR/$BASE.mk $NPROC DATESTAMP="$DATESTAMP" $@
