#!/bin/sh
THEDIR=`dirname $0`
FROMDB=$1
FROMNAME=$2
TODB=$3
TONAME=$4
if [ -z "$TONAME" ]; then
  echo "copyglycan.sh fromdb fromname todb toname"
  exit 1
fi
( cd $TODB; unzip -p ../$FROMDB.zip $FROMNAME.txt > $TONAME.txt )
$THEDIR/makedb.sh $TODB
