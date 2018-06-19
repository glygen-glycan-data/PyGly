#!/bin/sh
NAME="$1"
if [ -z "$NAME" ]; then
  echo "addimages.sh database.[gct|gdb]"
  exit 1
fi
case $NAME in
    /*) ;;
     *) NAME="$PWD/$NAME";;
esac
TMPDIR=`mktemp -d` || exit 1
THEDIR=`dirname $0`
JAR="$THEDIR/../pygly/GlycoCT2ImageBundle.jar"
cd "$TMPDIR"; 
unzip -q "$NAME" '*.txt'; 
ls *.txt | xargs -n 40 chmod +r 
ls *.txt | xargs -n 10 java -cp "$JAR" GlycoCT2Image scale 1.0 display compact
zip -q "$NAME" *.png
rm -rf $TMPDIR
