#!/bin/sh

DIR=`dirname $0`                                                                                                             
# MAIN="GlycoCT2Image"
# JAR="$DIR/../pygly/GlycoCT2ImageBundle.jar"
MAIN="GlycanImageCmdline"
JAR="$DIR/../pygly/GlycanBuilder2.jar"
JAR=`readlink -f $JAR`
JAVA=${JAVA:-java}

NOTATION=snfg
DISP=normalinfo
SCALE=1.0
ORIENT=RL
REDEND=true
OPAQUE=true
FORMAT=png
FORCE=false
OUTDIR="."
VERBOSE=0

help() {
  if [ "$1" != "" ]; then
    echo "" 1>&2
    echo "$1" 1>&2
  fi
  cat <<EOF | fmt -w80 -s 1>&2

bulkimg.sh [ options ] <glycan1.seq> <glycan2.seq> ...

Options:
  -n    Notation. One of snfg, cfg Default: $NOTATION.
  -d    Display. One of normalinfo, compact. Default: $DISP.
  -s    Scale. Default: $SCALE.
  -O    Orientation. Default: $ORIENT.
  -r    Reducing-end. Default: $REDEND.
  -w    White-background (opaque). Default: $OPAQUE.
  -f    Format. One of png, svg. Default: $FORMAT.
  -F    Force (overwrite). Default: $FORCE.
  -o    Ouput directory. Default: $OUTDIR
  -v    Verbose.
  -h    Help text.
EOF
  exit 1;
}
                                                                                                                             
while getopts "n:d:s:O:r:w:f:F:o:vh" o ; do 
        case $o in                                                                                                           
                n ) NOTATION="$OPTARG";;
                d ) DISP="$OPTARG";;
                s ) SCALE="$OPTARG";;
                O ) ORIENT="$OPTARG";;
                r ) REDEND="$OPTARG";;
                w ) OPAQUE="$OPTARG";;
                f ) FORMAT="$OPTARG";;
                F ) FORCE="$OPTARG";;
                o ) OUTDIR="$OPTARG";;
                v ) VERBOSE=1;;
                h ) help "";;
                * ) help "Invalid option: -$OPTARG"
        esac                                                                                                                 
done                                                                                                                         
                                                                                                                             
shift $(($OPTIND - 1))                                                                                                       

if [ $VERBOSE -eq 1 ]; then
  echo $JAVA -cp "$JAR" $MAIN notation "$NOTATION" display "$DISP" scale "$SCALE" orient "$ORIENT" redend "$REDEND" opaque "$OPAQUE" format "$FORMAT" force "$FORCE" outdir "$OUTDIR" $@ 1>&2
fi
exec $JAVA -cp "$JAR" $MAIN notation "$NOTATION" display "$DISP" scale "$SCALE" orient "$ORIENT" redend "$REDEND" opaque "$OPAQUE" format "$FORMAT" force "$FORCE" outdir "$OUTDIR" $@ 2>/dev/null | egrep -w -v '(org.glycoinfo|DEBUG|GlycanImageCmdline.main|org.eurocarbdb.application.glycanbuilder)'
