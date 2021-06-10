#!/bin/sh

DIR=`dirname $0`                                                                                                             
DIR=`readlink -f $DIR`

NOTATION=snfg
DISP=normalinfo
SCALE=1.0
ORIENT=RL
REDEND=true
OPAQUE=true
FORMAT=png
FORCE=false
OUTDIR="."
PATTERN="*.txt"
PROCS="2"
BATCH="20"
VERBOSE=""

help() {
  if [ "$1" != "" ]; then
    echo "" 1>&2
    echo "$1" 1>&2
  fi
  cat <<EOF | fmt -w80 -s 1>&2

allimg.sh [ options ] <glycanseq-directory>

Options:
  -n    Notation. One of snfg, cfg Default: $NOTATION.
  -d    Display. One of normalinfo, compact. Default: $DISP.
  -s    Scale. Default: $SCALE.
  -O    Orientation. Default: $ORIENT.
  -r    Reducing-end. Default: $REDEND.
  -w    White-background (opaque). Default: $OPAQUE.
  -f    Format. One of png, svg. Default: $FORMAT.
  -F    Force (overwrite). Default: $FORCE.
  -o    Ouput directory. Default: $OUTDIR.
  -P    Processors. Default: $PROCS.
  -p    Filename pattern. Default: $PATTERN.
  -N    Batch size. Default: $BATCH.
  -v    Verbose.
  -h    Help text.
EOF
  exit 1;
}
                                                                                                                             
while getopts "p:P:N:n:d:s:O:r:w:f:F:o:vh" o ; do 
        case $o in                                                                                                           
                p ) PATTERN="$OPTARG";;
                P ) PROCS="$OPTARG";;
                N ) BATCH="$OPTARG";;
                n ) NOTATION="$OPTARG";;
                d ) DISP="$OPTARG";;
                s ) SCALE="$OPTARG";;
                O ) ORIENT="$OPTARG";;
                r ) REDEND="$OPTARG";;
                w ) OPAQUE="$OPTARG";;
                f ) FORMAT="$OPTARG";;
                F ) FORCE="$OPTARG";;
                o ) OUTDIR="$OPTARG";;
                v ) VERBOSE=" -v";;
                h ) help "";;
                * ) help "Invalid option: -$OPTARG"
        esac                                                                                                                 
done                                                                                                                         

shift $(($OPTIND - 1))

mkdir -p "$OUTDIR"
if [ "$1" = "" ]; then
  FILES='.'
else
  FILES="$1"
fi
find "$FILES" -name "$PATTERN" -type f | sort | xargs -n "$BATCH" -P "$PROCS" $DIR/bulkimg.sh -n "$NOTATION" -d "$DISP" -s "$SCALE" -O "$ORIENT" -r "$REDEND" -w "$OPAQUE" -f "$FORMAT" -F "$FORCE" -o "$OUTDIR" $VERBOSE
