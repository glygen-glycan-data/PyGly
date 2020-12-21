#!/bin/sh

set -x

help() {
  if [ "$1" != "" ]; then
    echo "" 1>&2
    echo "$1" 1>&2
  fi
  cat <<EOF | fmt -w80 -s 1>&2

OpenSWATH2JSON.sh [ options ] <spectra>.mzML.gz ...

Options:
  -p    Parameters file. Default: ./params.txt.
  -f    FDR threshold (in percent). Default: 1%.
  -s    Score threshold. Default: 1.5.
  -o    Output directory. Default: .
  -h	Help text.

Parameter file sets the follwing variables:

  TRANSITIONS=<Transition-Library-File-With-Decoys>
  NDECOYS=<Number-of-Decoy-Replicates>

EOF
  exit 1;
}

DIR=`dirname $0`
DIR=`readlink -f $DIR`

VERBOSE=0
OUTDIR="."
FDR="1"
SCORE="1.5"
PARAMS="./params.txt"
while getopts "p:o:f:s:vh" o ; do
        case $o in
	        p ) PARAMS="$OPTARG";;
	        o ) OUTDIR="$OPTARG";;
	        f ) FDR="$OPTARG";;
	        s ) SCORE="$OPTARG";;
	        v ) VERBOSE=1;;
                h ) help "";;
                * ) help "Invalid option: -$OPTARG"
        esac
done

shift $(($OPTIND - 1)) 

if [ ! -f "$PARAMS" ]; then
  help "Parameter file \"$PARAMS\" could not be opened"
fi

PARAMS=`readlink -f "$PARAMS"`
. $PARAMS

if [ "$TRANSITIONS" = "" ]; then
  help "TRANSITIONS missing from parameter file $1"
fi

if [ ! -f "$TRANSITIONS" ]; then
  help "Transition librrary in TRANSITIONS does not exist"
fi

if [ "$NDECOYS" = "" ]; then
  NDECOYS=5
fi

function openswath2json() {
    BASE0=`basename "$1" .mzML.gz`
    BASE1=`basename "$1" .centroid.mzML.gz`
    apython $DIR/openswath2json.py --transitions "$TRANSITIONS" --ndecoys $NDECOYS --outdir "$OUTDIR/$BASE1" \
	                           --chromatograms="${BASE0}_chrom.mzML" --results "${BASE0}_table.tsv" \
                                   --fdr $FDR --score $SCORE
}

for f in "$@"; do
    openswath2json "$f"
done    
