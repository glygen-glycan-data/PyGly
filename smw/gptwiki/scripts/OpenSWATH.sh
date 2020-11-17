#!/bin/sh

help() {
  if [ "$1" != "" ]; then
    echo "" 1>&2
    echo "$1" 1>&2
  fi
  cat <<EOF | fmt -w80 -s 1>&2

OpenSWATH.sh [ options ] <params> <spectra>.mzML.gz ...

Options:
  -c <lc-calibration>   NRT Peptides (NRT), Endogenous
                        Peptides (END), Explicit Slope and
                        and Intercept (ESI), None (NONE).
                        Default: NRT.
  -h	Help text.

Parameter file sets the follwing variables:

  TRANSITIONS=<Transition-Library-File-With-Decoys>
  NDECOYS=<Number-of-Decoy-Replicates>
  WINDOWS=<SWATH-Precursor-Windows-File>
  NRTTRANSITIONS=<Transition-Library-File> (if using NRT calibration)
  TRANSFORMATION=<TrafoXML-File> (if using ESI calibration)
  RTWINDOW=<2*NRT-Tolerance> (if not using NONE calibration, default: 840)

EOF
  exit 1;
}

DIR=`dirname $0`
DIR=`readlink -f $DIR`

LCCAL="NRT"
VERBOSE=0
OUTDIR="."
while getopts ":c:o:vh" o ; do
        case $o in
	        c ) LCCAL="$OPTARG";;
	        o ) OUTDIR="$OPTARG";;
	        v ) VERBOSE=1;;
                h ) help "";;
                * ) help "Invalid option: -$OPTARG"
        esac
done

DONRT=0; DOEND=0; DOESI=0;
case $LCCAL in 
  NRT) DONRT=1;;
  END) DOEND=1;;
  ESI) DOESI=1;;
  NONE) DONONE=1;;
  *) help "Bad calibration (-c) parameter: $LCCAL";;
esac

shift $(($OPTIND - 1)) 

if [ "$1" = "" ]; then
  help "Parameter file not provided on the command-line"
fi

PARAMS="$1"
PARAMS=`readlink -f "$PARAMS"`
. $PARAMS

shift

if [ $DONRT = 1 -a "$NRTTRANSITIONS" = "" ]; then
    help "NRT Peptide calibration: NRTTRANSITIONS missing from parameter file $1"
fi

if [ $DOESI = 1 -a "$TRANSFORMATION" = "" ]; then
    help "Explicit Slope and Intercept calibration: TRANSFORMATION missing from parameter file $1"
fi

if [ "$TRANSITIONS" = "" ]; then
  help "TRANSITIONS missing from parameter file $1"
fi

if [ ! -f "$TRANSITIONS" ]; then
  help "Transition librrary in TRANSITIONS does not exist"
fi

if [ "$WINDOWS" = "" ]; then
  help "WINDOWS missing from parameter file $1"
fi

if [ ! -f "$WINDOWS" ]; then
  help "SWATH windows file in WINDOWS does not exist"
fi

if [ "$RTWINDOW" = "" ]; then
  RTWINDOW=840
fi

if [ "$NDECOYS" = "" ]; then
  NDECOYS=5
fi

OSW="/tools/openswath/bin/OpenSwathWorkflow.sh"

function openswath() {
    FILE="$1"
    shift
    BASE=`basename "$FILE" .mzML.gz`
    TMP=`mktemp -d -p /scratch -t openswath.XXXXXX`
    echo $OSW \
	-in "$FILE" \
	-tr "$TRANSITIONS" \
	-sort_swath_maps \
	-readOptions cache \
	-threads 4 \
	-tempDirectory "$TMP/" \
	-batchSize 1000 \
	-extra_rt_extraction_window 10000000 \
	-swath_windows_file "$WINDOWS" \
	-out_tsv "${OUTDIR}/${BASE}_table.tsv" \
	-out_chrom "${OUTDIR}/${BASE}_chrom.mzML" \
	"$@" \
	>> "${OUTDIR}/${BASE}.log"
    $OSW \
	-in "$FILE" \
	-tr "$TRANSITIONS" \
	-sort_swath_maps \
	-readOptions cache \
	-threads 4 \
	-tempDirectory "$TMP/" \
	-batchSize 1000 \
	-extra_rt_extraction_window 10000000 \
	-swath_windows_file "$WINDOWS" \
	-out_tsv "${OUTDIR}/${BASE}_table.tsv" \
	-out_chrom "${OUTDIR}/${BASE}_chrom.mzML" \
	"$@" \
	>> "${OUTDIR}/${BASE}.log" 2>&1
    rm -rf $TMP*
}

function openswath_nort() {
    openswath "$@" -Scoring:Scores:use_rt_score false -rt_extraction_window -1
}

function openswath_nrt() {
    openswath "$@" -rt_extraction_window "$RTWINDOW" -tr_irt "$NRTTRANSITIONS" 
}

function openswath_esi() {
    openswath "$@" -rt_extraction_window "$RTWINDOW" -rt_norm "$TRAFOXML"
}

function getcoeff() {
    apython $DIR/getoswcoeffs.py "$@"
}

function fdr() {
    apython $DIR/../analysis/fdr.py "$@"
}

function openswath_cal() {
    if [ $DONRT -eq 1 ]; then
	openswath_nrt "$@"
    elif [ $DOESI -eq 1 ]; then
	TRAFOXML="${OUTDIR}/${BASE}_cal.trafoXML"
	cp -f "$TRANSFORMATION" "$TRAFOXML"
	openswath_esi "$@"
    elif [ $DOEND -eq 1 ]; then
	openswath_nort "$@"
	BASE=`basename "$1" .mzML.gz`
	TABLE="${OUTDIR}/${BASE}_table.tsv"
	TRAFOXML="${OUTDIR}/${BASE}_cal.trafoXML"
        echo -e "Score\tFDR\tTarget\tDecoy" >> "${OUTDIR}/${BASE}.log"
	fdr "$TABLE" "$NDECOYS" | awk '$2 < 0.05' >> "${OUTDIR}/${BASE}.log" 2>&1
	( getcoeff "$TABLE" "$TRANSITIONS" "$NDECOYS" > "$TRAFOXML" ) 2>> "${OUTDIR}/${BASE}.log"
	if [ -s "$TRAFOXML" ]; then
	    openswath_esi "$@"
        fi
    elif [ $DONONE -eq 1 ]; then
	openswath_nort "$@"
    fi
}

for f in "$@"; do
    BASE=`basename "$f" .mzML.gz`
    LOG="${OUTDIR}/${BASE}.log"
    echo "" > "$LOG"
    openswath_cal "$f"
    TABLE="${OUTDIR}/${BASE}_table.tsv"
    echo -e "Score\tFDR\tTarget\tDecoy" >> "$LOG"
    fdr "$TABLE" "$NDECOYS" | awk '$2 < 0.01' >> "$LOG"
done    
