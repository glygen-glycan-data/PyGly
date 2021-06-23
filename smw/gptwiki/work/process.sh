#!/bin/sh

set -x

DIR=$HOME/projects/PyGly
DIR1=/data/projects/GlyGen/PyGly/smw/gptwiki/scripts
DATADIR=/data/projects/GlyGen/PyGly/smw/gptwiki/data

. ./params.txt

SAMPLE=`printf 'SA%06d' $SAMPLE`
METHOD=`printf 'ME%06d' $METHOD`
PROTDB="$DATADIR/${PROTDB:-human}.fasta" 
GLYDB="$DATADIR/${GLYDB:-MSCompDB}.gct" 
PRECTOL="${PRECTOL:-0.05} Da"
FRAGTOL="${FRAGTOL:-0.05} Da"
if [ "$SPLIT" -eq 1 ]; then
  SPLIT="--splitresults"
else
  SPLIT=""
fi

echo "PARAMETERS:"
echo "SAMPLE: $SAMPLE"
echo "METHOD: $METHOD"
echo "PROTDB: $PROTDB"
echo "GLYDB: $GLYDB"
echo "PRECTOL: $PRECTOL"
echo "FRAGTOL: $FRAGTOL"
echo "SPLIT: $SPLIT"
echo "EXTRAARGS: $EXTRAARGS"

function echocmd {
  echo "$@"
  $@
}
                                                                                                                            
for SPEC in $@; do
  if [ ! -f "$SPEC" ]; then
    continue
  fi
  case "$SPEC" in 
    *.mzML.gz) BASE=`basename $SPEC .mzML.gz`;;
    *.msp) BASE=`basename $SPEC .msp`;;
  esac
  rm -f $BASE.SA*.ME*.*.txt
  # $BASE.MSFragger-*.tsv 
  for RESULT in `ls $BASE.Byonic-*.xlsx $BASE.GPS-*.tsv $BASE.msp 2>/dev/null`; do
    if [ -f "$RESULT" ]; then
      SEARCHID=`echo "$RESULT" | sed 's/^.*\.\([^.]*\)\.[^.]*$/\1/'`
      TRANS="$BASE.$SAMPLE.$METHOD.$SEARCHID.txt"
      echo "Extracting transitions from $RESULT"
      case "$RESULT" in 
        *.GPS-*.tsv) 
          $DIR/GlycoPeptideTransitions.py -s "$SPEC" -r "$RESULT" -g "$GLYDB" -p "$PROTDB" --format GPS -e "$FRAGTOL" -E "$PRECTOL" $SPLIT $EXTRAARGS > "$TRANS" ;;
        *.MSFragger-*.tsv) 
          $DIR/GlycoPeptideTransitions.py -s "$SPEC" -r "$RESULT" -g "$GLYDB" -p "$PROTDB" --format MSFragger -e "$FRAGTOL" -E "$PRECTOL" $SPLIT $EXTRAARGS > "$TRANS" ;;
        *.Byonic-*.xlsx) 
          $DIR/GlycoPeptideTransitions.py -s "$SPEC" -r "$RESULT" -g "$GLYDB" -p "$PROTDB" --format Byonic -e "$FRAGTOL" -E "$PRECTOL" $SPLIT $EXTRAARGS > "$TRANS" ;;
        *.msp) 
          TRANS="$BASE.$SAMPLE.$METHOD.MSP.txt";
          $DIR/GlycoPeptideTransitions.py -s "$RESULT" -r "$RESULT" -g "$GLYDB" -p "$PROTDB" --format MSP -e "$FRAGTOL" -E "$PRECTOL" $SPLIT $EXTRAARGS > "$TRANS" ;;
      esac
    fi
  done
  rm -f $BASE.*.*.merge.txt
  TRANS="$BASE.$SAMPLE.$METHOD.merge.txt"
  echo "Merging transitions for $BASE"  
  $DIR/GlycoPeptideTransitions.py -s "$SPEC" -r "`ls -1 $BASE.$SAMPLE.$METHOD.*.txt | fgrep -v .merge.txt | tr '\n' ' ' | sed 's/  *$//' `" -g "$GLYDB" --format GPT -e "$FRAGTOL" -E "$PRECTOL" > "$TRANS"
  if [ -d "$BASE" ]; then
      rm -rf "$BASE"
  fi
  awk 'NR > 1 {print $13}' $TRANS | tr '(;' ' \n' | awk '{print $1}' | sort -n -u | $DIR/JSONSpectraExplode.py "$SPEC"
  $DIR/JSONGlycopeptideFragments.py -g "$GLYDB" -p "$TRANS" -d $BASE
  case $SPEC in 
    *.mzML.gz)
      awk 'NR > 1 {print $2,$7,$8}' $TRANS | sort -u | $DIR/JSONXIC.py -s "$SPEC" -e 0.005;
      awk 'NR > 1 {print $2,$7,$8}' $TRANS | sort -u | $DIR/JSONXIC.py -s "$SPEC" -e 0.01;
      awk 'NR > 1 {print $2,$7,$8}' $TRANS | sort -u | $DIR/JSONXIC.py -s "$SPEC" -e 0.05;
      awk 'NR > 1 {print $2,$7,$8}' $TRANS | sort -u | $DIR/JSONXIC.py -s "$SPEC" -e 0.1;
      $DIR1/fitlcpeakdda.py $BASE.$SAMPLE.$METHOD.merge.txt
      $DIR/XIC-iRT.py --spectrum "$SPEC" --irt $DATADIR/irtpeptides.txt -e 0.005 --matchedpeaks 2 > ${BASE}/${BASE}.irt.5.txt;
      $DIR1/irtregress.py $DATADIR/irtpeptides.txt ${BASE}/${BASE}.irt.5.txt > ${BASE}/${BASE}.irt.5.fit;
      $DIR/XIC-iRT.py --spectrum "$SPEC" --irt $DATADIR/irtpeptides.txt -e 0.01 --matchedpeaks 2 > ${BASE}/${BASE}.irt.10.txt;
      $DIR1/irtregress.py $DATADIR/irtpeptides.txt ${BASE}/${BASE}.irt.10.txt > ${BASE}/${BASE}.irt.10.fit;
      $DIR/XIC-iRT.py --spectrum "$SPEC" --irt $DATADIR/irtpeptides.txt -e 0.05 --matchedpeaks 2 > ${BASE}/${BASE}.irt.50.txt;
      $DIR1/irtregress.py $DATADIR/irtpeptides.txt ${BASE}/${BASE}.irt.50.txt > ${BASE}/${BASE}.irt.50.fit;
       ;;
    *.msp) 
       ;;                                                                                                                   
  esac
  cp "$SPEC" "$BASE"
done
