#!/bin/sh
set -x
for q in ../queries/*.sparql ; do
  e=`echo "$q" | sed -e 's/queries/export/' -e 's/sparql$/tsv/'`
  ./queryts.sh "$1" "$q" | ./queryts2tsv.py > "$e"
done
for q in ../queries/*.sparql2zip ; do
  e=`echo "$q" | sed -e 's/queries/export/' -e 's/sparql2zip$/tsv/'`
  f=`echo "$q" | sed -e 's/queries/export/' -e 's/sparql2zip$/zip/'`
  ./queryts.sh "$1" "$q" | ./queryts2seqzf.py "$f" > "$e"
done
