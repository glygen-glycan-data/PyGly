#!/bin/sh

# Exports computed directly from glycandata (or its filesystem cache)
# because the triple store queries take too long (or are too difficult
# to implement)

# ./exports.sh cache 
set -x 
./glycanprop.py "$1" > ../exports/glycan_properties.tsv
./monocomp.py "$1" > ../exports/monocomp.tsv
./subsumption.py "$1" > ../exports/subsumption.tsv
./motif.py "$1" allmotifs > ../exports/allmotifs.tsv
./motif.py "$1" allaligns > ../exports/allmotifaligns.tsv
