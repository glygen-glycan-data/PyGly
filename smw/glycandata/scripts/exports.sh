#!/bin/sh

# Exports computed directly from glycandata (or its filesystem cache)
# because the triple store queries take too long (or are too difficult
# to implement)

# ./exports.sh cache 
set -x 
./glycanprop.py "$@" > ../export/glycan_properties.tsv
./monocomp.py "$@" > ../export/monocomp.tsv
./subsumption.py "$@" > ../export/subsumption.tsv
./motif.py "$@" allmotifs > ../export/allmotifs.tsv
./motif.py "$@" allaligns > ../export/allmotifaligns.tsv
./byonic_database.py "$@" > ../export/byonic_glygen_human_nlinked.txt
./species.py "$@" > ../export/species_expanded.tsv
./glycoTree_exports.sh ../glycoTree
./glygen_retired_accessions.sh
