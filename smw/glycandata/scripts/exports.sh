#!/bin/sh

# Exports computed directly from glycandata (or its filesystem cache)
# because the triple store queries take too long (or are too difficult
# to implement)

# ./exports.sh cache 
set -x 
rm -rf glycandata*.tdb
./tordf.py "$@" glycandata | gzip -9 -c > ../export/glycandata.rdf.gz
rm -f ../export/glycandata.rdf.gz.[0-9][0-9]
split -b 40m -d ../export/glycandata.rdf.gz ../export/glycandata.rdf.gz.
./loadts.sh ../export/glycandata.rdf.gz
./queryall.sh glycandata.tdb
rm -f ../export/glycoctxml.zip.[0-9][0-9]
split -b 40m -d ../export/glycoctxml.zip ../export/glycoctxml.zip.
./glycanprop.py "$@" > ../export/glycan_properties.tsv
./monocomp.py "$@" > ../export/monocomp.tsv
./subsumption.py "$@" > ../export/subsumption.tsv
./motif.py "$@" allmotifs > ../export/allmotifs.tsv
./motif.py "$@" allaligns > ../export/allmotifaligns.tsv
./byonic_database.py "$@" > ../export/byonic_glygen_human_nlinked.txt
./species.py "$@" > ../export/species_expanded.tsv
./glycoTree_exports.sh "$@" ../glycoTree
python3 ./dumpiupacsyms.py > ../export/iupac_syms.tsv
./glygen_retired_accessions.sh
