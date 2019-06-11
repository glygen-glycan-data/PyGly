#!/bin/sh
set -x
# Specify whether to work with TEST or DEV
# SMWENV="TEST"
SMWENV=${1:-DEV}
export SMWENV

set -x
./loadsite.py ../wiki
./loadgtc.py ../data/glygen_accessions.txt ../data/glygen_new_accessions.txt
./loadgtc.py ../data/extra_accessions.txt
./loadunicarb.py ../data/uc2gtc.txt ../data/uc2taxa.txt
./loadgtc2pubchem.py ../data/GlyTouCan-PubChem_2019-04-19.csv
./loadglygen.py ../data/glygen_accessions.txt
./loadedlab.py
./loadclassification.py
./loadmonosDB.py ../data/GlyGen_glycans2monodbID.tsv
./loadgdb2gog.py ../data/gdb2gog.txt
./loadsubsump.py ../data/gnome_subsumption_raw.txt ../data/glygen_accessions.txt ../data/glygen_new_accessions.txt ../data/extra_accessions.txt
./loadnames.py EdwardsLab ../data/uckbcomp2glytoucan.txt ../data/shortuckbcomp2glytoucan.txt
./refresh.sh
./refresh.sh -
