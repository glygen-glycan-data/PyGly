#!/bin/sh
set -x
# Specify whether to work with TEST or DEV
# SMWENV="TEST"
SMWENV=${1:-DEV}
export SMWENV

set -x
./loadsite.py ../wiki
./loadgtc.py ../data/glygen_accessions.txt
./loadgtc.py ../data/extra_accessions.txt
./loadglygen.py ../data/glygen_accessions.txt
./loadedlab.py
./loadclassification.py
./loadmonosDB.py ../data/GlyGen_glycans2monodbID.tsv
./loadsubsump.py ../data/gnome_subsumption_raw.txt ../data/glygen_accessions.txt ../data/extra_accessions.txt
./loadnames.py EdwardsLab ../data/uckbcomp2glytoucan.txt ../data/shortuckbcomp2glytoucan.txt
./refresh.sh
./refresh.sh -
