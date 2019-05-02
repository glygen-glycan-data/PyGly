#!/bin/sh

# Specify whether to work with TEST or DEV
# SMWENV="TEST"
SMWENV="DEV"
export SMWENV

set -x
./loadsite.py ../wiki
./loadgtc.py ../data/glygen_accessions.txt
./loadgtc.py ../data/extra_accessions.txt
./loadglygen.py ../data/glygen_accessions.txt
./loadedlab.py
./loadclassification.py
./loadmonosDB.py ../data/GlyGen_glycans2monodbID.tsv
./refresh.sh
./refresh.sh -
