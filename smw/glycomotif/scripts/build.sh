#!/bin/sh

# Specify whether to work with PROD or DEV
SMWENV="DEV"
# SMWENV="PROD"
export SMWENV

set -x
./loadsite.py ../wiki
./clearmotifs.py
./makecoll.py
./loadgtcmotif.py
./loadccrcmotif.py ../data/MotifsMP2.xlsx
# ./dlglycoepitope.py > ../data/glycoepitope.txt
./loadglycoepimotif.py ../data/glycoepitope.txt
./loadoriginalglydinmotif.py
./loadglydinmotif.py
./loadallmotif.py
./addsameas.py