#!/bin/sh

# Specify whether to work with PROD or DEV
# SMWENV="TEST"
# SMWENV="PROD"
SMWENV="DEV"
export SMWENV

set -x
./loadsite.py ../wiki
./makecoll.py
./loadgtcmotif.py
./loadccrcmotif.py ../data/MotifsMP2.xlsx
# ./dlglycoepitope.py > ../data/glycoepitope.txt
./loadglycoepimotif.py ../data/glycoepitope.txt
./loadglydinmotif.py ../data/epitopes.xlsx
./loadunicarbmotif.py ../data/Unicarb.xlsx
./loadallmotif.py
./addsameas.py
rm -rf ./dumps
./dumpallseq.py
./json_substructure.py
./json_topology.py
./loadwidgetbool.py
./addtopology.py
./refresh.py
./refresh.py -
