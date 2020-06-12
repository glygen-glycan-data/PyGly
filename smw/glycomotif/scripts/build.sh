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
./loadglygenmotif.py ../data/GlyGen-Motif*.tsv
./loadallmotif.py
./addsameas.py
# ./loadmotifalign.py ../data/motif_alignment.tsv
./widgetdata.py ../data/topology.json ../data/nonredonly.json ../data/redonly.json
echo "please commit the new json files to GitHub"
./loadwidgetbool.py
./addtopology.py
./refresh.py
./refresh.py -

