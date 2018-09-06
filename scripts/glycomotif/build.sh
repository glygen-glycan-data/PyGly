#!/bin/sh

./clearmotifs.py
./makecoll.py
./loadgtcmotif.py
./loadccrcmotif.py data/MotifsMP2.xlsx
# ./dlglycoepitope.py > data/glycoepitope.txt
./loadglycoepitopemotif.py data/glycoepitope.txt
./loadglygenmotif.py
./addsameas.py
