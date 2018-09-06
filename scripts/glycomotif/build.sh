#!/bin/sh

./clearmotifs.py
./makecoll.py
./loadgtcmotif.py
./loadccrcmotif.py MotifsMP2.xlsx
./loadglygenmotif.py
./addsameas.py
