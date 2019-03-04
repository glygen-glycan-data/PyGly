#!/bin/sh

# Specify whether to work with TEST or DEV
# SMWENV="TEST"
SMWENV="DEV"
export SMWENV

set -x
./loadsite.py ../wiki
./loadgtc.py ../data/glytoucan_accessions.txt
./refresh.sh
./refresh.sh -
