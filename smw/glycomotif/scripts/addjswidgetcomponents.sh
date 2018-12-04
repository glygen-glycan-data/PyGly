#!/bin/sh

# Specify whether to work with PROD or DEV
SMWENV="TEST"
# SMWENV="PROD"
export SMWENV

# This script requires the glycancomparison.py

set -x
./dumpallseq.py
./addtopology.py
# Add component will require 3 other scripts
./json_substructure.py
./json_topology.py
./mvcomponents.py
