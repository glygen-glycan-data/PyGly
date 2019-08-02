#!/bin/sh
set -x
# Specify whether to work with TEST or DEV
# SMWENV="TEST"
SMWENV=${1:-DEV}
export SMWENV

CACHE=cache

if [ ! -d $CACHE ]; then
  ./todiskcache.py $CACHE
fi

# Load/Modify Disk cache...
./loadgtc.py $CACHE ../data/glygen_accessions.txt ../data/glygen_new_accessions.txt
./loadgtc.py $CACHE ../data/extra_accessions.txt
./loadunicarb.py $CACHE ../data/uc2gtc.txt ../data/uc2taxa.txt ../data/uc2pubmed.txt
./loadgtc2pubchem.py $CACHE ../data/GlyTouCan-PubChem_2019-04-19.csv
./loadglygen.py $CACHE ../data/glygen_accessions.txt
./loadedlab.py $CACHE 
./loadclassification.py $CACHE
./loadmonosDB.py $CACHE ../data/GlyGen_glycans2monodbID.tsv
./loadgdb2gog.py $CACHE ../data/gdb2gog.txt
./loadcanonres.py $CACHE ../data/canonres.csv ../data/canonres2gtc.csv
./loadsubsump.py $CACHE ../data/gnome_subsumption_raw.txt ../data/glygen_accessions.txt ../data/glygen_new_accessions.txt ../data/extra_accessions.txt
./loadspecies.py $CACHE ../data/gnome_subsumption_raw.txt
./loadnames.py $CACHE EdwardsLab ../data/uckbcomp2glytoucan.txt ../data/shortuckbcomp2glytoucan.txt

# Load to wiki
./loadsite.py --smwenv $SMWENV ../wiki
./tosite.py $CACHE
./refresh.py 
./refresh.py -
