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
./loadgtc.py $CACHE ../data/glygen_accessions.txt ../data/glygen_new_accessions.txt ../data/extra_accessions.txt
./loadunicarb.py $CACHE ../data/uc2gtc.txt ../data/uc2pubmed.txt ../data/uckbcomp2glytoucan.txt
./unicarbkb_taxid.py ../data/uc2gtc.txt ../data/uc2taxa.txt ../data/uckbcomp2glytoucan.txt > ../data/unicarbkb_taxid.txt
./glygen_taxid.py --unicarbkb "../export/unicarbkb.tsv" --glyconnect ../data/glyconnect2glytoucan.txt > ../data/glygen_taxid.txt
./loadtaxid.py $CACHE ../data/glygen*_taxid.txt ../data/unicarbkb_taxid.txt ../data/gptwiki_taxid.txt 
./loadgtc2pubchem.py $CACHE ../data/GlyTouCan-PubChem_2020-04-08.csv
./loadgtc2chebi.py $CACHE ../data/GlyTouCan-ChEBI_2019-08-23.tsv
./loadgtc2psimod.py $CACHE ../data/psimod2glytoucan.txt
./loadglygen.py $CACHE ../data/glygen_accessions.txt
./loadedlab.py $CACHE 
./loadclassification.py $CACHE
./loadmonosDB.py $CACHE ../data/allmonosaccharideDB.tsv
# ./loadgdb2gog.py $CACHE ../data/gdb2gog.txt
./loadcanonres.py $CACHE ../data/canonres.csv ../data/canonres2gtc.csv
./loadsubsump.py $CACHE ../data/gnome_subsumption_raw.txt ../data/glygen_accessions.txt ../data/glygen_new_accessions.txt ../data/extra_accessions.txt
./loadspecies.py $CACHE ../data/gnome_subsumption_raw.txt
./loadnames.py $CACHE UniCarbKB EdwardsLab ../data/uckbcomp2glytoucan.txt 
./loadnames.py $CACHE ShortUniCarbKB EdwardsLab ../data/shortuckbcomp2glytoucan.txt
./loadnames.py $CACHE Byonic EdwardsLab ../data/byonic2glytoucan.txt
./loadnames.py $CACHE ShortComp EdwardsLab ../data/shortcomp2glytoucan.txt

# Load to wiki
./loadsite.py --smwenv $SMWENV ../wiki
./tosite.py $CACHE
./refresh.py 
./refresh.py -
