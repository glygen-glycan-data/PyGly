#!/bin/sh

set -x

if [[ ! ("$(pwd)" =~ "/PyGly/scripts") ]]; then
  exit 1
fi

if [ ! -d ./GNOme ]; then
  exit 1
fi

set -x

GLYRES="./glyres.py"

function restriction () {
  $GLYRES $1 $2 | sort -u > ./GNOme/restrictions/GNOme_${3}.accessions.txt 
}

function glygentyperes () {
  $GLYRES GlyGen glycans_bytype $1 | sort -u > ./GNOme/restrictions/GNOme_${2}.accessions.txt 
}

function sandboxres () {
  $GLYRES GlycoTreeSandbox list $1 | sort -u > ./GNOme/restrictions/GNOme_${2}.accessions.txt 
}

function glycomotifres () {
  $GLYRES GlycoMotifNoCache getstruct GGM $1 | awk '{print $1}' | sort -u > ./GNOme/restrictions/GNOme_${2}.accessions.txt 
}

function getdata () {
  $GLYRES $1 $2 | bash -c '(read -r; echo "$REPLY"; sort)' > ./GNOme/data/$3
}

# Restrictions
restriction GlyGen allglycans GlyGen
restriction GlyCosmosNoCache allaccessions GlyCosmos

glygentyperes N-linked GlyGen_NGlycans
glygentyperes O-linked GlyGen_OGlycans

sandboxres mapped_N GlycoTree_NGlycans
sandboxres mapped_O GlycoTree_OGlycans

glycomotifres 001001 NGlycans

# GlyTouCan retired etc.
getdata GlyCosmosNoCache archived glytoucan_archived.txt
getdata GlyCosmosNoCache replaced glytoucan_replaced.txt

# Synonyms
python3 ../smw/glycandata/scripts/getbasecomplist.py \* > ./GNOme/data/basecomplist.txt
( cd ./GNOme/data; ../../../smw/glycandata/data/splitbasecomp.sh )
