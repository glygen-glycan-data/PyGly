#!/bin/sh

set -x

if [[ ! ("$(pwd)" =~ "/PyGly/scripts") ]]; then
  exit 1
fi

if [ ! -d ./GNOme ]; then
  exit 1
fi

set -x

GLYRES="python2 ../pygly/GlycanResource/main.py"

function restriction () {
  $GLYRES $1 $2 | sort -u > ./GNOme/restrictions/GNOme_${1}.accessions.txt 
}

function getdata () {
  $GLYRES $1 $2 | bash -c '(read -r; echo "$REPLY"; sort)' > ./GNOme/data/$3
}

# Restrictions
restriction GlyGen allglycans
restriction GlyCosmosNoCache allaccessions

# GlyTouCan retired etc.
getdata GlyCosmosNoCache archived glytoucan_archived.txt
getdata GlyCosmosNoCache replaced glytoucan_archived.txt

# Synonyms
python2 ../smw/glycandata/scripts/getbasecomplist.py \* > ./GNOme/data/basecomplist.txt
( cd ./GNOme/data; ../../../smw/glycandata/data/splitbasecomp.sh )
