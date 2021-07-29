#!/bin/sh

set -x

if [[ ! ("$(pwd)" =~ "/PyGly/scripts") ]]; then
  exit 1
fi

if [ ! -d ./GNOme ]; then
  exit 1
fi

set -x

# Force us to get from the github master branch (not locally!) to standard out
# We are running under a cloned copy of PyGly repository, which also
# contains smw/glycandata, so we can use git show.
GET="git show origin/master:smw/glycandata/"

function get() {
  ${GET}$1
}

function getprereq () {
  get $1 > ./GNOme/$2
}

function getaccprereq () {
  get $1 | awk 'NR > 1 {print $1}' | sort -u > ./GNOme/$2
}

# Restrictions
getaccprereq export/bcsdb.tsv restrictions/GNOme_BCSDB.accessions.txt
getprereq data/glygen_accessions.txt restrictions/GNOme_GlyGen.accessions.txt
getprereq data/glycosmos_allacc.txt restrictions/GNOme_GlyCosmos.accessions.txt

# GNOme raw file
getprereq data/gnome_subsumption_raw.txt data/gnome_subsumption_raw.txt

# GlyTouCan retired etc.
getprereq data/glytoucan_archived.txt data/glytoucan_archived.txt
getprereq data/glytoucan_replaced.txt data/glytoucan_replaced.txt

# Synonyms
getprereq data/byonic2glytoucan.txt data/byonic2glytoucan.txt
getprereq data/shortuckbcomp2glytoucan.txt data/shortuckbcomp2glytoucan.txt
getprereq data/shortcomp2glytoucan.txt data/shortcomp2glytoucan.txt

