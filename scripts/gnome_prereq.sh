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

function getclassificationprereq () {
  awk -F'\t' '$2 == "GlycanType" && $3 == "N-linked" && NR > 1 {print $1}' $1 | sort -u > ./GNOme/$2
}

function getglycotreeprereq () {
  zipinfo -1 $1 | awk '{split($0,a,"."); print a[1]}' | sort -u > ./GNOme/$2
}

# Restrictions
getaccprereq export/bcsdb.tsv restrictions/GNOme_BCSDB.accessions.txt
getprereq data/glygen_accessions.txt restrictions/GNOme_GlyGen.accessions.txt
getprereq data/glycosmos_allacc.txt restrictions/GNOme_GlyCosmos.accessions.txt
getclassificationprereq export/classification.tsv restrictions/GNOme_GlyGen_NGlycans.accessions.txt
getclassificationprereq export/classification.tsv  restrictions/GNOme_GlyGen_OGlycans.accessions.txt
getglycotreeprereq export/glycotree_nlinked_gct.zip restrictions/GNOme_GlycoTree_NGlycans.accessions.txt
getglycotreeprereq export/glycotree_olinked_gct.zip restrictions/GNOme_GlycoTree_OGlycans.accessions.txt

#N glycan restriction from alignment to N glycan core motif
curl -s --compressed https://glycomotif.glyomics.org/glycomotif/sparql/query -X POST --data 'query=PREFIX+glycomotif%3A+%3Chttp%3A%2F%2Fglyomics.org%2Fglycomotif%23%3E%0A%0ASELECT+%3Fstructure_gtc_acc%0AWHERE+%7B%0A++++%3Fcollection+a+glycomotif%3ACollection+.%0A++++%3Fcollection+glycomotif%3Aid+%3Fcollectionid+.%0A++++%0A++++%3Fmotif+a+glycomotif%3AMotif+.%0A++++%3Fmotif+glycomotif%3Aaccession+%3Faccession+.%0A++++%3Fmotif+glycomotif%3Aincollection+%3Fcollection+.%0A%0A++++%3Fmotif+glycomotif%3Aglytoucan+%3Fmotif_gtc_acc+.%0A%0A++++%3Falignment+glycomotif%3Amotif_accession+%3Fmotif_gtc_acc+.%0A++++%3Falignment+glycomotif%3Astructure_accession+%3Fstructure_gtc_acc+.%0A%0A++++%3Fmotif+glycomotif%3Aid+%22GGM.001001%22+.%0A%7D' -H 'Accept: application/sparql-results+json' | jq -r '.results.bindings[].structure_gtc_acc.value' | sort -u > ./GNOme/restrictions/GNOme_NGlycans.accessions.txt

# GNOme raw file
getprereq data/gnome_subsumption_raw.txt data/gnome_subsumption_raw.txt

# GlyTouCan retired etc.
getprereq data/glytoucan_archived.txt data/glytoucan_archived.txt
getprereq data/glytoucan_replaced.txt data/glytoucan_replaced.txt

# Synonyms
getprereq data/byonic2glytoucan.txt data/byonic2glytoucan.txt
getprereq data/shortuckbcomp2glytoucan.txt data/shortuckbcomp2glytoucan.txt
getprereq data/shortcomp2glytoucan.txt data/shortcomp2glytoucan.txt

