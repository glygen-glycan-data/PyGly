#!/bin/sh
GLYCOTREE="$1"
EXPORT="../export"
set -x
cp $GLYCOTREE/accessions.lst $EXPORT/glycotree_accessions.tsv
cp $GLYCOTREE/glycotree_annotated_glycans.tsv.gz $EXPORT
gunzip -c $GLYCOTREE/glycotree_glycan_caveats.tsv.gz > $EXPORT/glycotree_glycan_caveats.tsv
gzip -c $GLYCOTREE/portal/api/paths/allPaths.json > $EXPORT/glycotree_allpaths.json.gz
