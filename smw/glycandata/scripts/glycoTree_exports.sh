#!/bin/sh
GLYCOTREE="$1"
EXPORT="../export"
set -euxo pipefail
cp $GLYCOTREE/accessions.lst $EXPORT/glycotree_accessions.tsv
cp $GLYCOTREE/glycotree_annotated_glycans.tsv.gz $EXPORT
gunzip -c $GLYCOTREE/glycotree_glycan_caveats.tsv.gz > $EXPORT/glycotree_glycan_caveats.tsv
cp $GLYCOTREE/portal/api/paths/allPaths.json.gz $EXPORT/glycotree_allpaths.json.gz
