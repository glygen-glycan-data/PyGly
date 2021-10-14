#!/bin/sh
SMW="$1"
GLYCOTREE="$2"
EXPORT="../export"
# head -n 1 $GLYCOTREE/model/annotated_glycans.csv | tr ',' '\t' > $EXPORT/glycotree_annotated_glycans.tsv
# fgrep -v -w glytoucan_ac $GLYCOTREE/model/annotated_glycans.csv | tr ',' '\t' | sort -k1,1 >> $EXPORT/glycotree_annotated_glycans.tsv
# echo "GlyTouCanAccession" > $EXPORT/glycotree_accessions.tsv
# awk 'NR > 1 {print $1}' $EXPORT/glycotree_annotated_glycans.tsv | sort -u >> $EXPORT/glycotree_accessions.tsv
./glycoTree.sh "$SMW" "$GLYCOTREE"
rm -f $EXPORT/glycotree_gct.zip $EXPORT/glycotree_svg.zip
zip -j $EXPORT/glycotree_gct.zip $GLYCOTREE/data/glygengct/G*.txt
# zip -j $EXPORT/glycotree_svg.zip $GLYCOTREE/data/glygensvg/G*.svg
