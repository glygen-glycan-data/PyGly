#!/bin/sh
SMW="$1"
GLYCOTREE="$2"
EXPORT="../export"
./glycoTree.sh "$SMW" "$GLYCOTREE"
rm -f $EXPORT/glycotree*.zip
zip -j $EXPORT/glycotree_nlinked_gct.zip $GLYCOTREE/data/glygengct_nlinked/G*.txt $GLYCOTREE/data/glygengct_nlinked/residmap.txt
zip -j $EXPORT/glycotree_olinked_gct.zip $GLYCOTREE/data/glygengct_olinked/G*.txt $GLYCOTREE/data/glygengct_olinked/residmap.txt
