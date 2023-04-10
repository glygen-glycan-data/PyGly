#!/bin/sh

set -x

#N-linked
rm -f glycotree_nlinked_gct.zip
wget https://raw.githubusercontent.com/glygen-glycan-data/PyGly/GlyGen-GlycanData-Export-Current/smw/glycandata/export/glycotree_nlinked_gct.zip
rm -f nlinked/*.txt
cd nlinked
unzip ../glycotree_nlinked_gct.zip
cd ..

wget -O - -q 'https://sandbox.glyomics.org/api/list.php?par=&mode=all_N' | jq -r '.[].glytoucan_ac' > all_N.txt
wget -O - -q 'https://sandbox.glyomics.org/api/list.php?par=&mode=mapped_N' | jq -r '.[].glytoucan_ac' > mapped_N.txt

rm -f nlinked/*.json
for acc in `cat mapped_N.txt`; do
  wget -O "nlinked/$acc.json" -q "https://sandbox.glyomics.org/api/glycan-v5.php/$acc"
done

find nlinked -name "*.txt" | fgrep -v -f mapped_N.txt | xargs -n 10 rm -f

#O-linked
rm -f glycotree_olinked_gct.zip
wget https://raw.githubusercontent.com/glygen-glycan-data/PyGly/GlyGen-GlycanData-Export-Current/smw/glycandata/export/glycotree_olinked_gct.zip
rm -f olinked/*.txt
cd olinked
unzip ../glycotree_olinked_gct.zip
cd ..

wget -O - -q 'https://sandbox.glyomics.org/api/list.php?par=&mode=all_O' | jq -r '.[].glytoucan_ac' > all_O.txt
wget -O - -q 'https://sandbox.glyomics.org/api/list.php?par=&mode=mapped_O' | jq -r '.[].glytoucan_ac' > mapped_O.txt

rm -f olinked/*.json
for acc in `cat mapped_O.txt`; do
  wget -O "olinked/$acc.json" -q "https://sandbox.glyomics.org/api/glycan-v5.php/$acc"
done

find olinked -name "*.txt" | fgrep -v -f mapped_O.txt | xargs -n 10 rm -f
