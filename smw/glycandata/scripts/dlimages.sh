#!/bin/sh
set -x
DIR=`pwd`
mkdir -p images/cfg/extended/png images/snfg/extended/png images/cfg/extended/svg images/snfg/extended/svg
find images -name ".gtccache*" -exec rm -f {} \;
cat "$@" | (cd images/cfg/extended/png; $DIR/getimages.py cfg extended png )
cat "$@" | (cd images/snfg/extended/png; $DIR/getimages.py snfg extended png )
cat "$@" | (cd images/cfg/extended/svg; $DIR/getimages.py cfg extended svg )
cat "$@" | (cd images/snfg/extended/svg; $DIR/getimages.py snfg extended svg )
find images -name ".gtccache*" -exec rm -f {} \;
rm -f images*.zip images*.tbz*
cd images
tar cf - cfg/extended/png | bzip2 -c | split -b 40m -d - '../images-cfg-extended-png.tbz.'
tar cf - snfg/extended/png | bzip2 -c | split -b 40m -d - '../images-snfg-extended-png.tbz.'
tar cf - cfg/extended/svg | bzip2 -c | split -b 40m -d - '../images-cfg-extended-svg.tbz.'
tar cf - snfg/extended/svg | bzip2 -c | split -b 40m -d - '../images-snfg-extended-svg.tbz.'
echo -e "GlyTouCanAccession\tImage-Size\tImage-Notation\tImage-Style\tImage-Format" > "../images.tsv"
find . \( -name "*.svg" -o -name "*.png" \) -ls | fgrep '/extended/' | awk '{print $7,$11}' | tr '/.' '  ' | awk '{print $5,$1,$2,$3,$4}' | sort -k1,1 -k3,3 -k4,4 -k5,5 | tr ' ' '\t' >> "../images.tsv"
