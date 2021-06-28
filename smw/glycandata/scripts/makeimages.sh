#!/bin/sh

set -x
Xvfb :1 &
XSCR=$!
export DISPLAY=localhost:1.0
./dumpallseq.py cache images/wurcs
../../../scripts/allimg.sh -N 20 -P 2 -n snfg -d normalinfo -f png -o images/snfg/extended/png images/wurcs
../../../scripts/allimg.sh -N 20 -P 2 -n snfg -d normalinfo -f svg -o images/snfg/extended/svg images/wurcs
./dumpallseq.py cache images/glycoct
../../../scripts/allimg.sh -N 20 -P 2 -n snfg -d normalinfo -f png -o images/snfg/extended/png images/glycoct
../../../scripts/allimg.sh -N 20 -P 2 -n snfg -d normalinfo -f svg -o images/snfg/extended/svg images/glycoct
./dumpallseq.py cache images/genglycoct
../../../scripts/allimg.sh -N 20 -P 2 -n snfg -d normalinfo -f png -o images/snfg/extended/png images/genglycoct
../../../scripts/allimg.sh -N 20 -P 2 -n snfg -d normalinfo -f svg -o images/snfg/extended/svg images/genglycoct
kill $XSCR

rm -f images*.zip images*.tbz*
cd images
tar cf - snfg/extended/png | bzip2 -c | split -b 40m -d - '../images-snfg-extended-png.tbz.'
tar cf - snfg/extended/svg | bzip2 -c | split -b 40m -d - '../images-snfg-extended-svg.tbz.'
echo -e "GlyTouCanAccession\tImage-Size\tImage-Notation\tImage-Style\tImage-Format" > "../images.tsv"
find . \( -name "*.svg" -o -name "*.png" \) -ls | fgrep '/extended/' | awk '{print $7,$11}' | tr '/.' '  ' | awk '{print $5,$1,$2,$3,$4}' | sort -k1,1 -k3,3 -k4,4 -k5,5 | tr ' ' '\t' >> "../images.tsv"
