#!/bin/sh
set -x
DIR=`pwd`
mkdir -p images; 
mkdir -p images/cfg images/cfg/normal images/cfg/compact images/cfg/extended
cat "$@" | (cd images/cfg/normal;  $DIR/getimages.py cfg normal )
cat "$@" | (cd images/cfg/compact; $DIR/getimages.py cfg compact )
cat "$@" | (cd images/cfg/extended; $DIR/getimages.py cfg extended )
cd images
rm -f ../images.zip
zip -q -r "../images.zip" cfg
echo "GlyTouCanAccession" > "../images.txt"
find cfg -name "*.png" \! -empty | sed -n 's/\.png$//p' | sed 's/^.*\///' | sort -u >> ../images.txt
cd ..
echo -e "GlyTouCanAccession\tImage-Size\tImage-CRC\tImage-Notation\tImage-Style" > "images-crc.txt"
unzip -lv images.zip | awk '{print $1,$7,$8}' | tr '/.' '  ' | awk '{print $5,$1,$2,$3,$4}' | awk 'NF == 5' | egrep -v '(00000000|17981711)' | sort -k5,5 -k1,1 | tr ' ' '\t' >> images-crc.txt
$DIR/checkimg.py images-crc.txt
