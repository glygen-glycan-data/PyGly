#!/bin/sh
set -x

rm -f images*.zip images*.tbz*
mkdir -p images/snfg/extended/png
cat ../data/{glygen,glygen_new,extra}_accessions.txt | awk '{printf("/data/projects/GlyGen/APIFramework/src/Application/Glymage/work/snfg/extended/%s.png\n",$1)}' | xargs -i cp {} images/snfg/extended/png
mkdir -p images/snfg/extended/svg
cat ../data/{glygen,glygen_new,extra}_accessions.txt | awk '{printf("/data/projects/GlyGen/APIFramework/src/Application/Glymage/work/snfg/extended/%s.svg\n",$1)}' | xargs -i cp {} images/snfg/extended/svg
mkdir -p images/snfg/extended/json
cat ../data/{glygen,glygen_new,extra}_accessions.txt | awk '{printf("/data/projects/GlyGen/APIFramework/src/Application/Glymage/work/snfg/extended/%s.json\n",$1)}' | xargs -i cp {} images/snfg/extended/json
pwd

cd images
tar cf - snfg/extended/png | bzip2 -c | split -b 40m -d - '../images-snfg-extended-png.tbz.'
tar cf - snfg/extended/svg | bzip2 -c | split -b 40m -d - '../images-snfg-extended-svg.tbz.'
tar cf - snfg/extended/json | bzip2 -c | split -b 40m -d - '../images-snfg-extended-json.tbz.'
echo -e "GlyTouCanAccession\tImage-Size\tImage-Notation\tImage-Style\tImage-Format" > "../images.tsv"
find snfg/extended \( -name "*.svg" -o -name "*.png"  -o -name "*.json" \) -ls | fgrep '/extended/' | awk '{print $7,$11}' | tr '/.' '  ' | awk '{print $5,$1,$2,$3,$4}' | sort -k1,1 -k3,3 -k4,4 -k5,5 | tr ' ' '\t' >> "../images.tsv"
