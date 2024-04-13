#!/bin/sh

set -x

restriction_set_names_standard=(
  "BCSDB"
  "GlyGen"
  "GlyCosmos"
  "GlyGen_NGlycans"
  "GlyGen_OGlycans"
  "NGlycans"
)

restriction_set_names_ancestor=(
  "GlycoTree_NGlycans"
  "GlycoTree_OGlycans"
)


if [[ ! ("$(pwd)" =~ "/PyGly/scripts") ]]; then
    # Check for executing folder
    echo "Please execute from PyGly/scripts"
    exit 1
fi

if [ -z "$1" ]; then
    # Check for tag
    echo "Please provide the release tag number, eg V1.2.3"
    exit 1
fi

if [ -d ./GNOme ]; then
    echo "Please remove the GNOme directory"
    exit 1
fi

gh repo clone glygen-glycan-data/GNOme
# git clone git@github.com:glygen-glycan-data/GNOme.git
./gnome_prereq.sh
./gnome_compute.py writeowl \
  ./GNOme/data/gnome_subsumption_raw.txt \
  ./GNOme.owl \
  ./GNOme/data/mass_lookup_2decimal \
  ./GNOme/data/all_accession \
  -version $1 \
  -exact_sym1 ./GNOme/data/shortuckbcomp2glytoucan.txt \
  -exact_sym2 ./GNOme/data/shortcomp2glytoucan.txt \
  -byonic_sym ./GNOme/data/byonic2glytoucan.txt \
  -allExactSymOutput ./GNOme/data/exact_synonym.txt \
  -archive ./GNOme/data/glytoucan_archived.txt \
  -replace ./GNOme/data/glytoucan_replaced.txt \

./gnome_compute.py viewerdata ./GNOme.owl ./BrowserData.json
mv ./BrowserData.json ./GNOme/
./minify_json.py ./GNOme/BrowserData.json > ./GNOme/BrowserData.min.json
jq -r 'keys|@tsv'  ./GNOme/BrowserData.json | fmt -w 8 | sort -u > ./GNOme/valid-accessions.txt

for Restriction_set in "${restriction_set_names_standard[@]}"
do
  echo $Restriction_set
  ./gnome_compute.py writeresowl ./GNOme.owl ./GNOme/restrictions/GNOme_$Restriction_set.accessions.txt ./GNOme_$Restriction_set.owl
  ./gnome_compute.py viewerdata ./GNOme_$Restriction_set.owl ./GNOme/restrictions/$Restriction_set.BrowserData.json
  ./minify_json.py ./GNOme/restrictions/$Restriction_set.BrowserData.json > ./GNOme/restrictions/$Restriction_set.BrowserData.min.json
  jq -r 'keys|@tsv' ./GNOme/restrictions/$Restriction_set.BrowserData.json | fmt -w 8 | sort -u > ./GNOme/restrictions/$Restriction_set.valid-accessions.txt
done

for Restriction_set in "${restriction_set_names_ancestor[@]}"
do
  echo $Restriction_set
  ./gnome_compute.py writeresowl_with_ancestor_structures ./GNOme.owl ./GNOme/restrictions/GNOme_$Restriction_set.accessions.txt ./GNOme_$Restriction_set.owl
  ./gnome_compute.py viewerdata ./GNOme_$Restriction_set.owl ./GNOme/restrictions/$Restriction_set.BrowserData.json
  jq -r 'keys|@tsv' ./GNOme/restrictions/$Restriction_set.BrowserData.json | fmt -w 8 | sort -u > ./GNOme/restrictions/$Restriction_set.valid-accessions.txt
done

./gnome_compute.py UpdateTheme ./GNOme/restrictions ./GNOme/JS/theme/
for thjson in `ls ./GNOme/JS/theme/*.json | fgrep -v .min.json`; do
  echo "$thjson"
  thjson1=`basename "$thjson" .json`
  ./minify_json.py $thjson > $thjson1.min.json
done

./GNOme/convert.sh

cd ./GNOme
echo "Ready to commit & push"

git checkout -b "Branch_$1"
git add -A

git commit -m "Version $1"
git push origin "Branch_$1"

cd ..

rm -rf GNOme
