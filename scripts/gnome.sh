#!/bin/sh

set -x

restriction_set_names=(
  "BCSDB"
  "GlyGen"
  "GlyCosmos"
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

git clone git@github.com:glygen-glycan-data/GNOme.git
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


for Restriction_set in "${restriction_set_names[@]}"
do
  lowersetname=$(echo "$Restriction_set" | awk '{print tolower($0)}')
  echo $Restriction_set
  # python ../pygly/GNOme.py UpdateAcc $Restriction_set ./GNOme/restrictions/GNOme_$Restriction_set.accessions.txt ./GNOme/JS/"$lowersetname"_accession.json
  ./gnome_compute.py writeresowl ./GNOme.owl ./GNOme/restrictions/GNOme_$Restriction_set.accessions.txt ./GNOme_$Restriction_set.owl
  ./gnome_compute.py viewerdata ./GNOme_$Restriction_set.owl ./$Restriction_set.BrowserData.json
done

./gnome_compute.py UpdateTheme ./GNOme/restrictions ./GNOme/JS/theme/

cp ./GNOme/convert.sh ./
./convert.sh



mv *.BrowserData.json ./GNOme/restrictions/


cd ./GNOme
echo "Ready to commit & push"

git checkout -b "Branch_$1"
git add -A

git commit -m "Version $1"
git push origin "Branch_$1"

rm -rf GNOme
rm convert.sh
