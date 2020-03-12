#!/bin/sh

restriction_set_names=(
  "BCSDB"
  "GlyGen"
  "GlycanData"
)


if [[ ! ("$(pwd)" =~ "/PyGly/scripts") ]]; then
  # Check for executing folder
  exit
fi

if [ -z "$1" ]
  then
    # Check for tag
    echo "Please provide the release tag number, eg V1.2.3"
    exit
fi

git clone git@github.com:glygen-glycan-data/GNOme.git
# python27 ../pygly/GNOme.py writeowl ./GNOme/data/gnome_subsumption_raw.txt ./GNOme.owl ./GNOme/data/mass_lookup_2decimal $1
python27 ../pygly/GNOme.py writeowl ./GNOme/data/gnome_subsumption_raw.txt ./GNOme.owl ./GNOme/data/mass_lookup_2decimal -version $1 -exact_sym1 ./GNOme/data/shortuckbcomp2glytoucan.txt -byonic_sym ./GNOme/data/byonic2glytoucan.txt
python27 ../pygly/GNOme.py viewerdata ./GNOme.owl ./GNOme.browser.json ./GNOme.browser.composition.json
# python27 ../pygly/GNOme.py viewerdata2 ./GNOme.owl ./GNOme.browser2.json ./GNOme.browser.composition.json


for Restriction_set in "${restriction_set_names[@]}"
do
  lowersetname=$(echo "$Restriction_set" | awk '{print tolower($0)}')
  echo Restriction_set
  # python27 ../pygly/GNOme.py UpdateAcc $Restriction_set ./GNOme/restrictions/GNOme_$Restriction_set.accessions.txt ./GNOme/JS/"$lowersetname"_accession.json
  python27 ../pygly/GNOme.py writeresowl ./GNOme.owl ./GNOme/restrictions/GNOme_$Restriction_set.accessions.txt ./GNOme_$Restriction_set.owl
  python27 ../pygly/GNOme.py viewerdata ./GNOme_$Restriction_set.owl ./GNOme_$Restriction_set.browser.json ./GNOme_$Restriction_set.browser.composition.json
done

python27 ../pygly/GNOme.py UpdateTheme ./GNOme/restrictions ./GNOme/JS/theme/

cp ./GNOme/convert.sh ./
./convert.sh
mv ./GNOme.browser.* ./GNOme/

for Restriction_set in "${restriction_set_names[@]}"
do
	for file_ext in "owl" "obo" "json" "browser.json" "browser.composition.json"
  do
    	mv ./GNOme_$Restriction_set.$file_ext ./GNOme/restrictions/
  done
done


cd ./GNOme
echo "Ready to commit & push"

git checkout -b "Branch_$1"
git add -A

git commit -m "Version $1"
git push origin "Branch_$1"

rm -rf GNOme
rm convert.sh