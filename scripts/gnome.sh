#!/bin/sh

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
python ../pygly/GNOme.py writeowl ../smw/glycandata/data/gnome_subsumption_raw.txt ./GNOme.owl ./GNOme/id_mapping/mass_lookup_2decimal $1
python ../pygly/GNOme.py writeresowl ./GNOme.owl BCSDB ./GNOme_BCSDB.owl
python ../pygly/GNOme.py writeresowl ./GNOme.owl GlyGen ./GNOme_GlyGen.owl
exit
# TODO generation js data script here


cp ./GNOme/convert.sh ./
./convert.sh
mv ./GNOme.* ./GNOme/

for Restriction_set in "BCSDB" "GlyGen"
do
	for file_ext in "owl" "obo" "json" "browser.js"
  do
      # touch ./GNOme_$Restriction_set.$file_ext
    	mv ./GNOme_$Restriction_set.$file_ext ./GNOme/restrictions/
  done
done


cd ./GNOme
# echo $(pwd)

git checkout -b "Branch_$1"
git add -A

git commit -m "Version $1"
git push origin "Branch_$1"

rm -rf GNOme
rm convert.sh