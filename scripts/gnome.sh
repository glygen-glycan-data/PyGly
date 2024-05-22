#!/bin/sh

restriction_set_names_standard=(
  "BCSDB"
  "GlyGen"
  "GlyCosmos"
  "PubChemCID"
  "GlyGen_NGlycans"
  "GlyGen_OGlycans"
  "NGlycans"
)

restriction_set_names_ancestor=(
  "GlycoTree_NGlycans"
  "GlycoTree_OGlycans"
)


help() {
  if [ "$1" != "" ]; then
    echo "" 1>&2
    echo "$1" 1>&2
  fi
  cat <<EOF | fmt -w80 -s 1>&2

./gnome.sh [ options ] [ <version> ]

Options:
  -r <FILE> GNOme raw data file
  -g <FILE> GlyGen accessions
  -t <TAG> Version tag
  -v Verbose
  -h Help

EOF
  exit 1;
}

GLYGENACC=""
GNOMERAW=""
VERBOSE=0
TAG="-"
while getopts ":r:g:tvh" o ; do
        case $o in
                r ) GNOMERAW="$OPTARG";;
                g ) GLYGENACC="$OPTARG";;
                t ) TAG="$OPTARG";;
                v ) VERBOSE=1;;
                h ) help "";;
                * ) help "Invalid option: -$OPTARG"
        esac
done

shift $(($OPTIND - 1))

if [ "$1" != "" ]; then
    TAG="$1"
fi

if [ "$VERBOSE" == 1 ]; then set -x; fi

if [[ ! ("$(pwd)" =~ "/PyGly/scripts") ]]; then
    # Check for executing folder
    echo "Please execute from PyGly/scripts"
    exit 1
fi

if [ -d ./GNOme -a "$TAG" != "-" ]; then
    echo "Please remove the GNOme directory"
    exit 1
fi

if [ ! -d ./GNOme ]; then
    gh repo clone glygen-glycan-data/GNOme
    ./gnome_prereq.sh
fi
if [ "$GNOMERAW" != "" ]; then
  cp "$GNOMERAW" ./GNOme/data/gnome_subsumption_raw.txt
fi
if [ "$GLYGENACC" != "" ]; then
  fgrep -v GlyTouCan "$GLYGENACC" > ./GNOme/restrictions/GNOme_GlyGen.accessions.txt
fi
./gnome_compute.py writeowl \
  ./GNOme/data/gnome_subsumption_raw.txt \
  ./GNOme.owl \
  ./GNOme/data/mass_lookup_2decimal \
  ./GNOme/data/all_accession \
  -version "$TAG" \
  -exact_sym1 ./GNOme/data/shortuckbcomp2glytoucan.txt \
  -exact_sym2 ./GNOme/data/shortcomp2glytoucan.txt \
  -byonic_sym ./GNOme/data/byonic2glytoucan.txt \
  -allExactSymOutput ./GNOme/data/exact_synonym.txt \
  -archive ./GNOme/data/glytoucan_archived.txt \
  -replace ./GNOme/data/glytoucan_replaced.txt \

./gnome_compute.py viewerdata ./GNOme.owl ./GNOme/BrowserData.json
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
  ./minify_json.py ./GNOme/restrictions/$Restriction_set.BrowserData.json > ./GNOme/restrictions/$Restriction_set.BrowserData.min.json
  jq -r 'keys|@tsv' ./GNOme/restrictions/$Restriction_set.BrowserData.json | fmt -w 8 | sort -u > ./GNOme/restrictions/$Restriction_set.valid-accessions.txt
done

./gnome_compute.py UpdateTheme ./GNOme/restrictions ./GNOme/JS/theme/
for thjson in `ls ./GNOme/JS/theme/*.json | fgrep -v .min.json`; do
  echo "$thjson"
  thjson1=`basename "$thjson" .json`
  ./minify_json.py $thjson > ./GNOme/JS/theme/$thjson1.min.json
done

./GNOme/convert.sh

if [ "$TAG" = "-" ]; then
    exit 0;
fi

cd ./GNOme
echo "Ready to commit & push"

git checkout -b "Branch_$TAG"
git add -A

git commit -m "Version $TAG"
git push origin "Branch_$TAG"

cd ..

rm -rf GNOme
