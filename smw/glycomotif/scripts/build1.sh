#!/bin/sh

SMWENV="DEV"
export SMWENV

set -euxo pipefail 
# ./loadallmotif.py
# ./addsameas.py
# ./widgetdata.py ../data/topologydev.json ../data/nonredonlydev.json ../data/redonlydev.json
# ( cd ../data; git add *dev.json; git commit -m "new widgetdata"; git push )
# ./loadwidgetbool.py
# ./addenzyme.py
# ./loadimpc.py ../data/phenotypeHitsPerGene.csv.gz
# ./loadHPOphenotype.py ../data/genes_to_phenotype.txt
# ./loadHPOdisease.py ../data/phenotype.hpoa ../data/genes_to_disease.txt
# ./loadHPAtissue.py ../data/proteinatlas.tsv
# ./loadHPAct.py ../data/proteinatlas.tsv
./enzymealignments.py | gzip -c > ../data/enzymealignments-dev.tsv.gz
./addenzymealign.py ../data/enzymealignments-dev.tsv.gz
./computealignments2.py -v --workers yion:8,bion:8,proton:8,maldi:8 -o "../data/computealignments-dev.tsv" > ../data/computealignments-dev.log 2>&1
./computeclassification.py < ../data/computealignments-dev.tsv > ../data/computeclass-dev.tsv
./alignments_tsv2rdf.py ../data/computealignments-dev.tsv ../data/computeclass-dev.tsv | gzip -9 -c > ../data/glycomotifdev_motifalign.rdf.gz
cp ../data/glycomotifdev_motifalign.rdf.gz /data/projects/smw/docker/glycomotifdev;
ssh trypsin "cd /data/projects/smw/docker/glycomotifdev; ../bin/dumprdf.sh";
DIR="$PWD"
cd /data/projects/smw/docker/glycomotifdev; \
    TRIP1=`ls glycomotifdev.*.rdf.gz | tail -n 1`; \
    TRIP2=glycomotifdev_motifalign.rdf.gz; \
    ssh -n trypsin "cd /data/projects/smw/docker/glycomotifdev; ../bin/loadts.sh $TRIP1 $TRIP2"
cd "$DIR"
./refresh.py &  
./refresh.py - & 
wait

