#!/bin/sh
set -x

PYTHON=apython
PYGLY=../../../pygly
DATA=../data

$PYTHON $PYGLY/GlycanResource/main.py GlyGen allglycans > $DATA/glygen_accessions.txt
# $PYTHON $PYGLY/GlycanResource/main.py GlyGen alltimeglycans > $DATA/glygen_hist_accessions.txt

# Special/Other Requests
rm -f $DATA/glygen_req_accessions.txt
touch $DATA/glygen_req_accessions.txt

touch $DATA/glygen_manual_accessions.txt
cat $DATA/glygen_manual_accessions.txt >> $DATA/glygen_req_accessions.txt

# wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/current/glytoucan_ac.csv?token=ABCMDN6KPUJZYE6LIH62BR26RYZZ2' | awk -F, 'NR > 1 {print $1}' >> $DATA/glygen_req_accessions.txt

# All GlyTouCan motifs...
$PYTHON $PYGLY/GlycanResource/main.py GlycoMotif allmotifs GGM | awk '{print $2}' >> $DATA/glygen_req_accessions.txt

for taxid in 9606 10090 10116 10114 11103 11108 694009 2697049; do
  $PYTHON $PYGLY/GlycanResource/main.py GlyTouCanNoCache bytaxa $taxid >> $DATA/glygen_req_accessions.txt
done

$PYTHON $PYGLY/GlycanResource/main.py GlyGenSourceFile allgtc >> $DATA/glygen_req_accessions.txt

cat $DATA/glygen_req_accessions.txt | \
  fgrep -v -f $DATA/glygen_accessions.txt | \
    sort -u > $DATA/glygen_new_accessions.txt

$PYTHON get_extra_accessions.py $DATA/gnome_subsumption_raw.txt $DATA/glygen_accessions.txt $DATA/glygen_new_accessions.txt > $DATA/extra_accessions.txt

( apython ./getbasecomplist.py \* > $DATA/basecomplist.txt ) > $DATA/basecomplist.log 
( cd $DATA; ./splitbasecomp.sh )
./uckbdata.sh
$PYTHON $PYGLY/GlycanResource/main.py GlyCosmosNoCache replaced > $DATA/glytoucan_replaced.txt
$PYTHON $PYGLY/GlycanResource/main.py GlyCosmosNoCache archived > $DATA/glytoucan_archived.txt
 
