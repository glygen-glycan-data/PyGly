#!/bin/sh
set -x

PYTHON=python2
PYGLY=../../../pygly
DATA=../data

$PYTHON $PYGLY/GlycanResource/main.py GlyGen allglycans | sort > $DATA/glygen_accessions.txt

touch $DATA/glygen_exclude_accessions.txt
for accfile in $DATA/glygen_accessions.txt; do
    fgrep -v -f $DATA/glygen_exclude_accessions.txt ${accfile} > ${accfile}.new; mv -f ${accfile}.new $accfile
done

# Special/Other Requests
rm -f $DATA/glygen_req_accessions.txt
touch $DATA/glygen_req_accessions.txt

touch $DATA/glygen_manual_accessions.txt
cat $DATA/glygen_manual_accessions.txt >> $DATA/glygen_req_accessions.txt

# wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/current/glytoucan_ac.csv?token=ABCMDN6KPUJZYE6LIH62BR26RYZZ2' | awk -F, 'NR > 1 {print $1}' >> $DATA/glygen_req_accessions.txt

# All GlyTouCan motifs...
# echo "#GlycoMotif allmotifs GGM" >> $DATA/glygen_req_accessions.txt
$PYTHON $PYGLY/GlycanResource/main.py GlycoMotif allmotifs GGM | awk '{print $2}' >> $DATA/glygen_req_accessions.txt

for taxid in 9606 10090 10116 10114 111108 11116 11103 63746 694009 2697049 7227 4932 9823 44689; do
  # echo "#GlyTouCanNoCache bytaxa $taxid" >> $DATA/glygen_req_accessions.txt
  $PYTHON $PYGLY/GlycanResource/main.py GlyTouCanNoCache bytaxa $taxid >> $DATA/glygen_req_accessions.txt
  $PYTHON $PYGLY/GlycanResource/main.py GlyCosmosNoCache bytaxa $taxid >> $DATA/glygen_req_accessions.txt
done

# echo "#GlyGenSourceFile allgtc" >> $DATA/glygen_req_accessions.txt
$PYTHON $PYGLY/GlycanResource/main.py GlyGenSourceFile allsourcegtc | awk '{print $2}' >> $DATA/glygen_req_accessions.txt

cat $DATA/glygen_req_accessions.txt | \
  fgrep -v -f $DATA/glygen_accessions.txt | \
    sort -u > $DATA/glygen_new_accessions.txt

$PYTHON get_extra_accessions.py $DATA/gnome_subsumption_raw.txt $DATA/glygen_accessions.txt $DATA/glygen_new_accessions.txt | sort > $DATA/extra_accessions.txt

( python2 ./getbasecomplist.py \* > $DATA/basecomplist.txt ) > $DATA/basecomplist.log 
( cd $DATA; ./splitbasecomp.sh )
./uckbdata.sh
$PYTHON $PYGLY/GlycanResource/main.py GlyCosmosNoCache replaced | bash -c '(read -r; echo "$REPLY"; sort)' > $DATA/glytoucan_replaced.txt
$PYTHON $PYGLY/GlycanResource/main.py GlyCosmosNoCache archived | bash -c '(read -r; echo "$REPLY"; sort)' > $DATA/glytoucan_archived.txt
$PYTHON $PYGLY/GlycanResource/main.py GlyTouCanNoCache allaccessions | sort > $DATA/glytoucan_allacc.txt
$PYTHON $PYGLY/GlycanResource/main.py GlyCosmosNoCache allaccessions | sort > $DATA/glycosmos_allacc.txt

for accfile in $DATA/glygen_accessions.txt $DATA/glygen_new_accessions.txt $DATA/extra_accessions.txt; do
  fgrep -v -f $DATA/glytoucan_archived.txt $accfile > ${accfile}.new; mv -f ${accfile}.new $accfile
done
