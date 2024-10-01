#!/bin/sh
set -x

GLYRES=../../../scripts/glyres.py
DATA=../data

$GLYRES GlyGen allglycans | sort > $DATA/glygen_accessions.txt

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
$GLYRES GlycoMotifNoCache allmotifs GGM | awk '{print $2}' >> $DATA/glygen_req_accessions.txt

for taxid in 9606 10090 10116 10114 111108 11116 11103 3052230 63746 694009 2697049 7227 4932 9823 44689 9031 3702 9913; do
  $GLYRES GlyTouCanNoCache bytaxa $taxid >> $DATA/glygen_req_accessions.txt
  $GLYRES GlyCosmosNoCache bytaxa $taxid >> $DATA/glygen_req_accessions.txt
done

$GLYRES GlyGenSourceFile allsourcegtc | awk '{print $2}' >> $DATA/glygen_req_accessions.txt

cat $DATA/glygen_req_accessions.txt | \
  fgrep -v -f $DATA/glygen_accessions.txt | \
    sort -u > $DATA/glygen_new_accessions.txt

./get_extra_accessions.py $DATA/gnome_subsumption_raw.txt $DATA/glygen_accessions.txt $DATA/glygen_new_accessions.txt | sort > $DATA/extra_accessions.txt

( ./getbasecomplist.py \* > $DATA/basecomplist.txt ) > $DATA/basecomplist.log 
( cd $DATA; ./splitbasecomp.sh )
./uckbdata.sh
$GLYRES GlyCosmosNoCache replaced | bash -c '(read -r; echo "$REPLY"; sort)' > $DATA/glytoucan_replaced.txt
$GLYRES GlyCosmosNoCache archived | bash -c '(read -r; echo "$REPLY"; sort)' > $DATA/glytoucan_archived.txt
$GLYRES GlyTouCanNoCache allaccessions | sort > $DATA/glytoucan_allacc.txt
$GLYRES GlyCosmosNoCache allaccessions | sort > $DATA/glycosmos_allacc.txt

for accfile in $DATA/glygen_accessions.txt $DATA/glygen_new_accessions.txt $DATA/extra_accessions.txt; do
  fgrep -v -f $DATA/glytoucan_archived.txt $accfile > ${accfile}.new; mv -f ${accfile}.new $accfile
done
