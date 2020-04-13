#!/bin/sh
set -x

PYTHON=apython
PYGLY=../../../pygly
DATA=../data

# $PYTHON $PYGLY/GlyGen.py allacc | sort -u > $DATA/glygen_accessions.txt
# python27 glygenacc.py | sort -u > $DATA/glygen_accessions.txt
$PYTHON $PYGLY/GlycanResource/main.py GlyGen allglycans > $DATA/glygen_accessions.txt

# Special/Other Requests
rm -f $DATA/glygen_req_accessions.txt
touch $DATA/glygen_req_accessions.txt

touch $DATA/glygen_manual_accessions.txt
cat $DATA/glygen_manual_accessions.txt >> $DATA/glygen_req_accessions.txt

wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/current/glytoucan_ac.csv?token=ABCMDN6KPUJZYE6LIH62BR26RYZZ2' | awk -F, 'NR > 1 {print $1}' >> $DATA/glygen_req_accessions.txt

# Will's aligned GlycO alignments
awk -F, 'NR > 1 {print $1}' $DATA/canonres2gtc.csv | grep -v '^$' >> $DATA/glygen_req_accessions.txt

# All GlyTouCan motifs...
$PYTHON $PYGLY/GlycanResource/main.py GlyTouCan allmotifs | awk '{print $1}' >> $DATA/glygen_req_accessions.txt

for taxid in 9606 10090 10116 10114 11103; do
  $PYTHON $PYGLY/GlycanResource/main.py GlyTouCan bytaxa $taxid > $DATA/gtc.${taxid}.txt
  $PYTHON $PYGLY/GlycanResource/main.py UniCarbKB gtcbytaxa $taxid > $DATA/uc.${taxid}.txt
done

( ls $DATA/gtc.*.txt $DATA/uc.*.txt | egrep '(gtc|uc)\.[0-9][0-9]*\.txt$'; \
  ls "$DATA/glygen_req_accessions.txt" "$DATA/glygen_hist_accessions.txt" ) | \
  xargs -n 10 cat | \
    fgrep -v -f $DATA/glygen_accessions.txt | \
      sort -u > $DATA/glygen_new_accessions.txt

$PYTHON get_extra_accessions.py $DATA/gnome_subsumption_raw.txt $DATA/glygen_accessions.txt $DATA/glygen_new_accessions.txt > $DATA/extra_accessions.txt
