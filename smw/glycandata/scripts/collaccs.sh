#!/bin/sh
set -x

PYTHON=apython
PYGLY=../../../pygly
DATA=../data

# $PYTHON $PYGLY/GlyGen.py allacc | sort -u > $DATA/glygen_accessions.txt
# python27 glygenacc.py | sort -u > $DATA/glygen_accessions.txt
$PYTHON $PYGLY/GlycanResource.py GlyGen allglycans > $DATA/glygen_accessions.txt

# Special/Other Requests
rm -f $DATA/glygen_req_accessions.txt
touch $DATA/glygen_req_accessions.txt

# wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/GlyGen_1.5/Other/GlyTouCan_accessions.csv?token=ABCMDN5KR6POBV3VV2MOPN25WHPNC' | awk -F, 'NR > 1 {print $1}' > $DATA/glygen_req_accessions.txt

wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/GlyGen_post_1.5/gtc_glyconnect_glycoprotein.csv?token=ABCMDNYY7VUUN4LBH673PKK6JLFDY' | awk -F, 'NR > 1 {print $1}' >> $DATA/glygen_req_accessions.txt

wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/GlyGen_post_1.5/gtc_hcv_glycoprotein.csv?token=ABCMDN2E6FRKOJC3ISME4IC6JLFKS' | awk -F, 'NR > 1 {print $1}' >> $DATA/glygen_req_accessions.txt

wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/GlyGen_post_1.5/gtc_synthesized_glycan_updated.csv?token=ABCMDNYVNXNKJPKSPR4EUZ26JLFMA' | awk -F, 'NR > 1 {print $1}' >> $DATA/glygen_req_accessions.txt

# wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/GlyGen_1.5/UniCarbKB/human_2019_08_23_02_49_01.csv?token=ABCMDN2FNFHTHHQRPVFIPHS5WHP2E' | awk -F, 'NR > 1 {print $4}' | grep -v '^$' >> $DATA/glygen_req_accessions.txt
# wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/GlyGen_1.5/UniCarbKB/mouse_2019_08_23_02_49_14.csv?token=ABCMDN7Y2WSBTCR4G7NP5G25WHQBO' | awk -F, 'NR > 1 {print $4}' | grep -v '^$' >> $DATA/glygen_req_accessions.txt
# wget -q -O - 'https://raw.githubusercontent.com/GW-HIVE/data_share/master/GlyGen_1.5/UniCarbKB/rat_2019_08_23_02_49_25.csv?token=ABCMDN4Q4PM2P7RDGTABJNS5WHQIG' | awk -F, 'NR > 1 {print $4}' | grep -v '^$' >> $DATA/glygen_req_accessions.txt

# Will's aligned GlycO alignments
awk -F, 'NR > 1 {print $1}' $DATA/canonres2gtc.csv | grep -v '^$' >> $DATA/glygen_req_accessions.txt

# All GlyTouCan motifs...
$PYTHON $PYGLY/GlycanResource.py GlyTouCan allmotifs | awk '{print $1}' >> $DATA/glygen_req_accessions.txt

for taxid in 9606 10090 10116 10114 11103; do
  $PYTHON $PYGLY/GlycanResource.py GlyTouCan bytaxa $taxid > $DATA/gtc.${taxid}.txt
  $PYTHON $PYGLY/GlycanResource.py UniCarbKB gtcbytaxa $taxid > $DATA/uc.${taxid}.txt
done

( ls $DATA/gtc.*.txt $DATA/uc.*.txt | egrep '(gtc|uc)\.[0-9][0-9]*\.txt$'; \
  ls "$DATA/glygen_req_accessions.txt" "$DATA/glygen_hist_accessions.txt" ) | \
  xargs -n 10 cat | \
    fgrep -v -f $DATA/glygen_accessions.txt | \
      sort -u > $DATA/glygen_new_accessions.txt

$PYTHON get_extra_accessions.py $DATA/gnome_subsumption_raw.txt $DATA/glygen_accessions.txt $DATA/glygen_new_accessions.txt > $DATA/extra_accessions.txt
