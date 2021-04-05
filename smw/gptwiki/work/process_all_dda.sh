#!/bin/sh
set -x
DDA_DIRS=" \
  HEK_cells_nonenriched \
  IgG_Glycopeptide_HCD_V1 \
  MS_24_U01 \
  MS_31_U01 \
  MS_35_Frac_plasma_run01 \
  MS_35_Frac_plasma_run02 \
  MS_47_UO1_HEK \
  MS_48_UO1_Erlic \
  MS_50_ACE2 \
  MS_50_SPIKE \
  MS_51_plasma \
"

for workdir in $DDA_DIRS; do 
  ( cd $workdir; ../process.sh *.mzML.gz *.msp > process.log 2>&1 ; ../export.sh DEV *.mzML.gz *.msp > export.log 2>&1 )
done
