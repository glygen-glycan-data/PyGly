#!/bin/sh
if [ -s basecomplist.txt ]; then
   awk '$3 == "SHORTCOMP" {print $1"\t"$2}' basecomplist.txt | sort > shortcomp2glytoucan.txt
   awk '$3 == "SHORTUCKB" {print $1"\t"$2}' basecomplist.txt | sort > shortuckbcomp2glytoucan.txt
   awk '$3 == "UCKBCOMP" {print $1"\t"$2}' basecomplist.txt | sort > uckbcomp2glytoucan.txt
   awk '$3 == "BYONIC" {print $1"\t"$2}' basecomplist.txt | sort > byonic2glytoucan.txt
fi
