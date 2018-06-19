#!/bin/sh
THEDIR=`dirname $0`
THEDIR=`readlink -f $THEDIR`
if [ ! -f $THEDIR/taxonomy/ncbitaxonomytree.sqlite3 ]; then
  if [ ! -f $THEDIR/../data/taxdmp.zip ]; then
    wget -nv -nH --cut-dirs=2 --passive-ftp \
         -O $THEDIR/../data/taxdmp.zip \
         ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip
  fi
  $THEDIR/maketax.py $THEDIR/../data
  rm -f $THEDIR/../data/taxdmp.zip
fi
