#!/bin/sh
BINDIR=`dirname $0`
BINDIR=`readlink -f $BINDIR`
DATDIR=`readlink -f "$BINDIR/../data"`
rm -f $BINDIR/taxonomy/ncbitaxonomytree.sqlite3
rm -f $DATDIR/taxdmp.zip
rm -f $DATDIR/*.{gdb,gct,gdb.index,gct.index,xml.gz}
rm -rf dist build
find . -name "*~" -exec rm -f {} \;
find . -name "*.pyc" -exec rm -f {} \;
