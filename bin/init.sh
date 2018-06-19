#!/bin/sh
set -x
THEDIR=`dirname $0`
THEDIR=`readlink -f $THEDIR`
$THEDIR/maketax.sh
cd $THEDIR/../data
$THEDIR/makedb.sh Named
$THEDIR/GlyDbIndex.sh build CFGArray-Mammalian-v5.1.cfg
$THEDIR/GlyDbIndex.sh glycoct CFGArray-Mammalian-v5.1.cfg
mkdir CFGArray-Mammalian-v5.1
(cd CFGArray-Mammalian-v5.1; unzip ../CFGArray-Mammalian-v5.1.gct)
$THEDIR/makedb.sh CFGArray-Mammalian-v5.1
rm -f CFGArray-Mammalian-v5.1.cfg.index
rm -rf CFGArray-Mammalian-v5.1
$THEDIR/getGlycomeDb.sh
$THEDIR/partition_glycomedb.sh glycomedb.xml.gz
for gdb in *.gdb; do
  $THEDIR/addimages.sh $gdb
done
$THEDIR/buildall.sh
$THEDIR/countall.sh
