#!/bin/sh
THEDIR=`dirname $0`
GDB="$1"
if [ -z "$GDB" ]; then
  echo "partition_glycomedb.sh glycomedb.xml.gz"
  exit 1
fi
PARTGDB="$THEDIR/partition_glycomedb.py --glycomedb $GDB"

$PARTGDB --taxid 9606 --out GlycomeDb-Human.gdb 
$PARTGDB --taxid 9606 --out GlycO-GlycomeDb-Human.gdb --resource glyco
$PARTGDB --taxid 9606 --out KEGG-GlycomeDb-Human.gdb --resource kegg

$PARTGDB --taxid 40674 --out GlycomeDb-Mammalian.gdb 
$PARTGDB --taxid 40674 --out GlycO-GlycomeDb-Mammalian.gdb --resource glyco
$PARTGDB --taxid 40674 --out KEGG-GlycomeDb-Mammalian.gdb --resource kegg

