#!/bin/env python27

import sys
from collections import defaultdict

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()


def monosDB(infile):
    d = defaultdict(set)
    f = open(infile,'r')
    
    for line in f:
	if "Some monosaccharide" in line or "MonosaccharideDB_ID" in line:
	    continue
	k,v = line.split()
	d[k].add(v)
    return d

monosdb = monosDB(sys.argv[1])

for m in w.iterglycan():
    acc = m.get('accession')
    m.delete_annotations(source="EdwardsLab",property="MonosaccharideDB",type="CrossReference")
    if acc in monosdb:
        m.set_annotation(value=list(monosdb[acc]), property="MonosaccharideDB",source="EdwardsLab",type="CrossReference")
    if w.put(m):
        print acc


