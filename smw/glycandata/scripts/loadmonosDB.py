#!/bin/env python27

import sys

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

def monosDB():
    d = {}
    f = open('../data/GlyGen_glycans2monodbID.tsv','r').read().strip().split("\n")
    
    for line in f:
	k,v = line.split("\t")
	if v == "Some monosaccharide(s) has no ID in monosaccharide DB" or v == "MonosaccharideDB_ID":
	    continue
	if k in d:
	    d[k].append(v)
	else:
	    d[k] = [v]    
    return d

monosdb = monosDB()

for m in w.iterglycan():
    acc = m.get('accession')
    if acc in monosdb:
        m.delete_annotations(source="EdwardsLab",type="CrossReference")
        m.set_annotation(value=monosdb[acc], property="MonosaccharideDB",source="EdwardsLab",type="CrossReference")
    if w.put(m):
        print acc


