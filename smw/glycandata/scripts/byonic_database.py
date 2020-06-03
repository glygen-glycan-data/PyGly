#!/bin/env python27

import re, sys
from getwiki import GlycanData
from collections import defaultdict

import findpygly
from pygly.ElementMass import MonoisotopicElementMass
from pygly.CompositionTable import Composition
mt = MonoisotopicElementMass()
h2o = Composition(H=2,O=1).mass(mt)
 
w = GlycanData()

species = 'Human'

output = []
for acc in w.iterglycanid():
    # print acc
    g = w.get(acc)
    
    # Have a Byonic name...
    try:
        byonic = g.get_annotation_values(property="Byonic",type="Name",source="EdwardsLab")[0]
    except LookupError:
	continue

    # Prefer dHex version...
    if 'Fuc' in byonic:
	continue

    # is of the species by direct annotation or subsumption
    try:
        cat = g.get_annotation_value(property='%s Category'%(species,),type="Species",source="EdwardsLab")
    except LookupError:
	continue

    if cat not in ("Direct","Subsumption"):
	continue

    # is itself, or has a descendant annotated as N-linked
    descs = set([acc])
    try:
        descs.update(set(g.get_annotation_values(property="Descendant",type="Subsumption",source="GNOme")))
    except LookupError:
	pass
    nlinked = False
    for desc in descs:
	dg = w.get(desc)
	try:
	    type = dg.get_annotation_value(property="GlycanType",type="Classification",source="EdwardsLab")
	    if type == 'N-linked':
		nlinked = True
	        break
	except LookupError:
	    continue

    if not nlinked:
	continue
	
    # has molecular weight
    try:
        mw = float(g.get_annotation_value(property='UnderivitizedMW',type="MolWt",source="EdwardsLab"))
    except LookupError:
	continue

    mw -= h2o # account for attachment to the peptide...

    output.append(dict(accession=acc,mw=mw,byonic=byonic))

for d in sorted(output,key=lambda d: d.get('mw',1e+20)):
    print "%(byonic)s %% %(mw).6f, %(accession)s"%d
