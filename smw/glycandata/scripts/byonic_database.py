#!/bin/env python3.12

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

    # Is N-linked by subsumption or glycomotif
    nlinked = False
    for cls in g.get_all_annotations_values(property='GlycanType',type='Classification'):
        if cls == 'N-linked':
            nlinked = True
            break

    if not nlinked:
        continue

    # N-linked truncated cannot be the only reason it is annotated as N-linked
    nlinked_sc = 0
    nlinked_scnottrunc = 0
    for subcls in g.get_all_annotations_values(property='GlycanSubtype',type='Classification'):
        if subcls.startswith('N-linked'):
            nlinked_sc += 1
            if subcls != "N-linked truncated" and subcls != "N-linked monoantennary":
              nlinked_scnottrunc += 1

    if nlinked_sc > 0 and nlinked_scnottrunc == 0:
        continue

    # has molecular weight
    try:
        mw = float(g.get_annotation_value(property='UnderivitizedMW',type="MolWt",source="EdwardsLab"))
    except LookupError:
        continue

    mw -= h2o # account for attachment to the peptide...

    output.append(dict(accession=acc,mw=mw,byonic=byonic))

for d in sorted(output,key=lambda d: d.get('mw',1e+20)):
    print("%(byonic)s %% %(mw).6f, %(accession)s"%d)
