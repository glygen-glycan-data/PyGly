#!/bin/env python3.12

import sys, time
from collections import defaultdict

from getwiki import GlycanData
w = GlycanData()

import findpygly
from pygly.GlycanResource import PubChemDownload

pubchem = defaultdict(set)
chebi = defaultdict(set)

pcd = PubChemDownload(verbose=False)

for pcid,gtc in pcd.allgtc():
    pubchem[gtc].add(pcid)

for chebiid,gtc in pcd.allchebigtc():
    chebi[gtc].add(chebiid)

for m in w.iterglycan():
    start=time.time()
    acc = m.get('accession')
    m.delete_annotations(property="PubChem",type="CrossReference")
    m.delete_annotations(source="PubChem",property="ChEBI",type="CrossReference")
    if len(list(pubchem[acc])) > 0:
        m.set_annotation(value=list(pubchem[acc]), property="PubChem",source="PubChem",type="CrossReference")
    else:
        m.delete_annotations(property="PubChem",source="PubChem",type="CrossReference")
    if len(chebi[acc]) > 0:
        m.set_annotation(value=list(chebi[acc]), property="ChEBI",source="PubChem",type="CrossReference")
    else:
        m.delete_annotations(property="ChEBI",source="PubChem",type="CrossReference")
    if w.put(m):
        print("%s updated in %.2f sec"%(acc,time.time()-start,),file=sys.stderr)
    else:
        print("%s checked in %.2f sec"%(acc,time.time()-start,),file=sys.stderr)
