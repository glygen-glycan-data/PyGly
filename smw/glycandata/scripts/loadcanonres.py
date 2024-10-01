#!/bin/env python3.12

import sys, time
from collections import defaultdict
import csv, copy

from getwiki import GlycanData, Glycan
w = GlycanData()

gtc2enz = defaultdict(list)
glycotreeacc = set()

# glytoucan_ac,residue_name,residue_id,uniprot,gene_name,gene_id,parent_residue_id,type,species
for i,r in enumerate(csv.DictReader(open(sys.argv[1]))):
    gtcacc = r['glytoucan_ac']
    glycotreeacc.add(gtcacc)
    species = r['species']
    gtc2enz[(gtcacc,species)].append(copy.copy(r))

taxids={'mouse': '10090', 'human': '9606'}

for g in w.iterglycan():
    start = time.time()
    gtc = g.get('accession')

    g.delete_annotations(type='Enzyme')
    enzymes = set()
    for species in taxids:
        genes = set()
        for enz in gtc2enz[(gtc,species)]:
            genes.add(enz['gene_name'])
            enzymes.add(":".join([enz['uniprot'],enz['gene_name'],enz['gene_id'],taxids[enz['species']]]))
        g.set_annotation(value=list(genes),property='%s Enzyme'%(species.title()),type='Enzyme',source='GlycoTree')
    g.set_annotation(value=list(enzymes),property='EnzymeDetails',type='Enzyme',source='GlycoTree')
    if gtc in glycotreeacc:
        g.set_annotation(value=gtc,property="GlycoTree",type="CrossReference",source="GlycoTree")
    else:
        g.delete_annotations(property="GlycoTree",type="CrossReference",source="GlycoTree")

    if w.put(g):
        print("%s updated in %.2f sec"%(gtc,time.time()-start,),file=sys.stderr)
    else:
        print("%s checked in %.2f sec"%(gtc,time.time()-start,),file=sys.stderr),file=sys.stderr)
