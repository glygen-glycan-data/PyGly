#!/bin/env python2

from getwiki import GlycoMotifWiki
from collections import defaultdict
import sys, csv 

w = GlycoMotifWiki()
celltype = defaultdict(set)
e_id = dict()

seen_celltypes = set()

celltype_map = dict(filter(lambda t: len(t) == 2,map(lambda l: map(str.strip,l.split(':',1)),"""
Alveolar cells type 1: Alveolar cells
Alveolar cells type 2: Alveolar cells
""".splitlines())))

spec_map = dict(filter(lambda t: len(t) == 2,map(lambda l: map(str.strip,l.split(':',1)),"""
Cell type enriched: Enriched
Cell type enhanced: Enhanced
Group enriched: Enriched
""".splitlines())))

#
# loadHPAct.py ../data/proteinatlas.tsv
#
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t"): #fieldnames=header):
   spec = row["RNA single cell type specificity"].strip()
   spec = spec_map.get(spec,spec)
   if spec in ("Enriched","Enhanced"):
       for tntpm in row['RNA single cell type specific nTPM'].split(';'):
           thecelltype = tntpm.split(':',1)[0].strip()
           thecelltype = celltype_map.get(thecelltype,thecelltype)
           if thecelltype not in seen_celltypes:
               seen_celltypes.add(thecelltype)
               # print(thecelltype)
           celltype[row['Gene']].add(spec + ' in ' + thecelltype)
   else:
       celltype[row['Gene']].add(spec)
   e_id[row['Gene']] = row["Ensembl"]

for e in w.iterenzyme():
    gn = e.get('genename')
    if len(celltype.get(gn,[])) > 0:
       e.set("celltype_source_key",e_id[gn])
       e.set("celltype_source","HPA") 
       e.set("celltype", sorted(celltype[gn]))
    else:	
       e.delete("celltype_source_key")
       e.delete("celltype_source")
       e.delete("celltype")
    if w.put(e):
       print (gn)
