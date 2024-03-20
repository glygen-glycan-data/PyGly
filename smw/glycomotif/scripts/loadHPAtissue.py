#!/bin/env python2

from getwiki import GlycoMotifWiki
from collections import defaultdict
import sys, csv 

w = GlycoMotifWiki()
tissue = dict()
e_id = dict()

seen_tissues = set()

tissue_map = dict(filter(lambda t: len(t) == 2,map(lambda l: map(str.strip,l.split(':',1)),"""
stomach 1: stomach
skin 1: skin
endometrium 1: endometrium
""".splitlines())))

spec_map = dict(filter(lambda t: len(t) == 2,map(lambda l: map(str.strip,l.split(':',1)),"""
Tissue enriched: Enriched
Tissue enhanced: Enhanced
Group enriched: Enriched
""".splitlines())))

#
# loadHPAtissue.py ../data/proteinatlas.tsv
#
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t"): #fieldnames=header):
   tissue[row['Gene']] = []
   spec = row["RNA tissue specificity"].strip()
   spec = spec_map.get(spec,spec)
   if spec in ("Enriched","Enhanced"):
       for tntpm in row['RNA tissue specific nTPM'].split(';'):
           thetissue = tntpm.split(':',1)[0].strip()
           thetissue = tissue_map.get(thetissue,thetissue)
           if thetissue not in seen_tissues:
               seen_tissues.add(thetissue)
           tissue[row['Gene']].append(spec + ' in ' + thetissue)
   else:
       tissue[row['Gene']].append(spec)
   e_id[row['Gene']] = row["Ensembl"]




# print(tissue)


for e in w.iterenzyme():
    gn = e.get('genename')
    if len(tissue.get(gn,[])) > 0:
       e.set("tissue_source_key",e_id[gn])
       e.set("tissue_source","HPA") 
       e.set("tissue", sorted(tissue[gn]))
    else:	
       e.delete("tissue_source_key")
       e.delete("tissue_source")
       e.delete("tissue")
    if w.put(e):
       print (gn)
