#!/bin/env python2

from getwiki import GlycoMotifWiki
from collections import defaultdict
import sys, csv 

w = GlycoMotifWiki()
tissue = defaultdict(list)
e_id = defaultdict(list)


#header = ["ensembl_id", "gene_symbol", "tissue_type", "expression"]
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t"): #fieldnames=header):
   tissue[row['Gene']] = row["RNA tissue specificity"]
   e_id[row['Gene']] = row["Ensembl"]




print(tissue)


for e in w.iterenzyme():
    gn = e.get('genename')
    print(gn)
    if len(tissue.get(gn,[])) > 0:
       e.set("tissue_source_key",(tissue[gn]))
       e.set("tissue_source","HPA") 
       e.set("tissue", sorted(tissue[gn]))
       print(e)
    else:	
       e.delete("tissue_source_key")
       e.delete("tissue_source")
       e.delete("tissue")

"""
    if w.put(e):
       print (gn)
"""
