#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

e_id = defaultdict(list)
tis_type = defaultdict(list)
tissue = defaultdict(list)

header = ["ensembl_id", "gene_symbol", "tissue_type", "expression"]
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t", fieldnames=header):
	 tis_type[row['gene_symbol']] = [row["tissue_type"], row["expression"]]
         e_id[row['gene_symbol']] = row["ensembl_id"]



for k, v in tis_type.items():
      
      if "detected" in v[1]:
         part = v[1].partition(",")
         detected = part[0]
         detected = detected.replace("%", "")
         high = part[2]
         high = high.replace("%", "")
         if float(high.split("=")[1]) >= 80:
            tissue[k] = ["Ubiquitous, highly expressed"]
         elif float(detected.split("=")[1]) >= 80:
            tissue[k] = ["Ubiquitous"]

      if v[0] == "Group enriched":
         for exp in v[1].split(","): 
            tissue[k].append("Enriched in " + exp)
      elif v[0] == "Moderately enriched":
         tissue[k] = ["Enriched in " + v[1]]
      elif v[0] == "Highly enriched":
         tissue[k] = ["Highly enriched in " + v[1]]




#print(tissue)


for e in w.iterenzyme():
    gn = e.get('genename')
    species = e.get('species')
    if species != "Human":
       continue
    if len(tissue.get(gn,[])) > 0:
       e.set("tissue_source_key",(e_id[gn]))
       e.set("tissue_source","GTEx") 
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
