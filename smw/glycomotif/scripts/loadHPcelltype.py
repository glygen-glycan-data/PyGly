#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

e_id = defaultdict(list)
c_type = defaultdict(list)
celltype= defaultdict(list)

header = ["ensembl_id", "gene_symbol", "tissue_type", "expression"]
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t", fieldnames=header):
	 c_type[row['gene_symbol'].upper()] = [row["tissue_type"], row["expression"]]
         e_id[row['gene_symbol'].upper()] = row["ensembl_id"]



for k, v in c_type.items():
      
      if "detected" in v[1]:
         detected = v[1].partition(",")[0]
         detected = detected.replace("%", "")
         high = v[1].partition(",")[1]
         high = detected.replace("%", "")

         if float(detected.split("=")[1]) >= 80:
            celltype[k] = ["Ubiquitous, highly expressed"]
         elif float(high.split("=")[1]) >= 80:
            celltype[k] = ["Ubiquitous"]

      if v[0] == "Group enriched":
         for exp in v[1].split(","): 
            celltype[k].append("Enriched in " + exp)
      elif v[0] == "Moderately enriched":
         celltype[k] = ["Enriched in " + v[1]]
      elif v[0] == "Highly enriched":
         celltype[k] = ["Highly enriched in " + v[1]]



#print(tissue)


for e in w.iterenzyme():
    gn = e.get('genename')
    species = e.get('species')
    if species == "Mouse":
       continue
    if len(celltype.get(gn,[])) > 0:
       e.set("celltype_source_key",(e_id[gn]))
       e.set("celltype_source","GTEX")
       sorted_list = sorted(celltype[gn]) if isinstance(celltype[gn], list) else celltype[gn]
       e.set("celltype", sorted_list)
       e.delete("cell_source_key")
       e.delete("cell_source")
       e.delete("cell")
       print(e)
    else:	
       e.delete("celltype_source_key")
       e.delete("celltype")


    if w.put(e):
       print (gn)

