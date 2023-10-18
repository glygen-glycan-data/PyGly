#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

dis_id = defaultdict(list)
gdisease = defaultdict(list)
disease = defaultdict(list)
geneid = defaultdict(list)


# "ncbi_gene_id","gene_symbol","association_type",  "disease_id", "source"
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t"):
	 dis_id[row['gene_symbol'].upper()].append(row["disease_id"])	
         geneid[row['gene_symbol'].upper()] = row["ncbi_gene_id"].partition("e:")[2]


#database_id	disease_name	qualifier	hpo_id	reference	evidence	onset	frequency	sex	modifier	aspect	biocuration
for i, row in enumerate(csv.reader((open(sys.argv[2])), delimiter = "\t")):
	if i == 4:
	   header = row



for row in csv.DictReader((open(sys.argv[2])), delimiter = "\t", fieldnames=header):
	disease[row["database_id"]] = row["disease_name"]


print(dis_id)
for gene, v in dis_id.items():
   for ids in v:
       if ids in disease:
          #print(disease[ids])
          gdisease[gene].append(disease[ids])
       else:
         gdisease[gene] = []

#print(gdisease)

for e in w.iterenzyme():
    gn = e.get('genename')
    species = e.get('species')
    if species == "Mouse":
       continue
    if len(gdisease.get(gn,[])) > 0:
       e.set("disease_source_key",(geneid[gn]))
       e.set("disease_source","HPO")
       e.set("disease",(sorted(gdisease[gn])))
       print(e)
    else:	
       e.delete("disease")
       e.delete("disease_source_key")


    if w.put(e):
       print (gn)

