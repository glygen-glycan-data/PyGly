#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

dis_id = defaultdict(list)
gdisease = defaultdict(set)
disease = defaultdict(list)
geneid = defaultdict(list)



#database_id	disease_name	qualifier	hpo_id	reference	evidence	onset	frequency	sex	modifier	aspect	biocuration
for i, row in enumerate(csv.reader((open(sys.argv[1])), delimiter = "\t")):
	if i == 4:
	   header = row


counter = 0
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t", fieldnames=header):
   if counter < 6:
      counter += 1
      continue

   else:
      disease[row["database_id"]] = row["disease_name"]


# "ncbi_gene_id","gene_symbol","association_type",  "disease_id", "source"
for row in csv.DictReader((open(sys.argv[2])), delimiter = "\t"):
   geneid[row['gene_symbol']] = row["ncbi_gene_id"].partition(":")[2]

   if row["disease_id"] in disease:
      gdisease[row["gene_symbol"]].add(disease[row["disease_id"]])
         

#print(gdisease)


for e in w.iterenzyme():
    gn = e.get('genename')
    species = e.get('species')
    if species != "Human":
       continue
    if len(gdisease.get(gn,[])) > 0:
       e.set("disease_source_key",(geneid[gn]))
       e.set("disease_source","HPO")
       e.set("disease",(sorted(gdisease[gn])))
       print(e)
    else:	
       e.delete("disease")
       e.delete("disease_source_key")
       e.delete("disease_source")


    if w.put(e):
       print (gn)

