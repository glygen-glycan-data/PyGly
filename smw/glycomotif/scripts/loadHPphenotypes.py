#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

geneid = {}
phenotypes = defaultdict(set)
disease = defaultdict(list)

# "ncbi_gene_id","gene_symbol","hpo_id", "hpo_name","frequency",  "disease_id", 
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t"):
    disease[row['gene_symbol'].upper()] = row['disease_id']
    geneid[row['gene_symbol'].upper()] = row['ncbi_gene_id']
    phenotypes[row['gene_symbol'].upper()].add(row["hpo_name"])



for k, v in phenotypes.items():
    v = sorted(v)
    phenotypes[k] = v





for e in w.iterenzyme():
    gn = e.get('genename')
    species = e.get('species')
    if species == 'Mouse':
       continue
        
    if len(disease.get(gn,[])) > 0: 
       e.set("phenotype",(phenotypes[gn]), )
       e.set("phenotype_source_key",(geneid[gn]))
       e.set("phenotype_source","HPO")
       print(e)
    else:	
       e.delete("phenotype")
       e.delete("phenotype_source_key")


    if w.put(e):
       print (gn)

