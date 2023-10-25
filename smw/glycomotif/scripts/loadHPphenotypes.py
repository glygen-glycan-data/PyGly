#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

geneid = {}
phenotypes = defaultdict(set)

# "ncbi_gene_id","gene_symbol","hpo_id", "hpo_name","frequency",  "disease_id", 
for row in csv.DictReader((open(sys.argv[1])), delimiter = "\t"):
    geneid[row['gene_symbol']] = row['ncbi_gene_id']
    phenotypes[row['gene_symbol']].add(row["hpo_name"])




for e in w.iterenzyme():
    gn = e.get('genename')
    species = e.get('species')
    #if species == 'Mouse':
    if species != 'Human':
       continue
        
    if len(phenotypes.get(gn,[])) > 0: 
       e.set("phenotype",(sorted(phenotypes[gn])))
       e.set("phenotype_source_key",(geneid[gn]))
       e.set("phenotype_source","HPO")
       print(e)
    else:	
       e.delete("phenotype")
       e.delete("phenotype_source_key")
       e.delete("phenotype_source")


    if w.put(e):
       print (gn)

