#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

phenotypes = defaultdict(list)
mgiacc = dict()
# "Gene Symbol","MGI Gene Id","# Phenotype Hits","Phenotype Hits"
for row in csv.DictReader(gzip.open(sys.argv[1])):
    phenotypes[row['Gene Symbol']].extend(filter(None,[p.strip() for p in row["Phenotype Hits"].split("::")]))
    mgiacc[row['Gene Symbol']] = row["MGI Gene Id"]

for e in w.iterenzyme():
    gn = e.get('genename')
    if e.get('species') != 'Mouse':
       continue
    if len(phenotypes.get(gn,[])) > 0:
       e.set("phenotype",sorted(phenotypes[gn]))
       e.set("phenotype_source","IMPC")
       e.set("phenotype_source_key",mgiacc[gn])
    else:
       e.delete("phenotype")
       e.delete("phenotype_source")
       e.delete("phenotype_source_key")
    if w.put(e):
        print gn
