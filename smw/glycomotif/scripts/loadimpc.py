#!/bin/env python2

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

phenotypes = defaultdict(list)
# "Gene Symbol","MGI Gene Id","# Phenotype Hits","Phenotype Hits"
for row in csv.DictReader(gzip.open(sys.argv[1])):
    phenotypes[row['Gene Symbol']].extend([p.strip() for p in row["Phenotype Hits"].split("::")])

for e in w.iterenzyme():
    gn = e.get('genename')
    if len(phenotypes.get(gn,[])) > 0:
       e.set("phenotype",sorted(phenotypes[gn]))
    else:
       e.delete("phenotype")
    if w.put(e):
        print gn
