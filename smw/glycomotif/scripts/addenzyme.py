#!/bin/env python27

from getwiki import GlycoMotifWiki, Enzyme
import sys, re, glob, json

w = GlycoMotifWiki()

allenz = dict()

for jf in glob.glob(sys.argv[1]+'/*.json'):
    enzdata = json.loads(open(jf).read())
    for r in enzdata["residues"]:
        if r.get('glycotree',"none") == "none":
            continue
        for e in r['enzymes']:
            if e['gene_name'] in allenz:
                continue
            allenz[e['gene_name']] = dict(genename=e['gene_name'],
                                          species=('Human' if e['species'] == 'Homo sapiens' else 'Mouse'),
                                          uniprot=e['uniprot'])

for ed in allenz.values():
    enz = Enzyme(**ed)
    if w.put(enz):
        print("Enzyme %s updated."%(ed['genename'],))
