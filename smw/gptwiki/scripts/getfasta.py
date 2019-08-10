#!/bin/env python27

import sys, re
import urllib

acclist = sys.argv[1:]
for acc in acclist:
    link = "https://www.uniprot.org/uniprot/" + acc + ".fasta"
    lines = urllib.urlopen(link).readlines()

    desc = lines[0][1:]
    sequence = "".join(lines[1:])

    accessions,rest = desc.split(None,1)
    
    accession = accessions.split('|')[1]

    assert acc == accession

    restsplit = re.split(r' ([A-Z]+)=',rest)
    
    description = restsplit[0].strip()

    keyvaluepairs = dict()
    for i in range(1,len(restsplit),2):
        keyvaluepairs[restsplit[i]] = restsplit[i+1].strip()

    species = keyvaluepairs.get('OS')
    gene = keyvaluepairs.get('GN')
    
    f = open(acc + ".txt", "w")
    f.write("{{Protein\n")
    f.write("|accession=" + accession + "\n")
    f.write("|description=" + description + "\n")
    if gene:
        f.write("|gene=" + gene + "\n")
        f.write("|name=" + gene + "\n")
    f.write("|sequence=" + sequence)
    if species:
        f.write("|species=" + species + "\n}}")
    f.write('\n')
    f.close

