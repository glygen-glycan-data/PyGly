#!/bin/env python27

import sys
import urllib

acclist = sys.argv[1:]
for acc in acclist:
    link = "https://www.uniprot.org/uniprot/" + acc + ".fasta"
    lines = urllib.urlopen(link).readlines()
    seqlist = lines[1:]
    sequence = ""
    for i in seqlist:
        sequence += str(i)
    accession = lines[0][1:].split('|')[1]
    description = lines[0].split(' ', 1)[1].split('=')[0][:-3]
    species = lines[0].split(' ', 1)[1].split('=')[1][:-3]
    gene = lines[0].split(' ', 1)[1].split('=')[3][:-3]

    f = open("../wiki/pages/" + acc + ".txt", "w")
    f.write("{{Protein\n")
    f.write("|accession=" + accession + "\n")
    f.write("|description=" + description + "\n")
    f.write("|gene=" + gene + "\n")
    f.write("|sequence=" + sequence)
    f.write("|species=" + species + "\n}}")
    f.close

