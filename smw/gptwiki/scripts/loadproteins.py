#!/bin/env python27

from getwiki import GPTWiki, Protein

import sys, urllib, string, csv
import Bio.SeqIO

w = GPTWiki()

seen = set()
for praccfile in sys.argv[1:]:
  for pracc in open(praccfile):
    pracc = pracc.strip()
    if pracc in seen:
	continue
    # print >>sys.stderr, pracc
    seen.add(pracc)
    data = urllib.urlopen('http://www.uniprot.org/uniprot/'+pracc+'.xml')
    for seq_record in Bio.SeqIO.parse(data,'uniprot-xml'):
	desc = seq_record.description
	pracc1 = seq_record.id
	seq = str(seq_record.seq)
	gene = seq_record.annotations['gene_name_primary']
        sys.stdout.write(seq_record.format('fasta'))
	break
    name = gene
    seqlines = []
    for i in range(0,len(seq),60):
	seqlines.append(seq[i:i+60])
    seq = "\n".join(seqlines)
    p = Protein(accession=pracc1,description=desc,gene=gene,sequence=seq,species='Homo sapiens',name=name)
    if w.put(p):
	print >>sys.stderr, pracc1
