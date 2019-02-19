#!/bin/env python27

from getwiki import GPTWiki, Protein

import sys, urllib, string, csv
import Bio.SeqIO

w = GPTWiki()
praccfile = sys.argv[1]
fastafile = praccfile.rsplit('.',1)[0]+'.fasta'
fastawh = open(fastafile,'w')

seen = set()
for r in csv.DictReader(open(praccfile),dialect='excel-tab'):
    pracc = r['ProteinName'].strip()
    if pracc in seen:
	continue
    # print >>sys.stderr, pracc
    seen.add(pracc)
    data = urllib.urlopen('http://www.uniprot.org/uniprot/'+pracc+'.xml')
    for seq_record in Bio.SeqIO.parse(data,'uniprot-xml'):
	desc = seq_record.description
	pracc1 = seq_record.id
	gene = seq_record.annotations['gene_name_primary']
        fastawh.write(seq_record.format('fasta'))
	break
    name = gene
    p = Protein(accession=pracc1,description=desc,gene=gene,species='Homo sapiens',name=name)
    if w.put(p):
	print >>sys.stderr, pracc1
    
fastawh.close()
