#!/bin/env python27

import sys, time
from collections import defaultdict
import csv

from getwiki import GlycanData, Glycan
w = GlycanData()

taxids={'mouse': '10090', 'human': '9606'}

cr2enz = defaultdict(list)
headers = None
for i,r in enumerate(csv.reader(open(sys.argv[1]))):
    if i == 0:
	headers = map(str.lower,r)
	newheaders = []
	for h in headers:
	    if h not in newheaders:
		newheaders.append(h)
		continue
	    i = 1
	    h1 = h+str(i)
	    while h1 in newheaders:
		i += 1
		h1 = h+str(i)
	    newheaders.append(h1)
	headers = newheaders
    else:
	row = dict(zip(headers,r))
	# Uniprot,Protein RefSeq,DNA RefSeq,gene name,gene id,species
	for tag in ("","1","2","3"):
	    if row.get('species'+tag) in ("N/A",None,""):
		continue
	    if row.get('species'+tag) not in taxids:
		continue
	    if row.get('uniprot'+tag) in ("N/A",None,""):
		continue
	    species = row['species'+tag].strip()
	    cr2enz[(row['residue_name'],species)].append(dict(uniprot=row['uniprot'+tag].strip(),
				                         prrefseq=row['protein refseq'+tag].strip(),
					                 mrnarefseq=row['dna refseq'+tag].strip(),
					                 genename=row['gene name'+tag].strip(),
					                 geneid=row['gene id'+tag].strip(),
						         species=row['species'+tag].strip()))

gtc2cr = defaultdict(set)
for r in csv.DictReader(open(sys.argv[2])):
    # glycan_ID,order,residue
    gtc2cr[r['glycan_ID'].strip()].add(r['residue'].strip())

for g in w.iterglycan():
    start = time.time()
    gtc = g.get('accession')

    g.delete_annotations(type='Enzyme',source='GlycO')
    enzymes = set()
    for species in taxids:
	genes = set()
        for cr in gtc2cr[gtc]:
            for enz in cr2enz[(cr,species)]:
	        genes.add(enz['genename'])
                enzymes.add(":".join([enz['uniprot'],enz['genename'],enz['geneid'],taxids[enz['species']]]))
        g.set_annotation(value=list(genes),property='%s Enzyme'%(species.title()),type='Enzyme',source='GlycO')
    g.set_annotation(value=list(enzymes),property='EnzymeDetails',type='Enzyme',source='GlycO')

    if w.put(g):
        print >>sys.stderr, "%s updated in %.2f sec"%(gtc,time.time()-start,)
    else:
        print >>sys.stderr, "%s checked in %.2f sec"%(gtc,time.time()-start,)
