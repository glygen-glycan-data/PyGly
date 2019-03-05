#!/bin/env python27

import sys

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

import findpygly
from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

current = set()
for gtcacc in open(sys.argv[1]):
    gtcacc = gtcacc.strip()
    g = Glycan(accession=gtcacc,
               iupac=gtc.getseq(gtcacc,'iupac_extended'))
    g.add_annotation(value=gtc.getmass(gtcacc),
                     property='UnderivitizedMW',
                     source='GlyTouCan',type='MolWt')
    if gtcacc == 'G00031MO':
	g.add_annotation(value='O-linked',
			 property='GlycanType',
		 	 source='EdwardsLab',
			 type='Classification',
			 method='Glycan Type by Motif Match',
			 reference='https://glytoucan.org/Structures/Glycans/G00032MO')
	g.add_annotation(value='core 1',
			 property='GlycanSubtype',
		 	 source='EdwardsLab',
			 type='Classification',
			 method='Glycan Type by Motif Match',
			 reference='https://glytoucan.org/Structures/Glycans/G00032MO')
    if w.put(g):
	print >>sys.stderr, g.get('accession')
    current.add(gtcacc)

for m in w.iterglycan():
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))
