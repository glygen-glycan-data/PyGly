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
	       wurcs=gtc.getseq(gtcacc,'wurcs'),
	       glycoct=gtc.getseq(gtcacc,'glycoct'),
               iupac=gtc.getseq(gtcacc,'iupac_extended'))
    g.add_annotation(value=gtc.getmass(gtcacc),
                     property='UnderivitizedMW',
                     source='GlyTouCan',type='MolWt')
    try:
        mw = gtc.getGlycan(gtcacc).underivitized_molecular_weight()
        g.add_annotation(value=mw,
			 property='UnderivitizedMW',
			 source='EdwardsLab', type='MolWt')
    except (TypeError, ValueError):
        continue
    try:
        pmw = gtc.getGlycan(gtcacc).permethylated_molecular_weight()
        g.add_annotation(value=pmw,
			 property='PermethylatedMW',
			 source='EdwardsLab', type='MolWt')
    except (TypeError, ValueError):
        continue
    g.add_annotation(value=gtc.getmonocount(gtcacc),
		             property='MonosaccharideCount',
		             source='GlyTouCan',type='MonosaccharideCount')
    try: 
        comp = gtc.getGlycan(gtcacc).iupac_composition()
	for ckey in comp.keys():
            count = comp[ckey]
            if count > 0:
		if ckey=='Count':
		    g.add_annotation(value=count,
		         property='MonosaccharideCount',
		         source='EdwardsLab',type='MonosaccharideCount')
		else:
	            g.add_annotation(value=count,
		        property=ckey+'Count',
		        source='EdwardsLab',type='MonosaccharideCount')
    except (TypeError, ValueError):
        continue
    
    dic = {}    
    for xref in gtc.getcrossrefs(gtcacc):
        ref, c = xref.split(":")
	dic.setdefault(ref,[]).append(c)
    for key in dic:		
	g.add_annotation(value=dic[key],
		property=key,
		source='GlyTouCan',type='CrossReference')					 
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
