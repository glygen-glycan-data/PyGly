#!/bin/env python27

import sys

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

import findpygly
from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

current = set()

N_Glycan = ['G00026MO','G00028MO','G00029MO','G00030MO']

O_Glycan = ['G00031MO','G00032MO','G00033MO','G00034MO','G00035MO','G00036MO','G00037MO',
        'G00038MO','G00039MO','G00040MO','G00041MO','G00042MO','G00043MO','G00044MO']
	
subtype_dic = {'G00026MO':'core basic','G00028MO':'high mannose', 'G00029MO':'hybrid',
        'G00030MO':'complex','G00031MO':'core 1','G00032MO':'core 1 fuzzy','G00033MO':'core 2',
	'G00034MO':'core 2 fuzzy','G00035MO':'core 3','G00036MO':'core 3 fuzzy','G00037MO':'core 4',
	'G00038MO':'core 4 fuzzy','G00039MO':'core 5','G00040MO':'core 5 fuzzy','G00041MO':'core 6',
	'G00042MO':'core 6 fuzzy','G00043MO':'core 7','G00044MO':'core 7 fuzzy'}
	
for m in w.iterglycan():
    acc = m.get('accession')
    subtype_list = []
    count = 0    
    for motif in list(gtc.getmotif(acc)):
        count += 1
        acc, name = motif
	if acc in subtype_dic:
            subtype_list.append(subtype_dic[acc])
            m.add_annotation(value=subtype_list,
                    property='GlycanSubtype',
                    source='EdwardsLab', type='Classification')	    
        if acc in N_Glycan and count==1:
            m.add_annotation(value='N-linked',
                    property='GlycanType',
                    source='EdwardsLab', type='Classification')
        elif acc in O_Glycan and count==1:
            m.add_annotation(value='O-linked',
                    property='GlycanType',
                    source='EdwardsLab', type='Classification')
for m in w.iterglycan():
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:",m.get('pagename')
        w.delete(m.get('pagename'))