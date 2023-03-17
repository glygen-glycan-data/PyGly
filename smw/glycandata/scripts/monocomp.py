#!/bin/env python2

import re, sys
from getwiki import GlycanData
from collections import defaultdict
 
headers = """
accession Hex HexNAc dHex NeuAc NeuGc HexA HexN S P aldi Xxx X Count
""".split()

w = GlycanData()

print "\t".join(headers)
for acc in w.iterglycanid():
    g = w.get(acc)
    row = defaultdict(lambda: [0,False])
    row['accession'] = (acc,False)
    for ann in g.annotations(type="MonosaccharideCount",source="EdwardsLab"):
        try:
            value = [int(ann.get('value')),False]
        except ValueError:
            value = [int(ann.get('value')[:-1]),True]
        prop = ann.get('property')
        if prop.endswith('Count'):
            prop = prop[:-5]
        row[prop] = value    
    # print row
    row['Count'] = row['Monosaccharide']
    if 'Xxx' not in row:
        row['Xxx'] = [0,False]
    for k in list(row):
        if row[k][0] == 0:
            continue
	if k.startswith('Man') or k.startswith('Gal') or k.startswith('Glc') or k.startswith('Fuc'):
	    continue
	if k in ('accession','S','P','Me','Count','Monosaccharide','aldi','X','Xxx','IdoA'):
	    continue
	if k in ['Hex','HexNAc','dHex','NeuAc','NeuGc','HexA','HexN','Hex+aldi','dHex+aldi','HexNAc+aldi']:
	    continue
        assert k in ['Pent','Xyl','Sia','Kdn'], "Unexpected Xxx monosaccharide: "+k
        if k == 'Sia':
            extraxxx = row[k][0]-sum(map(lambda k: row[k][0],"NeuAc NeuGc Kdn".split()))
            if extraxxx > 0:
	        row['Xxx'][0] += extraxxx
                if not row['Kdn'][1] and not row['NeuAc'][1] and not row['NeuGc'][1]:
	            row['Xxx'][1] |= row[k][1]
            continue
        if k == 'Pent':
            extraxxx = row[k][0]-sum(map(lambda k: row[k][0],"Xyl".split()))
            if extraxxx > 0:
	        row['Xxx'][0] += extraxxx
                if not row['Xyl'][1]:
	            row['Xxx'][1] |= row[k][1]
            continue
        assert k in ['Xyl','Kdn'], "Unexpected Xxx monosaccharide: "+k
	row['Xxx'][0] += row[k][0]
	row['Xxx'][1] |= row[k][1]
    # S and P are accounted for
    row['X'][0] += row['Me'][0]
    row['X'][1] |= row['Me'][1]
    if row['Count'] == 0:
        continue
    if sum(map(lambda k: row[k][0],"Hex HexNAc dHex NeuAc NeuGc HexA HexN Xxx".split())) != row['Count'][0]:
        continue
    print "\t".join(map(lambda h: "%s%s"%(row[h][0],"+" if row[h][1] else ""),headers))
