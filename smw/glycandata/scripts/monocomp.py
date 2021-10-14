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
    row['Xxx'] = [0,False]
    for k in row:
	if k.startswith('Man') or k.startswith('Gal') or k.startswith('Glc') or k.startswith('Fuc'):
	    continue
	if k in ('accession','Sia','S','P','Me','Count','Monosaccharide','aldi','IdoA','X','Xxx'):
	    continue
	if k in ['Hex','HexNAc','dHex','NeuAc','NeuGc','HexA','HexN']:
	    continue
        assert k in ['dHex+aldi','Hex+aldi','HexNAc+aldi','Pent','Xyl'], "Unexpected Xxx monosaccharide: "+k
        # print(k)
	row['Xxx'][0] += row[k][0]
	row['Xxx'][1] |= row[k][1]
    # S and P are accounted for
    row['X'][0] += row['Me'][0]
    row['X'][1] |= row['Me'][1]
    print "\t".join(map(lambda h: "%s%s"%(row[h][0],"+" if row[h][1] else ""),headers))
