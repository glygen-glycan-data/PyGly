#!/bin/env python27

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
    row = defaultdict(int)
    row['accession'] = acc
    for ann in g.annotations(type="MonosaccharideCount",source="EdwardsLab"):
        value = int(ann.get('value'))
        prop = ann.get('property')
        if prop.endswith('Count'):
            prop = prop[:-5]
        row[prop] = value    
    row['Count'] = row['Monosaccharide']
    row['Xxx'] = row['Count'] - sum(map(lambda k: row[k],['Hex','HexNAc','dHex','NeuAc','NeuGc','HexA','HexN']))
    row['X'] += row['Me'] # S and P are accounted for...
    print "\t".join(map(lambda h: str(row.get(h,0)),headers))
