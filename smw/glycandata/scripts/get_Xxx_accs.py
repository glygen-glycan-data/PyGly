#!/bin/env python27

import sys

from getwiki import GlycanDataWiki
w = GlycanDataWiki()

Xxx_accs = set()

for m in w.iterglycan():
    try:
        comp_ann = list(m.annotations(property='XxxCount',type='MonosaccharideCount',source='EdwardsLab'))[0]
        Xxx_accs.add(comp_ann.get('id'))
    except:
        continue

wh = open('../data/Xxx_count.txt','w')
for acc in sorted(Xxx_accs):
    print >> wh, acc
wh.close()

        
