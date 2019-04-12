#!/bin/env python27

import sys

from getwiki import GlycanDataWiki
w = GlycanDataWiki()

u_mw = set()

for m in w.iterglycan():
    gtcmw_ann = list(m.annotations(property='UnderivitizedMW',type='MolWt',source='GlyTouCan'))[0]
    try:
        umw_ann = list(m.annotations(property='UnderivitizedMW',type='MolWt',source='EdwardsLab'))[0]
    except:   
        u_mw.add(gtcmw_ann.get('id'))
        continue

umw_wh = open('../data/missing_umw.txt','w')
for acc in sorted(u_mw):
    print >> umw_wh,acc
umw_wh.close()
        
