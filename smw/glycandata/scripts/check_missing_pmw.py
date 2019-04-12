#!/bin/env python27

import sys

from getwiki import GlycanDataWiki
w = GlycanDataWiki()

p_mw = set()

for m in w.iterglycan():
    gtcmw_ann = list(m.annotations(property='UnderivitizedMW',type='MolWt',source='GlyTouCan'))[0]
    try:
        pmw_ann = list(m.annotations(property='PermethylatedMW',type='MolWt',source='EdwardsLab'))[0]
    except:   
        p_mw.add(gtcmw_ann.get('id'))
        continue

pmw_wh = open('../data/missing_pmw.txt','w')
for acc in sorted(p_mw):
    print >> pmw_wh,acc
pmw_wh.close()
        
