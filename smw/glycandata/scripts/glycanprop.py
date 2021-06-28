#!/bin/env python2

import re, sys
from getwiki import GlycanData
 
headers = """
glytoucan_acc glytoucan_type glycan_mass glycan_permass
base_composition composition topology monosaccharides
""".split()

w = GlycanData()

print "\t".join(headers)
for acc in w.iterglycanid():
    g = w.get(acc)
    row = dict(glytoucan_acc=acc)
    for key,prop in zip(['glytoucan_type','base_composition','composition','topology'],
                        ['Level','BaseComposition','Composition','Topology']):
        try:
            row[key] = g.get_annotation_value(property=prop,
                                              type="Subsumption",
                                              source='GNOme')
        except LookupError:
            pass

    for key,prop in zip(['glycan_mass','glycan_permass'],['UnderivitizedMW','PermethylatedMW']):
        try:
            row[key] = g.get_annotation_value(property=prop,
                                              type="MolWt",
                                              source='EdwardsLab')
        except LookupError:
            pass

    for key,prop in zip(['monosaccharides'],['MonosaccharideCount']):
        try:
            row[key] = g.get_annotation_value(property=prop,
                                              type="MonosaccharideCount",
                                              source='EdwardsLab')
        except LookupError:
            pass
        
    print "\t".join(map(lambda h: row.get(h,""),headers))
