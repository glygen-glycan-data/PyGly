#!/bin/env python2

import re, sys
from getwiki import GlycanData
from collections import defaultdict
 
headers = """
GlyTouCanAccession Relationship RelatedAccession
""".split()

w = GlycanData()

print "\t".join(headers)
for acc in w.iterglycanid():
    g = w.get(acc)
    fullydets = set()
    for prop in sorted(['Ancestor','Descendant','SubsumedBy','Subsumes','BaseComposition','Composition','Topology','MissingScore','FullyDetermined','Level']):
        for ann in g.annotations(property=prop,type="Subsumption",source="GNOme"):
	    if prop in ('Composition','Topology','BaseComposition','MissingScore','Level'):
		value = ann.get('value',None)
		if value:
	            row = dict(GlyTouCanAccession=acc,Relationship=prop,RelatedAccession=value)
                    print "\t".join(map(lambda h: row.get(h,""),headers))
	    else:
	        for value in sorted(ann.get('value',[])):
	            row = dict(GlyTouCanAccession=acc,Relationship=prop,RelatedAccession=value)
                    print "\t".join(map(lambda h: row.get(h,""),headers))
