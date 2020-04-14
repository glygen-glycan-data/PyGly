#!/bin/env python27

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
    for prop in ('Ancestor','Descendant','SubsumedBy','Subsumes'):
        for ann in g.annotations(property=prop,type="Subsumption",source="GNOme"):
	    for value in sorted(ann.get('value',[])):
	        row = dict(GlyTouCanAccession=acc,Relationship=prop,RelatedAccession=value)
                print "\t".join(map(lambda h: row.get(h,""),headers))
