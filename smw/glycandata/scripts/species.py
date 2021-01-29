#!/bin/env python27

import re, sys
from getwiki import GlycanData
from collections import defaultdict
 
headers = filter(None,"""
GlyTouCan AC
Species Name
AnnotationCategory
Source
Source ID
tax_id
""".splitlines())

w = GlycanData()

print "\t".join(headers)
for acc in w.iterglycanid():
    g = w.get(acc)
    validspecies = set()
    speciesannotations = list(g.annotations(type="Species",source="EdwardsLab"))
    for ann in speciesannotations:
	prop = ann.get('property')
	if 'Category' in prop:
	    species_name = prop.rsplit(None,1)[0]
	    if ann.get('value') in ("Direct","Subsumption"):
	        validspecies.add(species_name)
    for ann in speciesannotations:
	prop = ann.get('property')
	if 'Evidence' in prop:
	    species_name = prop.rsplit(None,1)[0]
	    if species_name not in validspecies:
		continue
	    for value in sorted(ann.get('value',[])):
	        if 'Subsumption of ' in value:
		    via = value.rsplit(None,1)[1]
		    print "\t".join([acc,species_name,"Subsumption","GNOme",via,""])
		elif ' TaxID ' in value:
		    source = value.split(None,1)[0]
		    taxid = value.rsplit(None,1)[1]
		    if ':' in source:
			source,sourceid = source.split(':')
		    else:
			sourceid = ""
		    print "\t".join([acc,species_name,"Direct",source,sourceid,taxid])
