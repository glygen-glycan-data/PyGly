#!/bin/env python3.12

import re, sys
from getwiki import GlycanData
from collections import defaultdict
 
headers = filter(None,"""
GlyTouCan AC
Tissue
Species
AnnotationCategory
Source
Source ID
""".splitlines())

tospecies = {
	"TissueInHuman": "Human",
        "TissueInMouse": "Mouse",
        "TissueInZebraFish": "ZebraFish",
        "TissueInZebrafish": "ZebraFish",
        "TissueInFruitFly": "FruitFly"
}

w = GlycanData()

print("\t".join(headers))
for acc in w.iterglycanid():
    g = w.get(acc)
    tissueannotations = list(g.annotations(type="Tissue"))
    for ann in tissueannotations:
        sp = tospecies[ann.get("property")]
        for vi in ann.get('value'):
            if ann.get('source') != "GNOme":
                print("\t".join([acc,vi,sp,"Direct",ann.get('source'),ann.get('sourceid')]))
            else:
                print("\t".join([acc,vi,sp,"Subsumption",ann.get('source'),ann.get('sourceid')]))
