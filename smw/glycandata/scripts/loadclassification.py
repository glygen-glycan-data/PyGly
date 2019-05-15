#!/bin/env python27

import sys
from operator import itemgetter

from getwiki import GlycanDataWiki, Glycan
w = GlycanDataWiki()

motif_rules_data = """
G00026MO	N-linked	
G00028MO	N-linked	high mannose
G00029MO	N-linked	hybrid
G00030MO	N-linked	complex
G00031MO	O-linked	core 1
G00032MO	O-linked	core 1
G00033MO	O-linked	core 2
G00034MO	O-linked	core 2
G00035MO	O-linked	core 3
G00036MO	O-linked	core 3
G00037MO	O-linked	core 4
G00038MO	O-linked	core 4
G00039MO	O-linked	core 5
G00040MO	O-linked	core 5
G00041MO	O-linked	core 6
G00042MO	O-linked	core 6
G00043MO	O-linked	core 7
G00044MO	O-linked	core 7
"""

motifrules = dict()
for l in motif_rules_data.splitlines():
    if not l.strip():
	continue
    sl = l.split('\t')
    assert len(sl) == 3
    motifrules[sl[0]] = (sl[1],sl[2])

for m in w.iterglycan():
    acc = m.get('accession')

    classifications = set()

    count = 0; mancount = 0; xxxcount = 0;
    for ann in m.annotations(type='MonosaccharideCount',source='EdwardsLab'):
	if ann.get('property') == "MonosaccharideCount":
	    count = int(ann.get('value',[0])[0])
	elif ann.get('property') == "ManCount":
	    mancount = int(ann.get('value',[0])[0])
	elif ann.get('property') == "XxxCount":
	    xxxcount = int(ann.get('value',[0])[0])

    m.delete_annotations(source='EdwardsLab', type='Classification')

    try:
        motifann = list(m.annotations(property='Motif',type='Motif',source='GlyTouCan'))[0]
    except IndexError:
	w.put(m)
	continue

    for motifacc in motifann.get('value',[]):
	if motifacc in motifrules:
	    classifications.add(motifrules[motifacc])

    # Add logic for multi-classifications...
    types = list(set(map(itemgetter(0),classifications)))
    subtypes = list(set(filter(None,map(itemgetter(1),classifications))))

    if len(types) == 1 and len(subtypes) > 1:
	# resolve these as needed!
	thetype = types[0]
	if thetype == "N-linked" and set(subtypes) == set(["hybrid","high mannose"]):
	    subtypes = ["hybrid"]
        if thetype == "O-linked" and set(subtypes) == set(["core 2", "core 1","core 6"]):
            subtypes = ["core 2"]
        if thetype == "O-linked" and set(subtypes) == set(["core 4","core 3","core 6"]):
            subtypes = ["core 4"]

    # Make various checks on the details of the subtypes and types values. 
    if len(types) > 1:
	print >>sys.stderr, "Glycan %s has more than one type: %s."%(acc,", ".join(sorted(types)))
	continue
    if len(subtypes) > 1:
	print >>sys.stderr, "Glycan %s has more than one subtype: %s."%(acc,", ".join(sorted(subtypes)))
        continue

    # Check for additional issues...
    if types == ["N-linked"] and subtypes == ["high mannose"] and count > 0 and (mancount + 2) != count:
        subtypes = []

    if len(types) > 0:
        m.set_annotation(value=types[0], property='GlycanType',
                         source='EdwardsLab', type='Classification')
        if len(subtypes) > 0:
            m.set_annotation(value=subtypes[0], property='GlycanSubtype',
                             source='EdwardsLab', type='Classification')
    if w.put(m):
	print acc
