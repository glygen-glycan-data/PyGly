#!/bin/env python27

import sys
from operator import itemgetter

from getwiki import GlycanData, Glycan

from pygly.Monosaccharide import SuperClass, Stem, Config, Mod, Substituent

from pygly.GlyNLinkedFilter import GlyNLinkedFilter
mnlc = GlyNLinkedFilter(None).test1

w = GlycanData()

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

def gal(m):
    if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
        if not m.has_mods() and not m.has_substituents():
	    return True
    return False

def glcnac(m):
    if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
        if not m.has_mods() and m.is_nacetylated():
	    return True
    return False

def galnac(m):
    if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
        if not m.has_mods() and m.is_nacetylated():
	    return True
    return False

def fuc(m):
    if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,) and m.config() == (Config.l,):
        if m.count_mod(Mod.d) == 1 and m.count_mod() == 1:
	    return True
    return False

def sia(m):
    if m.superclass() == SuperClass.NON and m.stem() == (Stem.gro,Stem.gal):
	if set(map(itemgetter(1),m.mods())) == set([Mod.a,Mod.keto,Mod.d]):
	    found = False
	    for s in m.substituents():
	        if s.name() in (Substituent.nAcetyl,Substituent.nglycolyl):
		    found = True
		    break
	    if len(m.substituents()) == 1 and found:
		return True
    return False

def acceptcore1(g):
    r = g.root()
    galseen = False
    for c in r.children():
	if gal(c) and not galseen:
	    galseen = True
	    continue
	if fuc(c) or sia(c):
	    continue
	return False
    return True

def acceptcore36(g):
    r = g.root()
    glcnacseen = False
    for c in r.children():
	if glcnac(c) and not glcnacseen:
	    glcnacseen = True
	    continue
	if fuc(c) or sia(c):
	    continue
	return False
    return True

def acceptcore57(g):
    r = g.root()
    galnacseen = False
    for c in r.children():
	if galnac(c) and not galnacseen:
	    galnacseen = True
	    continue
	if fuc(c) or sia(c):
	    continue
	return False
    return True

debug = False
def iterglycan():
    global debug
    if len(sys.argv) > 1:
	for acc in sys.argv[1:]:
	    m = w.get(acc)
	    if m:
		debug = True
		yield m
    else:
	for m in w.iterglycan():
	    yield m
    
for m in iterglycan():
    acc = m.get('accession')

    classifications = set()

    try:
        count = int(m.get_annotation_value('MonosaccharideCount',source='EdwardsLab'))
    except LookupError:
	count = 0

    try:
        mancount = int(m.get_annotation_value('ManCount',source='EdwardsLab'))
    except LookupError:
	mancount = 0

    try:
        glcnaccount = int(m.get_annotation_value('GlcNAcCount',source='EdwardsLab'))
    except LookupError:
	glcnaccount = 0

    try:
        fuccount = int(m.get_annotation_value('FucCount',source='EdwardsLab'))
    except LookupError:
	fuccount = 0

    try:
        xxxcount = int(m.get_annotation_value('XxxCount',source='EdwardsLab'))
    except LookupError:
	xxxcount = 0

    m.delete_annotations(source='EdwardsLab', type='Classification')

    try:
        motifann = list(m.annotations(property='Motif',type='Motif',source='GlyTouCan'))[0]
    except IndexError:
	if not debug:
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
        if thetype == "O-linked" and set(subtypes) == set(["core 2", "core 1","core 6","core 3"]):
            subtypes = ["core 2"]
        if thetype == "O-linked" and set(subtypes) == set(["core 2", "core 1","core 6"]):
            subtypes = ["core 2"]
        if thetype == "O-linked" and set(subtypes) == set(["core 4","core 3","core 6"]):
            subtypes = ["core 4"]

    # Make various checks on the details of the subtypes and types values. 
    if len(types) > 1:
	print >>sys.stderr, "Glycan %s has more than one type: %s."%(acc,", ".join(sorted(types)))
	continue
    if len(types) == 1 and len(subtypes) > 1:
	print >>sys.stderr, "Glycan %s has more than one subtype: %s."%(acc,", ".join(sorted(subtypes)))
	subtypes = []

    # Check for additional issues...
    if types == ["N-linked"] and subtypes == ["high mannose"] and count > 0 and (mancount + 2) != count:
        subtypes = []

    # Paucimannose
    if subtypes == [] and count <= 6:
	g = m.getGlycan()
	if (g and not g.undetermined() and mnlc(g) and count <= 4) or (types == ["N-linked"] and count in (5,6)):
	    if count == (mancount + glcnaccount + fuccount) and glcnaccount == 2 and (count == 5 and fuccount == 0) or (count == 6 and fuccount == 1):
	        types = ["N-linked"]
	        subtypes = ['paucimannose']

    if types == ["O-linked"] and "core 1" in subtypes:
	g = m.getGlycan()
	if g and not acceptcore1(g):
	    subtypes.remove("core 1")

    if types == ["O-linked"] and "core 3" in subtypes:
	g = m.getGlycan()
	if g and not acceptcore36(g):
	    subtypes.remove("core 3")

    if types == ["O-linked"] and "core 5" in subtypes:
	g = m.getGlycan()
	if g and not acceptcore57(g):
	    subtypes.remove("core 5")

    if types == ["O-linked"] and "core 6" in subtypes:
	g = m.getGlycan()
	if g and not acceptcore36(g):
	    subtypes.remove("core 6")

    if types == ["O-linked"] and "core 7" in subtypes:
	g = m.getGlycan()
	if g and not acceptcore57(g):
	    subtypes.remove("core 7")

    if len(types) > 0:
        m.set_annotation(value=types[0], property='GlycanType',
                         source='EdwardsLab', type='Classification')
        if len(subtypes) > 0:
            m.set_annotation(value=subtypes[0], property='GlycanSubtype',
                             source='EdwardsLab', type='Classification')
    if not debug:
	if w.put(m):
	    print acc
    else:
	    print m
