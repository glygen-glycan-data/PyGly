#!/bin/env python27

import sys
from operator import itemgetter

from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

from pygly.GNOme import SubsumptionGraph

gnome = SubsumptionGraph()
gnome.loaddata(sys.argv[1])
sys.argv.pop(1)

debug = False
def iterglycan():
    if len(sys.argv) > 1:
	for acc in sys.argv[1:]:
	    m = w.get(acc)
	    if m:
		debug = True
		yield m
    else:
	for m in w.iterglycan():
	    yield m

ecstr = """
GTC	        GlyTouCan TaxID %s
UCKB	        UniCarbKB TaxID %s
noXyl		No Xyl
noXylNeuGc	No Xyl or NeuGc
sub		Subsumption of %s
"""
ecd = dict()
for l in ecstr.splitlines():
    sl = l.split(None,1)
    if len(sl) < 2:
	continue
    ecd[sl[0]] = sl[1].strip()

def ec(*args):
    key = args[0]
    if len(args) == 1:
        return ecd[key]
    rest = tuple(args[1:])
    return ecd[key]%rest

species = defaultdict(dict)

for m in iterglycan():
    acc = m.get('accession')

    print "Pass1:",acc

    try:
        gtctaxids = set(m.get_annotation_values(property="Taxonomy",source="GlyTouCan",type="Taxonomy"))
    except LookupError:
	gtctaxids = set()
    try:
        uckbtaxids = set(m.get_annotation_values(property="Taxonomy",source="UniCarbKB",type="Taxonomy"))
    except LookupError:
	uckbtaxids = set()

    comp = defaultdict(int)
    for an in m.annotations(type="MonosaccharideCount"):
	prop = an.get('property')
	if prop.endswith('Count'):
	    prop = prop[:-5]
	if prop == "Monosaccharide":
	    prop = "Count"
	comp[prop] = int(an.get('value'))

    evidence = set()
    for taxid in (9606,):
        if str(taxid) in gtctaxids:
            evidence.add(ec('GTC',taxid))
        if str(taxid) in uckbtaxids:
            evidence.add(ec('UCKB',taxid))
    if len(evidence) > 0 and comp['Xyl'] == 0 and  comp['NeuGc'] == 0:
	evidence.add(ec('noXylNeuGc'))

    if len(evidence) > 1 and ec('noXylNeuGc') in evidence:
        species[acc]['human'] = (True,True,sorted(evidence))
    else:
        species[acc]['human'] = (False,True,sorted(evidence))

    evidence = set()
    for taxid in (10090,):
        if str(taxid) in gtctaxids:
            evidence.add(ec('GTC',taxid))
        if str(taxid) in uckbtaxids:
            evidence.add(ec('UCKB',taxid))
    if len(evidence) > 0 and comp['Xyl'] == 0:
	evidence.add(ec('noXyl'))

    if len(evidence) > 1 and ec('noXyl') in evidence:
        species[acc]['mouse'] = (True,True,sorted(evidence))
    else:
        species[acc]['mouse'] = (False,True,sorted(evidence))

    evidence = set()
    for taxid in (10116,):
        if str(taxid) in gtctaxids:
            evidence.add(ec('GTC',taxid))
        if str(taxid) in uckbtaxids:
            evidence.add(ec('UCKB',taxid))
    if len(evidence) > 0 and comp['Xyl'] == 0:
	evidence.add(ec('noXyl'))

    if len(evidence) > 1 and ec('noXyl') in evidence:
        species[acc]['rat'] = (True,True,sorted(evidence))
    else:
        species[acc]['rat'] = (False,True,sorted(evidence))


for acc in sorted(species):

    print "Pass2:",acc

    any = False
    for sp in species[acc]:
        if species[acc][sp][0] and species[acc][sp][1]:
            any = True
            break

    if not any:
        continue
        
    for anc in gnome.ancestors(acc):

        if anc not in species:
            continue

        for sp in species[acc]:

            if species[acc][sp][0] and species[acc][sp][1]:

                isa,direct,evidence = species[anc][sp]

                evidence = sorted(set(list(evidence) + [ec('sub',acc)]))

                if not isa:
                    species[anc][sp] = (True,False ,evidence)
                else:
                    species[anc][sp] = (True,direct,evidence)

                    
for acc in sorted(species):

    print "Pass3:",acc

    m.delete_annotations(type="Species",source="EdwardsLab")

    m.set_annotation(property="Human Evidence",
                     value=species[acc]['human'][2],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Human",
                     value='true' if species[acc]['human'][0] else 'false',
                     type="Species",
                     source="EdwardsLab")

    m.set_annotation(property="Mouse Evidence",
                     value=species[acc]['mouse'][2],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Mouse",
                     value='true' if species[acc]['mouse'][0] else 'false',
                     type="Species",
                     source="EdwardsLab")
    
    m.set_annotation(property="Rat Evidence",
                     value=species[acc]['rat'][2],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Rat",
                     value='true' if species[acc]['rat'][0] else 'false',
                     type="Species",
                     source="EdwardsLab")
    
    if w.put(m):
        print acc

