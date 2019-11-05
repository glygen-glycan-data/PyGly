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

def iterglycan():
    if len(sys.argv) > 1:
	seen = set()
	for acc in sys.argv[1:]:
	    if acc in seen:
		continue
	    m = w.get(acc)
	    if m:
		seen.add(acc)
		yield m
	    for desc in gnome.descendants(acc):
		if desc in seen:
		    continue
		m = w.get(desc)
		if m:
		    seen.add(desc)
		    yield m
		
    else:
	for m in w.iterglycan():
	    yield m

ecstr = """
GTC	        GlyTouCan TaxID %s
UCKB	        UniCarbKB:%s TaxID %s
noXyl		No Xyl
noXylNeuGc	No Xyl or NeuGc
sub		Subsumption of %s
compcl		Composition closure of %s via %s
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
    print >>sys.stderr, "Pass 1:", acc

    try:
        gtctaxids = set(m.get_annotation_values(property="Taxonomy",source="GlyTouCan",type="Taxonomy"))
    except LookupError:
	gtctaxids = set()

    uckbtaxids = defaultdict(set)
    for an in m.annotations(property="Taxonomy",source="UniCarbKB",type="Taxonomy"):
	sid = an.get('sourceid')
	if sid:
	    uckbtaxids[sid].update(an.get('value'))

    comp = defaultdict(int)
    for an in m.annotations(type="MonosaccharideCount"):
	prop = an.get('property')
	if prop.endswith('Count'):
	    prop = prop[:-5]
	if prop == "Monosaccharide":
	    prop = "Count"
	comp[prop] = int(an.get('value'))

    evidence = set()
    for taxid in map(str,(9606,)):
        if taxid in gtctaxids:
            evidence.add(ec('GTC',taxid))
	for sid in uckbtaxids:
	    if taxid in uckbtaxids[sid]:
                evidence.add(ec('UCKB',sid,taxid))
    if len(evidence) > 0 and comp['Xyl'] == 0 and  comp['NeuGc'] == 0:
	evidence.add(ec('noXylNeuGc'))

    if len(evidence) > 1 and ec('noXylNeuGc') in evidence:
        species[acc]['human'] = ('Direct',sorted(evidence))
    else:
        species[acc]['human'] = ('False',sorted(evidence))

    evidence = set()
    for taxid in map(str,(10090,)):
        if taxid in gtctaxids:
            evidence.add(ec('GTC',taxid))
	for sid in uckbtaxids:
	    if taxid in uckbtaxids[sid]:
                evidence.add(ec('UCKB',sid,taxid))
    if len(evidence) > 0 and comp['Xyl'] == 0:
	evidence.add(ec('noXyl'))

    if len(evidence) > 1 and ec('noXyl') in evidence:
        species[acc]['mouse'] = ('Direct',sorted(evidence))
    else:
        species[acc]['mouse'] = ('False',sorted(evidence))

    evidence = set()
    for taxid in map(str,(10116,10114)):
        if taxid in gtctaxids:
            evidence.add(ec('GTC',taxid))
	for sid in uckbtaxids:
	    if taxid in uckbtaxids[sid]:
                evidence.add(ec('UCKB',sid,taxid))
    if len(evidence) > 0 and comp['Xyl'] == 0:
	evidence.add(ec('noXyl'))

    if len(evidence) > 1 and ec('noXyl') in evidence:
        species[acc]['rat'] = ('Direct',sorted(evidence))
    else:
        species[acc]['rat'] = ('False',sorted(evidence))
    
    hcvspecies_descendents = set([11103,33745, 31653, 356418, 44022, 44023,
        128819, 128820, 128821, 693426, 693427, 745709, 1094898, 1094899,
        1094900, 1208062, 33746, 31654, 356390, 356419, 40271, 31649, 11113,
        356411, 31650, 11115, 356412, 31651, 356413, 40361, 44021, 356387,
        56388, 356466, 356414, 693428, 693429, 760561, 1094901, 41856,
        31646, 11104, 11108, 63746, 31647, 11105, 11116, 31642, 31645, 333284,
        20174, 421877, 421879, 31648, 356386, 356410, 40360, 44019, 44020,
        484894, 42182, 31655, 356391, 356420, 356427, 356423, 356465, 356421,
        356467, 356424, 356468, 356422, 356469, 356425, 378506, 378507,
        413255, 413256, 413257, 438880, 438881, 467337, 467338, 467339,
        467354, 569610, 1193076, 356113, 11106, 11107, 11109, 11110, 11111,
        11112, 11114, 11117, 31643, 31644, 36390, 329389, 357985, 357986,
        357987, 357988, 357989, 357990, 357991, 668553, 668554, 668555,
        758864, 761189, 761190, 761191, 761192, 1173523, 1208063, 356114,
        2791, 357355, 42792, 356115, 356116, 356417, 356426, 356415, 356416,
        759939, 1094895, 1405107, 578319, 487624, 578303, 578306, 595609,
        595610, 595611, 595612, 679182, 761967, 945057, 945058, 945059,
        45060, 945061, 945062, 945063, 945064, 945065, 945066, 945067,
        945068, 945069, 945070, 945071, 945072, 945073, 1006431, 1006432,
        1006433, 1006434, 1006435, 1006436, 1006437, 1006438, 1006439,
        1006440, 1006441, 1006442, 1006443, 1006444, 1053139, 1053140,
        053141, 1053142, 1053143, 1053144, 1053145, 1053146, 1053147,
        1259832, 1406420, 1406421, 1406422, 1406423, 1406424, 1406425,
        1406426, 1406427])

    evidence = set()
    for taxid in map(str,hcvspecies_descendents):
        if taxid in gtctaxids:
            evidence.add(ec('GTC',taxid))
	for sid in uckbtaxids:
	    if taxid in uckbtaxids[sid]:
                evidence.add(ec('UCKB',sid,taxid))

    if len(evidence) > 0:
        species[acc]['hcv'] = ('Direct',sorted(evidence))
    else:
        species[acc]['hcv'] = ('False',sorted(evidence))

for acc in sorted(species):
    print >>sys.stderr, "Pass 2:", acc
    for sp in species[acc]:
        if species[acc][sp][0] != 'Direct':
	    continue
        for anc in gnome.ancestors(acc):
            if anc not in species:
                continue
            category,evidence = species[anc][sp]
            if category == 'Direct':
		continue
            evidence = sorted(set(list(evidence) + [ec('sub',acc)]))
            species[anc][sp] = ('Subsumption',evidence)

for acc in sorted(species):
    print >>sys.stderr, "Pass 3:", acc
    for sp in species[acc]:
        if species[acc][sp][0] != 'Direct':
	    continue
        comp = gnome.get_composition(acc)
	if comp and comp in species:
            for desc in gnome.has_composition(comp):
	        if desc not in species:
		    continue
                category,evidence = species[desc][sp]
                if category in ('Direct','Subsumption'):
		    continue
                evidence = sorted(set(list(evidence) + [ec('compcl',acc,comp)]))
                species[desc][sp] = ('Composition',evidence)

for acc in sorted(species):
    m = w.get(acc)

    m.delete_annotations(type="Species",source="EdwardsLab")

    m.set_annotation(property="Human Evidence",
                     value=species[acc]['human'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Human",
                     value='true' if species[acc]['human'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['human'][0] != 'False':
        m.set_annotation(property="Human Category",
                         value=species[acc]['human'][0],
                         type="Species",
                         source="EdwardsLab")

    m.set_annotation(property="Mouse Evidence",
                     value=species[acc]['mouse'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Mouse",
                     value='true' if species[acc]['mouse'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['mouse'][0] != 'False':
        m.set_annotation(property="Mouse Category",
                         value=species[acc]['mouse'][0],
                         type="Species",
                         source="EdwardsLab")
    
    m.set_annotation(property="Rat Evidence",
                     value=species[acc]['rat'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Rat",
                     value='true' if species[acc]['rat'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['rat'][0] != 'False':
        m.set_annotation(property="Rat Category",
                         value=species[acc]['rat'][0],
                         type="Species",
                         source="EdwardsLab")

    m.set_annotation(property="HCV Evidence",
                     value=species[acc]['hcv'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="HCV",
                     value='true' if species[acc]['hcv'] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['hcv'][0] != 'False':
        m.set_annotation(property="HCV Category",
                         value=species[acc]['hcv'][0],
                         type="Species",
                         source="EdwardsLab")

    if w.put(m):
        print acc

