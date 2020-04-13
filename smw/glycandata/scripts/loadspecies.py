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
noXyl		No N-linked Xyl
noXylNeuGc	No N-linked Xyl and no NeuGc
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
        gtctaxids = set(map(int,m.get_annotation_values(property="Taxonomy",source="GlyTouCan",type="Taxonomy")))
    except LookupError:
	gtctaxids = set()

    uckbtaxids = defaultdict(set)
    for an in m.annotations(property="Taxonomy",source="UniCarbKB",type="Taxonomy"):
	sid = an.get('sourceid')
	if sid:
	    uckbtaxids[sid].update(map(int,an.get('value')))

    comp = defaultdict(int)
    for an in m.annotations(type="MonosaccharideCount"):
	prop = an.get('property')
	if prop.endswith('Count'):
	    prop = prop[:-5]
	if prop == "Monosaccharide":
	    prop = "Count"
	comp[prop] = int(an.get('value'))

    try:
        glycantype = m.get_annotation_value(property="GlycanType",source='EdwardsLab', type='Classification')
    except LookupError:
	glycantype = None
    try:
        glycansubtype = m.get_annotation_value(property="GlycanSubtype",source='EdwardsLab', type='Classification')
    except LookupError:
	glycansubtype = None

    human_taxids = set(map(int,"9606 63221 741158".split()))

    evidence = set()
    for taxid in human_taxids:
        if taxid in gtctaxids:
            evidence.add(ec('GTC',taxid))
	for sid in uckbtaxids:
	    if taxid in uckbtaxids[sid]:
                evidence.add(ec('UCKB',sid,taxid))
    if len(evidence) > 0 and (glycantype != "N-linked" or comp['Xyl'] == 0) and comp['NeuGc'] == 0:
	evidence.add(ec('noXylNeuGc'))

    if len(evidence) > 1 and ec('noXylNeuGc') in evidence:
        species[acc]['human'] = ('Direct',sorted(evidence))
    else:
        species[acc]['human'] = ('False',sorted(evidence))

    mouse_taxids = set(map(int,"""
    10090 1385377 35531 179238 947985 1266728 10091 10092 80274 57486 477816 39442 1879032 477815 46456 116058 1643390
    """.split()))

    evidence = set()
    for taxid in mouse_taxids:
        if taxid in gtctaxids:
            evidence.add(ec('GTC',taxid))
	for sid in uckbtaxids:
	    if taxid in uckbtaxids[sid]:
                evidence.add(ec('UCKB',sid,taxid))
    if len(evidence) > 0 and (glycantype != "N-linked" or comp['Xyl'] == 0):
	evidence.add(ec('noXyl'))

    if len(evidence) > 1 and ec('noXyl') in evidence:
        species[acc]['mouse'] = ('Direct',sorted(evidence))
    else:
        species[acc]['mouse'] = ('False',sorted(evidence))

    rat_taxids = set(map(int,"""
    10116 10114 2651915 2651916 323373 1258734 1258735 397343 10017 10018
    10019 10020 94247 94248 108146 323374 323375 323376 323377 323378
    323379 108145 323382 2338489 2338490 2338491 2338492 2338493 2338494
    2338495 2338496 2338497 2338498 1258742 50722 1258736 1258741 104298
    105255 147438 147439 108144 1258737 1258738 1258739 549621 549622
    549623 549624 549625 549626 284671
    """.split()))

    evidence = set()
    for taxid in rat_taxids:
        if taxid in gtctaxids:
            evidence.add(ec('GTC',taxid))
	for sid in uckbtaxids:
	    if taxid in uckbtaxids[sid]:
                evidence.add(ec('UCKB',sid,taxid))
    if len(evidence) > 0 and (glycantype != "N-linked" or comp['Xyl'] == 0):
	evidence.add(ec('noXyl'))

    if len(evidence) > 1 and ec('noXyl') in evidence:
        species[acc]['rat'] = ('Direct',sorted(evidence))
    else:
        species[acc]['rat'] = ('False',sorted(evidence))
    
    hcv_taxids = set(map(int,"""
	11103 945067 578306 945063 1173523 1544726 484894 945073 356386
	356387 356388 356390 356391 2340905 2340906 2340907 1405107
	1972279 1972280 356410 356411 356412 356413 356414 356415 356416
	356417 356418 356419 356420 356421 356422 413255 413256 413257
	356426 356427 1006434 758864 1006435 1006436 761189 438880 438881
	357986 357987 357988 357989 357990 357991 761191 761967 356465
	356466 356467 356468 356469 1006441 1006432 1006442 759939 2491016
	1006444 378506 378507 1006439 595609 595610 595611 595612 1544902
	1406421 329389 693426 693427 693428 693429 1006437 1193076 1544901
	42182 487624 1006438 1053139 36390 745709 1006443 1094895 760561
	1094898 1094899 1094900 1094901 1406420 761190 1208062 1208063
	63746 1406424 1406422 569610 1406423 679182 578319 356113 356114
	356115 356116 357985 1406426 467354 761192 1406425 1406427 42791
	42792 1006440 945065 128819 128820 128821 668553 1259832 668554
	668555 420174 40271 1006431 11104 11105 11106 11107 11108 11109
	11110 11111 11112 11113 11114 11115 11116 11117 1006433 41856
	467337 467338 467339 31642 31643 31644 31645 31646 31647 31648
	31649 31650 31651 945060 31653 31654 31655 40360 40361 945066
	356423 945068 945069 945070 945071 945072 356424 356425 44021
	945057 44023 945058 33745 33746 945059 1053140 1053141 1053142
	1053143 1053144 1053145 1053146 1053147 945061 333284 945062
	357355 945064 44019 44020 421877 44022 421879 578303
    """.split()))

    evidence = set()
    for taxid in hcv_taxids:
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

