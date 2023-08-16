#!/bin/env python2

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
dir	        %s:%s TaxID %s
dirns		%s TaxID %s
noXylAlt	No N-linked Xyl and no Alt
noXylAltNeuGc	No N-linked Xyl and no NeuGc, Alt
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

species2taxid = dict()
species2taxid['human'] = set(map(int,"""
    9606 63221 741158
    """.split()))
species2taxid['mouse'] = set(map(int,"""
    10090 1385377 35531 179238 947985 1266728 10091 10092 80274 57486
    477816 39442 1879032 477815 46456 116058 1643390
    """.split()))
species2taxid['dros'] = set(map(int,"""
    7227
""".split()))
species2taxid['rat'] = set(map(int,"""
    10116 10114 2651915 2651916 323373 1258734 1258735 397343 10017 10018
    10019 10020 94247 94248 108146 323374 323375 323376 323377 323378
    323379 108145 323382 2338489 2338490 2338491 2338492 2338493 2338494
    2338495 2338496 2338497 2338498 1258742 50722 1258736 1258741 104298
    105255 147438 147439 108144 1258737 1258738 1258739 549621 549622
    549623 549624 549625 549626 284671
    """.split()))
species2taxid['hcv'] = set(map(int,"""
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
species2taxid['sars'] = set(map(int,"""
	694009 2697049 284672 305408 511433 1415851 305409 292360 305405
	227860 227861 1283332 698398 241183 1699360 1699361 293933 293934
	293935 293936 305416 698419 228404 285946 228406 228407 228409
	511430 228415 1503296 1503299 421444 1503301 1503302 237639 237640
	385685 333387 235088 249070 285267 233044 285947 231513 231514
	231515 231516 231517 231518 231519 231520 231521 231522 264379
	229992 229993 741997 741998 741999 742000 742001 273522 742003
	742004 742005 742006 267383 267384 267385 230522 267387 267388
	267389 267390 267391 267392 267393 267394 267395 267396 267397
	267398 260743 239241 239242 239243 262338 239758 227984 240273
	385683 385684 321149 385686 385687 388737 305402 1415834 227998
	349343 349344 235173 511431 1415852 338605 338606 291613 511432
	264370 264371 264372 264373 264374 264375 264376 264377 264378
	391355 264380 264381 264382 264383 264384 264385 253634 264387
	264388 304858 243925 353145 255194 312027 267386 321147 249062
	249063 249064 249065 249066 249067 249068 249069 252654 249071
	249072 249073 249074 249075 249076 249077 249078 249079 249080
	249081 249082 249083 249084 249085 249086 249087 249088 249089
	305410 305411 305412 305413 305414 305415 264386 305417 305418
	389166 389167 264988 264989 264990 296837 267399 1431340 742002
	305407 633137 633140 627442 1283333 299335 2042697 2042698 242743
	305404 229206 1487703 2697049 347536 247149 442736 258963 281976
	258964 255729 258965 1503300 1508227 258966 296838 296839 296840
	296841 258445 233872 347537 235410 235411 235412 235413 235414
	258967 258968 258969 258970 256753 258972 627440 344702 258971
	265124 240549 240550 240551 240552 240553 230471 255730 349342
	265125 305403 511429 260550 293319 293320 293321 258507 258508
	266147 285945 266148 1503303 241628 241629 248485 253435 305401
	260069 285948 228330 285949 305406 722424 253433 253434 228607
    """.split()))

# Saccharomyces cerevisiae 
species2taxid['yeast'] = set(map(int,"""
        4932 1294336 1294337 1294338 1294339 1294340 1294341 1294342
        1294343 1294344 1294345 1294346 1294347 1294348 1294349 1294350
        1138861 1294352 1294353 1294354 1294355 1294356 1294357 1294358
        1216345 927256 927257 927258 1294363 1294364 1294365 1294366
        1294367 1294368 1294369 1294370 947035 1294372 1294373 1294374
        1294375 1294376 1294377 1294378 1294379 1294380 1294381 1294382
        929629 1294384 1294385 1294386 1294387 1294388 1348153 538975
        1158204 1158205 2305205 947040 1337540 1266529 2305207 1337541
        468558 1337443 307796 502869 545124 1337538 1294351 1337436
        1337437 1337438 1337439 643680 1337441 1337442 1177187 1337444
        1337445 1337446 1337447 1337448 1337449 1337450 1337451 1337452
        1337453 1337454 1337455 1337456 1337457 1337458 1337459 1337460
        1337461 1293430 1337463 1337464 1337465 1337466 1337467 1337468
        1337469 1337470 1337471 1337472 1337473 1337474 1337475 1337476
        1337477 1337478 1337479 721032 1337440 1337482 1294359 1337484
        1337485 1337486 580239 580240 1294360 1337490 1337491 1337492
        1290389 1218710 1294361 1337496 464025 1337498 1337499 1337500
        1294362 1337502 1337503 1337504 1337505 1337506 1337507 1337508
        1337509 1337510 1337511 1337512 1337513 1337514 1292971 1337516
        889517 1337518 1337519 1337520 1337521 1337522 1337523 1337524
        1337525 1337526 1337527 1337528 1337529 1337530 1337531 559292
        1337533 1337534 1337535 1337536 764097 764098 764099 764100 764101
        764102 1337543 1337544 1352823 1337549 1337552 1337553 1337554
        1294371 1337556 1337557 1337558 1337559 1337560 1337561 1337562
        1337563 1337494 1337501 1337542 1337555 1182966 1182967 1182968
        1434269 1390931 1331972 471510 614664 614665 2305204 1204498
        1294383 1337462 466209 1337649 1337650 1337644 1337645 1337646
        1337647 2305203 929585 929586 471859 1337652 471861 1337654
        1337481 1352824 1337480 1149757 1337653 1196866 1337483 1337651
        658763 1419746 285006 717647 1390929 1390930 1097555 1390932
        1220494 1337487 1095001 1216859 947036 947037 947038 947039 538976
        947041 947042 947043 947044 947045 947046 1337489 1234807 929587
        1337493 765312 462209 462210 1330326 1418121 1296266 1337495
        41870 1294317 1337517 1337497 2305201 1416897 1144731 1382555
        1294321 1162671 1162672 1162673 1162674 1337537 1337532 1337655
        1391042 1337515 1337488 1247190 1201112 1337539 1227742 1294303
        1294304 1294305 1294306 1294307 1294308 1294309 1294310 1294311
        1294312 1294313 1294314 1294315 1294316 1087981 1294318 1294319
        1294320 574961 1294322 1294323 1294324 1294325 1294326 1294327
        1294328 1294329 1294330 1294331 1294332 1294333 1294334 1294335
""".split()))

# Dictyostelium discoideum
species2taxid['slimemold'] = set(map(int,"""
        44689 352472 366501 1592886
""".split()))

# Pig
species2taxid['pig'] = set(map(int,"""
9823 2485929 9825 1611878 1611879 1611880 1170810 291050 415978 310260 310261 490583 309913 309914 375579 375578
""".split()))

species = defaultdict(dict)

for m in iterglycan():
    acc = m.get('accession')
    print >>sys.stderr, "Pass 1:", acc

    taxidanns = defaultdict(set)
    for an in m.annotations(property="Taxonomy",type="Taxonomy"):
	source = an.get('source')
	sid = an.get('sourceid')
	for taxid in map(int,an.get('value')):
	    taxidanns[taxid].add((source,sid))

    hasmono = defaultdict(lambda: False)
    if m.has_annotations(property="HasMonosaccharide"):
        for mono in m.get_annotation_values(property="HasMonosaccharide"):
	    hasmono[mono] = True

    try:
        glycantype = m.get_annotation_value(property="GlycanType",source='EdwardsLab', type='Classification')
    except LookupError:
	glycantype = None
    try:
        glycansubtype = m.get_annotation_value(property="GlycanSubtype",source='EdwardsLab', type='Classification')
    except LookupError:
	glycansubtype = None

    for sp in species2taxid:
        evidence = set()
        for taxid in species2taxid[sp]:
	    for source,sid in taxidanns[taxid]:
		if sid != None:
		    evidence.add(ec('dir',source,sid,taxid))
		else:
		    evidence.add(ec('dirns',source,taxid))
        direct = False
	if sp == 'human':
            if len(evidence) > 0 and (glycantype != "N-linked" or not hasmono['Xyl']) and not hasmono['NeuGc'] and not hasmono['Alt']:
	        evidence.add(ec('noXylAltNeuGc'))
		direct = True
	elif sp in ('mouse','rat'):
    	    if len(evidence) > 0 and (glycantype != "N-linked" or not hasmono['Xyl']) and not hasmono['Alt']:
	        evidence.add(ec('noXylAlt'))
                direct = True
	else:
            if len(evidence) > 0:
                direct = True
        if direct:
           species[acc][sp] = ('Direct',sorted(evidence))
        else:
           species[acc][sp] = ('False',sorted(evidence))

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
                     value='true' if species[acc]['hcv'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['hcv'][0] != 'False':
        m.set_annotation(property="HCV Category",
                         value=species[acc]['hcv'][0],
                         type="Species",
                         source="EdwardsLab")

    m.set_annotation(property="SARS Evidence",
                     value=species[acc]['sars'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="SARS",
                     value='true' if species[acc]['sars'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['sars'][0] != 'False':
        m.set_annotation(property="SARS Category",
                         value=species[acc]['sars'][0],
                         type="Species",
                         source="EdwardsLab")

    m.set_annotation(property="FruitFly Evidence",
                     value=species[acc]['dros'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="FruitFly",
                     value='true' if species[acc]['dros'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['dros'][0] != 'False':
        m.set_annotation(property="FruitFly Category",
                         value=species[acc]['dros'][0],
                         type="Species",
                         source="EdwardsLab")

    m.set_annotation(property="Yeast Evidence",
                     value=species[acc]['yeast'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Yeast",
                     value='true' if species[acc]['yeast'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['yeast'][0] != 'False':
        m.set_annotation(property="Yeast Category",
                         value=species[acc]['yeast'][0],
                         type="Species",
                         source="EdwardsLab")

    m.set_annotation(property="Pig Evidence",
                     value=species[acc]['pig'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="Pig",
                     value='true' if species[acc]['pig'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['pig'][0] != 'False':
        m.set_annotation(property="Pig Category",
                         value=species[acc]['pig'][0],
                         type="Species",
                         source="EdwardsLab")

    m.set_annotation(property="SlimeMold Evidence",
                     value=species[acc]['slimemold'][1],
                     type="Species",
                     source="EdwardsLab")
    m.set_annotation(property="SlimeMold",
                     value='true' if species[acc]['slimemold'][0] in ('Direct','Subsumption') else 'false',
                     type="Species",
                     source="EdwardsLab")
    if species[acc]['slimemold'][0] != 'False':
        m.set_annotation(property="SlimeMold Category",
                         value=species[acc]['slimemold'][0],
                         type="Species",
                         source="EdwardsLab")

    if w.put(m):
        print acc

