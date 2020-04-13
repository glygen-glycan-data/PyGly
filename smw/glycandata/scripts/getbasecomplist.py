import sys, re
from collections import defaultdict

# from getwiki import GlycanData
import findpygly
from pygly.GlycanFormatter import GlycoCTFormat
from pygly.GlycanResource import GlyTouCan

# w = GlycanData()
glycoctformat = GlycoCTFormat()

basecomp = {'x-HEX-x:x':'Hex',
            'x-HEX-x:x||(2d:1)n-acetyl':'HexNAc',
            'x-dgro-dgal-NON-x:x|1:a|2:keto|3:d||(5d:1)n-acetyl':'NeuAc',
            'x-dgro-dgal-NON-x:x|1:a|2:keto|3:d||(5d:1)n-glycolyl':'NeuGc',
            'x-HEX-x:x|6:d':'dHex', 
            'x-lgal-HEX-x:x|6:d':'Fuc',
	    'x-PEN-x:x':'Pent',
	    'x-dgro-dgal-NON-x:x|1:a|2:keto|3:d':'KDN',
	    'x-HEX-x:x|6:a':'HexA',
	    'phosphate':'P',
	    'sulfate':'S'}

badskel = set("""
axxxxh-1x
""".split())

# NeuAc, NeuGc, KDN, Fuc, Hex, HexNAc, dHex, HexA, Pent
expskel = set("""
  AUd21122h_5*NCC/3=O
  AUd21122h_5*NCCO/3=O
  AUd21122h
  u1221m
  uxxxxh
  uxxxxh_2*NCC/3=O
  uxxxxm
  uxxxxA
  uxxxh
""".split())

#
# G45924NL is floating in-link substituent - IMO this is not a correct representation
# G99993XU ditto
# G94401PZ ditto
# G10633KT ditto
# G48302UE ditto
# 

badacc = set("""
""".split())

# f = open('../data/basecomplist1.txt','w')

# Make sure we get the latest version of everything
gtc = GlyTouCan(usecache=False)

def accessions():
    if len(sys.argv) == 2 and sys.argv[1] == "-":
	for acc in sys.stdin:
	    yield acc.strip()
    elif len(sys.argv) == 2 and sys.argv[1] == "*":
	for acc in gtc.allaccessions():
	    yield acc
    else:
	for acc in sys.argv[1:]:
	    yield acc

allskel = set()

seen = defaultdict(set)

for acc in sorted(accessions()):

    if acc in badacc:
	continue

    glycan = gtc.getGlycan(acc,format='wurcs')
    if not glycan:
        continue

    wurcs = gtc.getseq(acc,format='wurcs')
    allgood = True
    for skel in wurcs.split('/[', 1)[1].split(']/')[0].split(']['):
	if skel in badskel:
	    allgood = False
	    break

    if not allgood:
	continue

    # hack! detect bad subst in link composition!
    badlink = 0
    for link in wurcs.split(']/')[1].split('/',1)[1].split('_'):
        # a?|b?}-{a?|b?
	if re.search(r'^([a-zA-Z]\?(\|[a-zA-Z]\?)+)\}-\{\1$',link):
	    continue
	# a?|b?}*OSO/3=O/3=O
	if re.search(r'^([a-zA-Z]\?(\|[a-zA-Z]\?)+)\}\*',link):
	    continue
	badlink += 1

    if badlink > 0:
	continue

    node_count = 0
    for m in glycan.all_nodes(undet_subst=True):
        node_count += 1

    if node_count == 0:
	continue

    if glycan.has_root() and node_count != 1:
        continue

    comp = defaultdict(int)

    for m in glycan.all_nodes(undet_subst=True):
        try:
            glycoctsym = glycoctformat.mtoStr(m)
        except KeyError:
            break
        if glycoctsym not in basecomp:
            break
        comp[basecomp[glycoctsym]] += 1

    if node_count != sum(comp.values()):
	continue

    if comp['Fuc'] > 0 and comp['dHex'] > 0:
	continue

    newskel = set(wurcs.split('/[', 1)[1].split(']/')[0].split(']['))
    if len(newskel - expskel) > 0:
	print >>sys.stderr, "Unexpected skeleton codes:",", ".join(newskel - expskel)
	sys.exit(1)

    byonic_good = False
    if node_count == sum(map(lambda k: comp.get(k,0),
                             ("HexNAc","Hex","Fuc","dHex","NeuAc","NeuGc","Pent","S","P"))):
	byonic_good = True

    uckb_good = False
    if node_count == sum(map(lambda k: comp.get(k,0),
                             ("HexNAc","Hex","dHex","NeuAc","NeuGc","Pent","S","P","KDN","HexA"))):
	uckb_good = True

    short_good = False
    if node_count == sum(map(lambda k: comp.get(k,0),
                             ("HexNAc","Hex","Fuc","NeuAc"))):
	short_good = True

    if not byonic_good and not uckb_good and not short_good:
	continue

    if comp['Fuc'] > 0:

        if byonic_good:

            byonic_string = ''
            for k in ['HexNAc','Hex','Fuc','NeuAc','NeuGc','Pent','S','P']:
                if comp[k] != 0:
		    if k == 'P':
                        byonic_string += 'Phospho' + '(' + str(comp[k]) + ')'
		    elif k == 'S':
                        byonic_string += 'Sulpho' + '(' + str(comp[k]) + ')'
		    else:
                        byonic_string += k + '(' + str(comp[k]) + ')'
	    if ("BYONIC",byonic_string) not in seen:
	        print "\t".join([acc,byonic_string,"BYONIC"])
	    seen[("BYONIC",byonic_string)].add(acc)

	if short_good:

	    short_string = ''
            comp['H'] = comp['Hex']
            comp['N'] = comp['HexNAc']
            comp['F'] = comp['Fuc']
            comp['Sia'] = comp['NeuAc']
            for k in ['H','N','F','Sia']:
                if comp[k] != 0:
		    short_string += k[0]
		    if comp[k] > 1:
		        short_string += str(comp[k])
	    if ("SHORTCOMP",short_string) not in seen:
	        print "\t".join([acc,short_string,"SHORTCOMP"])
	    seen[("SHORTCOMP",short_string)].add(acc)

    else:

	byonic_string = ''
        for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P']:
            if comp[k] > 0:
		if k == 'P':
                    byonic_string += 'Phospho' + '(' + str(comp[k]) + ')'
		elif k == 'S':
                    byonic_string += 'Sulpho' + '(' + str(comp[k]) + ')'
		else:
                    byonic_string += k + '(' + str(comp[k]) + ')'
        if byonic_good:
	    if ("BYONIC",byonic_string) not in seen:
	        print "\t".join([acc,byonic_string,"BYONIC"])
	    seen[("BYONIC",byonic_string)].add(acc)

	uckb_string = ''
	shortuckb_string = ''
        for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P','KDN','HexA']:
            uckb_string += k + str(comp[k])
	    if comp[k] > 0:
                shortuckb_string += k + str(comp[k])
        if uckb_good:
	    if ("UCKBCOMP",uckb_string) not in seen:
	        print "\t".join([acc,uckb_string,"UCKBCOMP"])
	    seen[("UCKBCOMP",uckb_string)].add(acc)
	    if ("SHORTUCKB",shortuckb_string) not in seen:
	        print "\t".join([acc,shortuckb_string,"SHORTUCKB"])
	    seen[("SHORTUCKB",shortuckb_string)].add(acc)

        if short_good and comp['dHex'] == 0:

	    short_string = ''
            comp['H'] = comp['Hex']
            comp['N'] = comp['HexNAc']
            comp['Sia'] = comp['NeuAc']
            for k in ['H','N','Sia']:
                if comp[k] > 0:
	    	    short_string += k[0]
		    if comp[k] > 1:
		        short_string += str(comp[k])
	    if ("SHORTCOMP",short_string) not in seen:
                print "\t".join([acc,short_string,"SHORTCOMP"])
	    seen[("SHORTCOMP",short_string)].add(acc)

for thetype,thestr in seen:
    if len(seen[(thetype,thestr)]) > 1:
	print >>sys.stderr, "WARNING: Repeat %s:%s for %s."%(thetype,thestr,", ".join(sorted(seen[(thetype,thestr)])))
