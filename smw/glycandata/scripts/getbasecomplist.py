import sys, re
from collections import defaultdict

# from getwiki import GlycanData
import findpygly
from pygly.Glycan import Glycan
from pygly.GlycanFormatter import GlycoCTFormat
from pygly.GlycanResource import GlyTouCan, GlyCosmos

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
# gtc = GlyTouCan(usecache=False)
gco = GlyCosmos(usecache=False)
replace = gco.replace()
allgco = gco.allaccessions()

gtc = None

def accessions():
    global gtc
    if len(sys.argv) == 2 and sys.argv[1] == "-":
        gtc = GlyTouCan(usecache=False)
        for acc in sys.stdin:
            yield acc.strip()
    elif len(sys.argv) == 2 and sys.argv[1] == "*":
        gtc = GlyTouCan(usecache=False)
        for acc in gtc.allaccessions():
            if acc in replace:
                continue
            yield acc
    else:
        gtc = GlyTouCan(prefetch=False)
        for acc in sys.argv[1:]:
            yield acc

allskel = set()
names = defaultdict(set)

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
        if link == "":
            continue
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

    special_node_count_one = False
    if glycan.has_root() and node_count == 1:
        r = glycan.root()
        floating_subst = []
        for sl in r.substituent_links():
            try:
                glycoctsym = glycoctformat.mtoStr(sl.child())
            except KeyError:
                continue
            if sl.parent_pos() != None:
                continue
            if glycoctsym in basecomp:
                sub = sl.child()
                sub.set_connected(False)
                r.remove_substituent_link(sl)
                floating_subst.append(sub)
        glycan.set_root(None)
        glycan.set_undetermined([r]+floating_subst)
        node_count = len(floating_subst)+1
        special_node_count_one = True

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
        if not special_node_count_one:
            print("Unexpected skeleton codes:",", ".join(newskel - expskel),file=sys.stderr)
            sys.exit(1)
        else:
            newelts = newskel.pop().split('_')
            newskel = []
            for e in newelts:
                if e not in ('?*OPO/3O/3=O','?*OSO/3=O/3=O'):
                    newskel.append(e)
            newskel = "_".join(newskel)
            if newskel not in expskel:
                print("Unexpected skeleton codes:",newskel,file=sys.stderr)
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
                             ("HexNAc","Hex","Fuc","NeuAc","NeuGc"))):
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
            names[("BYONIC",byonic_string)].add(acc)

        if short_good:

            short_string = ''
            comp['H'] = comp['Hex']
            comp['N'] = comp['HexNAc']
            comp['F'] = comp['Fuc']
            comp['Sia'] = comp['NeuAc']
            comp['G'] = comp['NeuGc']
            for k in ['H','N','F','Sia','G']:
                if comp[k] != 0:
                    short_string += k[0]
                    if comp[k] > 1:
                        short_string += str(comp[k])
            names[("SHORTCOMP",short_string)].add(acc)

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
            names[("BYONIC",byonic_string)].add(acc)

        uckb_string = ''
        shortuckb_string = ''
        for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P','KDN','HexA']:
            uckb_string += k + str(comp[k])
            if comp[k] > 0:
                shortuckb_string += k + str(comp[k])
        if uckb_good:
            names[("UCKBCOMP",uckb_string)].add(acc)
            names[("SHORTUCKB",shortuckb_string)].add(acc)

        if short_good and comp['dHex'] == 0:

            short_string = ''
            comp['H'] = comp['Hex']
            comp['N'] = comp['HexNAc']
            comp['Sia'] = comp['NeuAc']
            comp['G'] = comp['NeuGc']
            for k in ['H','N','Sia','G']:
                if comp[k] > 0:
                    short_string += k[0]
                if comp[k] > 1:
                    short_string += str(comp[k])
            names[("SHORTCOMP",short_string)].add(acc)

for thetype,thestr in names:
    if len(names[(thetype,thestr)]) == 1:
        print(names[(thetype,thestr)].pop(),thestr,thetype)
    else:
        gcoaccs = filter(lambda acc: acc in allgco, names[(thetype,thestr)])
        if len(gcoaccs) == 1:
            print(gcoaccs[0],thestr,thetype)
        elif len(gcoaccs) > 1:
            print("WARNING: Can't pick  correct accession for %s:%s, possibilities: %s."%(thetype,thestr,", ".join(sorted(names[(thetype,thestr)]))),file=sys.stderr)
            print(sorted(gcoaccs)[0],thestr,thetype)
        else:
            print("WARNING: Can't pick  correct accession for %s:%s, possibilities: %s."%(thetype,thestr,", ".join(sorted(names[(thetype,thestr)]))),file=sys.stderr)
            print(sorted(names[(thetype,thestr)])[0],thestr,thetype)
