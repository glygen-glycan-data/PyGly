import sys
from collections import defaultdict

# from getwiki import GlycanData
import findpygly
from pygly.GlycanFormatter import GlycoCTFormat
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCan

# w = GlycanData()
glycoctformat = GlycoCTFormat()

basecomp = {'x-HEX-x:x':'Hex',
            'x-HEX-x:x||(2d:1)n-acetyl':'HexNAc',
            'x-dgro-dgal-NON-x:x|1:a|2:keto|3:d||(5d:1)n-acetyl':'NeuAc',
            'x-dgro-dgal-NON-x:x|1:a|2:keto|3:d||(5d:1)n-glycolyl':'NeuGc',
            'x-HEX-x:x|6:d':'dHex', 
            'x-lgal-HEX-x:x|6:d':'Fuc',
	    'x-PEN-x:x':'Pent',
	    'phosphate':'P',
	    'sulfate':'S'}

# f = open('../data/basecomplist1.txt','w')

gtc = GlyTouCan()

def accessions():
    if len(sys.argv) == 2 and sys.argv[1] == "-":
	for acc in sys.stdin:
	    yield acc.strip()
    elif len(sys.argv) == 2 and sys.argv[1] == "*":
	for acc in gtc.allacc():
	    yield acc
    else:
	for acc in sys.argv[1:]:
	    yield acc

for acc in accessions():
    glycan = gtc.getGlycan(acc,format='wurcs')
    if not glycan:
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

    byonic_good = False
    if node_count == sum(map(lambda k: comp.get(k,0),("HexNAc","Hex","Fuc","dHex","NeuAc","NeuGc"))):
	byonic_good = True

    uckb_good = False
    if node_count == sum(map(lambda k: comp.get(k,0),("HexNAc","Hex","dHex","NeuAc","NeuGc","P","S","Pent"))):
	uckb_good = True

    short_good = False
    if node_count == sum(map(lambda k: comp.get(k,0),("HexNAc","Hex","Fuc","NeuAc"))):
	short_good = True

    if not byonic_good and not uckb_good and not short_good:
	continue

    byonic_string = ''
    similar_string = ''
    full_string = ''
    short_string = ''

    if comp['Fuc'] !=0:
        for k in ['HexNAc','Hex','Fuc','NeuAc','NeuGc']:
            if comp[k] != 0:
                byonic_string += k + '(' + str(comp[k]) + ')'
        byonic_line = acc + '\t' + byonic_string + "\tBYONIC\n"
        if byonic_string:
            print byonic_line.strip()
            # f.write(byonic_line)
        comp['H'] = comp['Hex']
        comp['N'] = comp['HexNAc']
        comp['F'] = comp['Fuc']
        comp['Sia'] = comp['NeuAc'] + comp['NeuGc']
        for k in ['H','N','F','Sia']:
            if comp[k] != 0:
                if comp[k] == 1:
                    if k == 'Sia':
                        short_string += 'S'
                    else:
                        short_string += k
                else:
                    if k == 'Sia':
                        short_string += 'S' + str(comp['Sia'])
                    else:
                        short_string += k + str(comp[k])
        short_line = acc + '\t' + short_string + '\tSHORTCOMP\n'
        if short_good:
            print short_line.strip()
            # f.write(short_line)
    else:
        for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','P','S']:
            if comp[k] != 0:
                byonic_string += k + '(' + str(comp[k]) + ')'
                similar_string += k + str(comp[k])
        byonic_line = acc + '\t' + byonic_string + '\tBYONIC\n'
        similar_line = acc + '\t' + similar_string + '\tSHORTUCKB\n'
        if byonic_good:
            print byonic_line.strip()
            # f.write(byonic_line)
        if uckb_good:
            print similar_line.strip()
            # f.write(similar_line)
        for k in ['KDN','HexA']:
            comp[k] = 0
        for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P','KDN','HexA']:
            full_string += k + str(comp[k])
        full_line = acc + '\t' + full_string + '\tUCKBCOMP\n'
        if uckb_good:
            print full_line.strip()
            # f.write(full_line)
        if comp['dHex'] == 0:
            comp['H'] = comp['Hex']
            comp['N'] = comp['HexNAc']
            comp['Sia'] = comp['NeuAc'] + comp['NeuGc']
            for k in ['H','N','Sia']:
                if comp[k] != 0:
                    if comp[k] == 1:
                        if k == 'Sia':
                            short_string += 'S'
                        else:
                            short_string += k
                    else:
                        if k == 'Sia':
                            short_string += 'S' + str(comp['Sia'])
                        else:
                            short_string += k + str(comp[k])
            short_line = acc + '\t' + short_string + '\tSHORTCOMP\n'
            if short_good:
                print short_line.strip()
                # f.write(short_line)
# f.close()

