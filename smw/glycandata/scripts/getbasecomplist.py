import sys
from collections import defaultdict

from getwiki import GlycanData, Glycan
from pygly.GlycanFormatter import GlycoCTFormat

w = GlycanData()
glycoctformat = GlycoCTFormat()

basecomp = {'x-HEX-x:x':'Hex',
            'x-HEX-x:x||(2d:1)n-acetyl':'HexNAc',
            'x-dgro-dgal-NON-x:x|1:a|2:keto|3:d||(5d:1)n-acetyl':'NeuAc',
            'x-dgro-dgal-NON-x:x|1:a|2:keto|3:d||(5d:1)n-glycolyl':'NeuGc',
            'x-HEX-x:x|6:d':'dHex', 
            'x-lgal-HEX-x:x|6:d':'Fuc'}

f = open('../data/basecomplist1.txt','w')

for g in w.iterglycan():
    acc = g.get('accession')
    glycan = g.getGlycan()
    if not glycan:
        continue
#    if glycan.has_root():
#        continue
    comp = defaultdict(int)

    node_count = 0
    for m in glycan.all_nodes():
        node_count += 1

    for m in glycan.all_nodes():
        try:
            glycoctsym = glycoctformat.mtoStr(m)
        except KeyError:
            continue
        if glycan.has_root() and node_count != 1:
            continue
        if glycoctsym not in basecomp:
            comp.clear()
            break
        else:
            comp[basecomp[glycoctsym]] += 1
            for k in ['Pent','S','P','KDN','HexA']:
                comp[k] = 0
    byonic_string = ''
    similar_string = ''
    full_string = ''
    short_string = ''

    if comp['Fuc'] !=0:
        for k in ['HexNAc','Hex','Fuc','NeuAc','NeuGc']:
            if comp[k] != 0:
                byonic_string += k + '(' + str(comp[k]) + ')'
        byonic_line = acc + '\t' + byonic_string + '\n'
        if byonic_string:
            print byonic_line
            f.write(byonic_line)
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
        short_line = acc + '\t' + short_string + '\n'
        if short_string:
            print short_line
            f.write(short_line)
    else:
        for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc']:
            if comp[k] != 0:
                byonic_string += k + '(' + str(comp[k]) + ')'
                similar_string += k + str(comp[k])
        byonic_line = acc + '\t' + byonic_string + '\n'
        similar_line = acc + '\t' + similar_string + '\n'
        if byonic_string:
            print byonic_line
            f.write(byonic_line)
        if similar_string:
            print similar_line
            f.write(similar_line)
        for k in ['Pent','S','P','KDN','HexA']:
                comp[k] = 0
        for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P','KDN','HexA']:
            full_string += k + str(comp[k])
            full_line = acc + '\t' + full_string + '\n'
        if byonic_string:
            print full_line
            f.write(full_line)
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
            short_line = acc + '\t' + short_string + '\n'
            if short_string:
                print short_line
                f.write(short_line)
f.close()

