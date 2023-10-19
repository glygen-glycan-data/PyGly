#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys, re, csv, gzip
from collections import defaultdict

w = GlycoMotifWiki()

alignmap = dict(core="Core",
                substr="Substructure",
                nonred="Nonreducing-End",
                whole="Whole-Glycan")

enzdata = defaultdict(dict)
laststruct = None
for r in csv.DictReader(gzip.open(sys.argv[1]),dialect='excel-tab'):
    mgtc = r['Motif']
    if r['MotifResidue'] == '-':
        continue
    if r['Structure'] != laststruct:
        print(r['Structure'])
        laststruct = r['Structure']
    align = alignmap[r['AlignmentType']]
    nlinked = r['GlycoTreeID'].startswith('N')
    he = {}
    for e in r['HumanEnzyme'].split(','):
        if e == "-":
            continue
        if e not in he:
            he[e] = set([r['MotifResidue']])
        else:
            he[e].add(r['MotifResidue'])
    me = {}
    for e in r['MouseEnzyme'].split(','):
        if e == "-":
            continue
        if e not in he:
            me[e] = set([r['MotifResidue']])
        else:
            me[e].add(r['MotifResidue'])
    if 'human' not in enzdata[mgtc,align]:
        enzdata[mgtc,align]['human'] = {}
    if 'mouse' not in enzdata[mgtc,align]:
        enzdata[mgtc,align]['mouse'] = {}
    for e in he:
        if e not in enzdata[mgtc,align]['human']:
            enzdata[mgtc,align]['human'][e] = he[e]
        else:
            enzdata[mgtc,align]['human'][e].update(he[e])
    for e in me:
        if e not in enzdata[mgtc,align]['mouse']:
            enzdata[mgtc,align]['mouse'][e] = me[e]
        else:
            enzdata[mgtc,align]['mouse'][e].update(me[e])
    if 'sandbox_nlinked' not in enzdata[mgtc,align]:
        enzdata[mgtc,align]['sandbox_nlinked'] = set()
        enzdata[mgtc,align]['sandbox_olinked'] = set()
    if nlinked:
        enzdata[mgtc,align]['sandbox_nlinked'].add(r['Structure'])
    else:
        enzdata[mgtc,align]['sandbox_olinked'].add(r['Structure'])

for m in w.itermotif(collection="GGM"):
    mgtc = m.get('glytoucan')
    align = m.get('alignment')[0]
    he = enzdata[mgtc,align].get('human',{})
    if len(he) > 0:
        m.set("humanenzyme",[ "%s:%s"%(t[0],",".join(t[1])) for t in he.items() ])
    else:
        m.delete("humanenzyme")
    me = enzdata[mgtc,align].get('mouse',{})
    if len(me) > 0:
        m.set("mouseenzyme",[ "%s:%s"%(t[0],",".join(t[1])) for t in me.items() ])
    else:
        m.delete("mouseenzyme")
    if len(enzdata[mgtc,align].get('sandbox_nlinked',[])) > 0:
        m.set("sandbox_nlinked", sorted(enzdata[mgtc,align]['sandbox_nlinked'])[:10])
    else:
        m.delete("sandbox_nlinked")
    if len(enzdata[mgtc,align].get('sandbox_olinked',[])) > 0:
        m.set("sandbox_olinked", sorted(enzdata[mgtc,align]['sandbox_olinked'])[:10])
    else:
        m.delete("sandbox_olinked")
    m.delete("sandbox")
    m.delete("enzyme")
    if w.put(m):
        print >>sys.stderr, m.get('id')
