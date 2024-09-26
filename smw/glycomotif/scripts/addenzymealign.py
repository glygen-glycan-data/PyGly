#!/bin/env python27

from getwiki import GlycoMotifWiki
import sys, re, csv, gzip, copy
from collections import defaultdict

w = GlycoMotifWiki()

alignmap = dict(core="Core",
                substr="Substructure",
                nonred="Nonreducing-End",
                whole="Whole-Glycan")

def alignmentiter(filename):
    k = None
    for r in csv.DictReader(gzip.open(filename),dialect='excel-tab'):
        mgtc = r['Motif']
        struct = r['Structure']
        altype = r['AlignmentType']
        alind = r['AlignmentIndex']
        if (mgtc,struct,altype,alind) != k:
            if k != None and len(rows) > 0:
                core = list(filter(lambda gtid: gtid == "NA" or gtid.startswith('OC'),[ r1['GlycoTreeID'] for r1 in rows ]))[0]
                for r1 in rows:
                    r1['GlycoTreeCore'] = core
                    yield r1
	    k = (mgtc,struct,altype,alind)
            rows = []
        rows.append(dict(r.items()))
    if len(rows) > 0:
        core = list(filter(lambda gtid: gtid == "NA" or gtid.startswith('OC'),[ r['GlycoTreeID'] for r in rows ]))[0]
        for r in rows:
            r['GlycoTreeCore'] = core
            yield r
        
allcores = set()
enzdata = defaultdict(dict)
laststruct = None
for r in alignmentiter(sys.argv[1]):
    if r['MotifResidue'] == '-':
        continue
    if r['Structure'] != laststruct:
        print(r['Structure'])
        laststruct = r['Structure']
    mgtc = r['Motif']
    align = alignmap[r['AlignmentType']]
    core = r['GlycoTreeCore']
    allcores.add(core)
    nlinked = (core == "NA")
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
    if ('human',core) not in enzdata[mgtc,align]:
        enzdata[mgtc,align]['human',core] = {}
    if ('mouse',core) not in enzdata[mgtc,align]:
        enzdata[mgtc,align]['mouse',core] = {}
    for e in he:
        if e not in enzdata[mgtc,align]['human',core]:
            enzdata[mgtc,align]['human',core][e] = he[e]
        else:
            enzdata[mgtc,align]['human',core][e].update(he[e])
    for e in me:
        if e not in enzdata[mgtc,align]['mouse',core]:
            enzdata[mgtc,align]['mouse',core][e] = me[e]
        else:
            enzdata[mgtc,align]['mouse',core][e].update(me[e])
    if ('sandbox',core) not in enzdata[mgtc,align]:
        enzdata[mgtc,align]['sandbox',core] = set()
    enzdata[mgtc,align]['sandbox',core].add(r['Structure'])

for m in w.itermotif(collection="GGM"):
    mgtc = m.get('glytoucan')
    align = m.get('alignment')[0]
    for species in ('human','mouse'):
      for core in allcores:
        if (species,core) not in enzdata[mgtc,align]:
            continue
        enz = enzdata[mgtc,align][species,core]
        key = "%s_%s_enzymes"%(species,core)
        if len(enz) > 0:
            m.set(key,[ "%s:%s"%(t[0],",".join(t[1])) for t in enz.items() ])
        else:
            m.delete(key)
    for core in allcores:
      if len(enzdata[mgtc,align].get(('sandbox',core),[])) > 0:
        m.set("sandbox_"+core, sorted(enzdata[mgtc,align]['sandbox',core])[:10])
      else:
        m.delete("sandbox_"+core)
    m.delete("humanenzyme")
    m.delete("mouseenzyme")
    m.delete("sandbox_nlinked")
    m.delete("sandbox_olinked")
    m.delete("sandbox")
    m.delete("enzyme")
    if w.put(m):
        print >>sys.stderr, m.get('id')
