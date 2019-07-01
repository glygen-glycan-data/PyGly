#!/bin/env python27
import sys

from collections import defaultdict

import findpygly

# use the current raw GNOme subsumption dump
# Back-fill with GlyTouCan if needed.

from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan(usecache=True)

compdict = dict()
topodict = dict()
bcompdict = dict()
outedges = defaultdict(set)
inedges = defaultdict(set)
sec = None
for l in open(sys.argv[1]):
    if l.startswith('#'):
        sec = l.split()[1]
        continue
    if sec == "NODES":
        sl = l.split()
        acc = sl[0]
        if sl[3] != "None":
            topodict[acc] = (sl[3].strip("*"),sl[3].endswith("*"))
        else:
            topodict[acc] = (None,False)
        if sl[4] != "None":
            compdict[acc] = (sl[4].strip("*"),sl[4].endswith("*"))
        else:
            compdict[acc] = (None,False)
        if sl[5] != "None":
            bcompdict[acc] = (sl[5].strip("*"),sl[5].endswith("*"))
        else:
            bcompdict[acc] = (None,False) 
        continue                                                                                                  
    if sec == "EDGES":
	sl = l.split()
        f = sl[0].strip(':')
	for t in sl[1:]:
	    outedges[f].add(t)
	    inedges[t].add(f)
	continue

hascompdict = defaultdict(set)
for acc,c in compdict.items():
    hascompdict[c].add(acc)

dummy = sys.argv.pop(1)

def getanc(acc):
    anc = set()
    for p in inedges[acc]:
	anc.add(p)
	anc.update(getanc(p))
    return anc

def getdesc(acc):
    desc = set()
    for c in outedges[acc]:
	desc.add(c)
	desc.update(getdesc(c))
    return desc
 
def gettopo(acc):
    if acc in topodict:
	return topodict.get(acc)[0]
    return gtc.gettopo(acc)

def getcomp(acc):
    if acc in compdict:
	return compdict.get(acc)[0]
    return gtc.getcomp(acc)
    
def getbcomp(acc):
    if acc in bcompdict:
	return bcompdict.get(acc)[0]
    return gtc.getbasecomp(acc)

def hascomp(acc):
    if acc in hascompdict:
	return hascompdict[acc]
    return gtc.hascomp(acc)

accs = set()
for f in sys.argv[1:]:
    accs.update(open(f).read().split())

newaccs = set()

for acc in accs:
    newaccs.update(getanc(acc))
    topo = gettopo(acc)
    comp = getcomp(acc)
    if comp:
        newaccs.update(getdesc(comp))
    elif topo:
        newaccs.update(getdesc(topo))
newaccs = (newaccs-accs)

for acc in sorted(newaccs):
    print acc

