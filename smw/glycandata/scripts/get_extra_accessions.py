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

hascompdict = defaultdict(set)
for acc,c in compdict.items():
    hascompdict[c].add(acc)

dummy = sys.argv.pop(1)
 
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
    topo = gettopo(acc)
    comp = getcomp(acc)
    bcomp = getbcomp(acc)
    if topo:
        newaccs.add(topo)
    if comp:
        for newacc in hascomp(comp):
            newaccs.add(newacc)    
        newaccs.add(comp)
    if bcomp:
        newaccs.add(bcomp)
newaccs = (newaccs-accs)

for acc in sorted(newaccs):
    print acc

