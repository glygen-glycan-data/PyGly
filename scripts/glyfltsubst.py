#!/bin/env python3
import sys, os.path
import findpygly
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCanNoPrefetch, GlyTouCan
from pygly.Monosaccharide import constantString, Substituent

if len(sys.argv) > 1:
    gtc = GlyTouCanNoPrefetch()
    iterable = sys.argv[1:]
else:
    gtc = GlyTouCan()
    iterable = gtc.allaccessions()
for acc in iterable:
    try:
        g = gtc.getGlycan(acc)
    except:
        continue
    if not g:
        continue
    if not g.has_root():
        continue
    found = False
    for n in g.all_nodes(subst=False,undet_subst=True):
        if not n.is_monosaccharide():
            print(acc,n.name(),constantString(Substituent,n.name()))
            found = True
            break
    if not found:
        pass #print(acc)
