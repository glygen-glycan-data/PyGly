#!/bin/env python3
import sys, os.path
import findpygly
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCanNoPrefetch, GlyTouCan
from pygly.Monosaccharide import constantString, Substituent
from pygly.WURCS20MonoFormatter import WURCS20MonoFormat
from pygly.alignment import MonosaccharideEqual, SubstituentEqual, LinkageEqual
        

wp = WURCS20MonoFormat()
mx = wp.get("ha122h-2b_2-5")
print(mx)

monoeq = MonosaccharideEqual(substcmp=SubstituentEqual(),sublinkcmp=LinkageEqual())

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
    for l in g.all_links():
        if l.parent_pos() == set([1]) and l.child_pos() == set([1]) and (l.child() == g.root() or l.parent() == g.root()):
            found = True
            break
    if not found:
        for m in g.all_nodes(subst=False):
            if m == mx:
                found = True
                break
    if found:
        print(acc)
