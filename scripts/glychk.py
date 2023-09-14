#!/bin/env python3
import sys, os.path
import findpygly
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCanNoPrefetch

if len(sys.argv) > 1:
    gtc = GlyTouCanNoPrefetch()
    iterable = sys.argv[1:]
else:
    gtc = GlyTouCanNoCache()
    iterable = gtc.allaccessions()
for acc in iterable:
    res = {"acc": acc}
    for fmt in ('wurcs','glycoct'):
        g = None
        try:
            g = gtc.getGlycan(acc,format=fmt)
        except:
            pass
        good = None
        badres = None
        msg = None
        if g:
            good=True
            for m in g.all_nodes():
                try:
                    m.has_valid_positions()
                except RuntimeError as e:
                    good = False
                    badres = m.id()
                    msg = e.args[0]
                    break
        res[fmt] = (good,badres if badres else "",msg if msg else "")
    print("\t".join(map(str,[res["acc"]] + list(res["wurcs"]) + list(res["glycoct"]))))
