#!/bin/env python27

import findpygly
from pygly.GlycanResource import GlyTouCan
import sys, re

gtc = GlyTouCan(retries=0)

symbol2wurcs_definition = """
NeuAc   AUd21122h_5*NCC/3=O   3
NeuGc   AUd21122h_5*NCCO/3=O  3
Fuc     u1221m                3
Hex     uxxxxh                2
HexNAc  uxxxxh_2*NCC/3=O      1
dHex    uxxxxm                3
Pent    uxxxh                 3
"""
symbol2wurcs={}
wurcsorder={}
for l in symbol2wurcs_definition.splitlines():
    if not l.strip():
	continue
    sl = l.split()
    symbol2wurcs[sl[0]] = sl[1]
    wurcsorder[sl[1]] = int(sl[2])

for l in sys.stdin:
    l = l.strip()
    sl = re.split(r'\((\d+)\)',l)
    sl = map(str.strip,sl)
    comp = {}
    badsym = False
    for i in range(0,len(sl)-1,2):
	if int(sl[i+1]) == 0:
	    continue
	if sl[i] not in symbol2wurcs:
	    badsym = True
	    break
	comp[symbol2wurcs[sl[i]]] = int(sl[i+1])
    if badsym:
	continue
    skels = list(comp)
    total = sum(comp.values())
    if total == 0:
	continue
    uniq = len(skels)
    wurcsseq = "WURCS=2.0/%s,%s,%s/" % (uniq, total, "0+")
    wurcsseq += "".join(map(lambda sk: "[%s]" % (sk,), skels)) + "/"
    inds = []
    for i, sk in enumerate(sorted(skels,key=wurcsorder.get)):
        inds.extend([str(i + 1)] * comp[sk])
    wurcsseq += "-".join(inds)
    wurcsseq += "/"
    wurcsseq = gtc.fixcompwurcs(wurcsseq)
    try:
        status = gtc.register(wurcsseq)
    except:
        status = "issue"
    print l
    print status
    print wurcsseq
    print "-" * 100 + "\n"
