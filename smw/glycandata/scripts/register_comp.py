#!/bin/env python2

import findpygly
# from pygly.GlyTouCan import GlyTouCan
from pygly.GlycanResource import GlyTouCan as GTC
import sys, re, time
import urllib
import hashlib

gtc = GTC(prefetch=False)

# Need wurcs skeleton codes for UniCarbKB and Byonic symbols...

symbol2wurcs_definition = """
NeuAc   AUd21122h_5*NCC/3=O   1
NeuGc   AUd21122h_5*NCCO/3=O  2
Fuc     u1221m                6
Hex     uxxxxh                5
HexNAc  uxxxxh_2*NCC/3=O      4
dHex    uxxxxm                3
Pent    uxxxh                 10
P	*OPO/3O/3=O	      -1
Phospho	*OPO/3O/3=O	      -1
S	*OSO/3=O/3=O	      -1
Sulpho	*OSO/3=O/3=O	      -1
"""
symbol2wurcs={}
wurcsorder={}
for l in symbol2wurcs_definition.splitlines():
    if not l.strip():
	continue
    sl = l.split()
    symbol2wurcs[sl[0]] = sl[1]
    wurcsorder[sl[1]] = int(sl[2])

for lineno,l in enumerate(sys.stdin):
    l = l.strip()
    l0 = l
    if l.startswith('comp_'):
	l = l[5:]
    if '(' in l:
        sl = re.split(r'\s*\((\s*\d+\s*)\s*\)',l)
    else:
	sl = re.split(r'\s*(\d+)\s*',l)
    sl = map(str.strip,sl)
    comp = {}
    subst = {}
    badsym = False
    for i in range(0,len(sl)-1,2):
	if int(sl[i+1]) == 0:
	    continue
	if sl[i] not in symbol2wurcs:
	    badsym = True
	    break
	skel = symbol2wurcs[sl[i]]
	cnt = int(sl[i+1])
	if wurcsorder.get(skel) < 0:
	    subst[skel] = cnt
	else:
	    comp[skel] = cnt
    if badsym:
	print lineno+1,l0,"Bad symbol:",sl[i]
	continue
    skels = sorted(comp,key=wurcsorder.get)
    total = sum(comp.values())
    if total == 0:
	continue
    if total == 1 and len(subst) > 0:
	newskel = skels[0]
	for subst,cnt in subst.items():
	    for i in range(cnt):
	        newskel += "_?" + subst;
	comp = {newskel:1}
	skels = [ newskel ]
	subst = {}
    uniq = len(skels)
    wurcsseq = "WURCS=2.0/%s,%s,%s/" % (uniq, total, "0+")
    wurcsseq += "".join(map(lambda sk: "[%s]" % (sk,), skels)) + "/"
    inds = []
    for i, sk in enumerate(skels):
        inds.extend([str(i + 1)] * comp[sk])
    wurcsseq += "-".join(inds)
    wurcsseq += "/"
    wurcsseq = gtc.fixcompwurcs(wurcsseq,subst)
    thehash = hashlib.sha256(wurcsseq).hexdigest().lower()
    hash,acc,error = gtc.gethashedseq(seq=wurcsseq)
    if not hash:
	hash = gtc.register(wurcsseq)
	if hash:
           print lineno+1,l0,hash
	else:
           print lineno+1,l0,None,wurcsseq,thehash
	time.sleep(60)
    elif not acc:
	if error:
            print lineno+1,l0,repr(error)
	else:
            print lineno+1,l0,hash
    else:
	print lineno+1,l0,acc
