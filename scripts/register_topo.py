#!/bin/env python2

import findpygly
# from pygly.GlyTouCan import GlyTouCan
from pygly.alignment import GlycanEqual
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format
from pygly.manipulation import Topology
from pygly.GlycanResource import GlyTouCan, GlyTouCanNoPrefetch, GlyTouCanNoCache
import sys, re, time
import urllib
import hashlib

# gtc = GlyTouCanNoPrefetch()
gtc = GlyTouCanNoCache()
topo = Topology()
wp = WURCS20Format()
gp = GlycoCTFormat()
geq = GlycanEqual()

alpha2ind = dict()
ind2alpha = dict()
for i in range(ord('a'),ord('z')+1):
    alpha2ind[chr(i)] = i-ord('a')+1
    ind2alpha[i-ord('a')+1] = chr(i)
for i in range(ord('A'),ord('Z')+1):
    alpha2ind[chr(i)] = i-ord('A')+1+26
    ind2alpha[i-ord('A')+1+26] = chr(i)

def totopo(wseq):
    # print(wseq)
    head,cnts,rest = wseq.split('/',2)
    monos,rest = rest.split(']/',1)
    inds,links = rest.split('/',1)

    cnts = list(map(int,cnts.split(',')))
    inds = list(map(int,inds.split('-')))
    monos = monos.lstrip('[').rstrip(']').split('][')
    links = links.split('_')

    # print(head,cnts,monos,inds,links)

    if any(['~' in l for l in links]):
        return None

    for i in range(len(monos)):
        m = monos[i]
        sm = m.split('_')
        try:
            base,anomer = sm[0].split('-')
        except ValueError:
            continue
        if not re.search(r'^\d[abx]$',anomer):
            continue
        anomer = anomer[0] + 'x'
        monos[i] = base+'-'+anomer+'_'+'_'.join(sm[1:])

    indsmap = dict()
    newmonos = list(monos)
    for i in range(len(newmonos)-1,-1,-1):
        try:
            ind = newmonos[:i].index(newmonos[i])
        except ValueError:
            continue
        del newmonos[i]

    for i in range(len(monos)):
        indsmap[i+1] = newmonos.index(monos[i])+1
    inds = [ indsmap.get(ind,ind) for ind in inds ]

    monos = newmonos
    for i in range(len(links)):
        l = links[i].split('-')
        if not re.search(r'^([a-zA-Z]\d\|)*[a-zA-Z]\d$',l[0]):
            continue
        if not re.search(r'^([a-zA-Z]\d\|)*[a-zA-Z]\d$',l[1]):
            continue
        l0 = l[0].split('|')
        l0mi = min(alpha2ind[li[0]] for li in l0)
        l0ma = max(alpha2ind[li[0]] for li in l0)
        l1 = l[1].split('|')
        l1ma = max(alpha2ind[li[0]] for li in l1)
        l1mi = min(alpha2ind[li[0]] for li in l1)
        if l0ma < l1mi:
           for j in range(len(l0)):
               l0[j] = l0[j][0]+'?'
           l[0] = "|".join(sorted(set(l0)))
        elif l1ma < l0mi:
           for j in range(len(l1)):
               l1[j] = l1[j][0]+'?'
           if len(set(l1))==1:
               l[0] = "|".join(sorted(set(l1)))
               l[1] = "|".join(l0)
           else:
               l[1] = "|".join(sorted(set(l1)))
        else:
            return None
        links[i] = '-'.join(l)

    cnts = [len(monos),len(inds),len(links)]
    # print(head,cnts,monos,inds,links)

    cnts = ",".join(map(str,cnts))
    monos = '[' + ']['.join(monos) + ']'
    inds = '-'.join(map(str,inds))
    links = '_'.join(links)
    return '/'.join([head,cnts,monos,inds,links])

seenhash = set()
count = 1
for acc0 in sys.stdin:
    acc0 = acc0.strip()
    wurcsseq = gtc.getseq(acc0,format='wurcs')
    if ',0/' in wurcsseq or '_1*N_' in wurcsseq:
        continue
    # print(acc0,wurcsseq)
    gly = wp.toGlycan(wurcsseq)
    if not gly:
        continue
    topowseq = totopo(wurcsseq)
    if not topowseq:
        continue
    wtgly = wp.toGlycan(topowseq)
    if not wtgly:
        continue
    tgly = topo(gly)
    if not geq.eq(wtgly,tgly):
        # print(wtgly.glycoct())
        # print(tgly.glycoct())
        # sys.exit(1)
        continue
    # newseq = tgly.glycoct()
    newseq = topowseq
    thehash = hashlib.sha256(newseq).hexdigest().lower()
    if thehash in seenhash:
        continue
    seenhash.add(thehash)
    hash,acc,error = gtc.gethashedseq(seq=newseq)
    if not hash:
	hash = gtc.register(newseq)
	if hash:
           print count,acc0,hash,"*"
	else:
           print count,acc0,"???","*"
	time.sleep(10)
    elif not acc:
	if error:
            print count,acc0,repr(error)
	else:
            print count,acc0,hash
    else:
	print count,acc0,acc
    count += 1
