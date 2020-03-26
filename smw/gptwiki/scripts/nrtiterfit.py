#!/bin/env python27

from getwiki import GPTWiki

import sys, urllib, string, re, os, json, copy
from collections import defaultdict
from operator import itemgetter

from optparse import OptionParser

parser = OptionParser()
parser.add_option("--cachefile",type='string',dest='cachefile',default=None,
                  help="GPTwiki Transition Group cachefile. Default: No cache.")
parser.add_option("--choose",type="int",dest='choose',default=10,
		  help="Number of bad transition groups to choose per iteration. Default: 10.")
parser.add_option("--thresh",type="float",dest='thresh',default=5,
		  help="Residual threshold to determine outliers. Default: 5.")
parser.add_option("--heuristic",type="choice",dest="heuristic",default='weighted',
		  help="Bad TG selection heuristic. One of uniform, weighted, worst, all. Default: weighted",
                  choices=["uniform","weighted","worst","all"])
parser.add_option("--progress",action="store_true",dest='progress',default=False,
		  help="Display status at each iteration. Default: False.")
parser.add_option("--iterations",type="int",dest='iterations',default=1,
		  help="Number of times to iteratively fit. Default: 1.")
parser.add_option("--upload",action="store_true",dest='upload',default=False,
		  help="Upload TG status and peptide nrt to GPTwiki. Default: False.")

opts,args = parser.parse_args()

if opts.cachefile and os.path.exists(opts.cachefile):
    data = json.loads(open('gptwikitgs.cache').read())
    rows = data['rows']
    tgs = data['tgs']
    origpepnrt = data['pepnrt']
else:
    monos = "NHFS"
    tgs = dict()
    origpepnrt = dict()
    rows = []
    w = GPTWiki(quiet=True)
    for tg in w.itertransgroups():
        pepid = tg.get('peptide')
        p = w.get(pepid)
        pepacc = p.get('id')
        pepseq = p.get('sequence')
        pepname = p.get('name')
        pepnrt = p.get('nrt')
        if not tg.has('nrt'):
	    continue
        nrt = float(tg.get('nrt'))
        glyacc = p.get('glycan')[0][0]
        g = w.get(glyacc)
        gsym = g.get('sym')
        mcnt = {}
        for mono in monos:
	    mcnt[mono] = 0
            m = re.search(mono+r'(\d+)',gsym)
            if m:
	        mcnt[mono] = int(m.group(1))
            elif mono in gsym:
	        mcnt[mono] = 1
        mcnt['Total'] = sum(mcnt.values())
        mcnt['Ox'] = min(pepname.count('[Ox]'),1)
        rows.append((tg.get('id'),pepid,pepseq,gsym,mcnt,nrt))
        tgs[tg.get('id')] = dict(pepseq=pepseq,pepid=pepid)
	if pepnrt != None:
            origpepnrt[pepid] = float(pepnrt)

    if opts.cachefile:
        h = open(opts.cachefile,'w')
        h.write(json.dumps({'rows': rows, 'tgs': tgs, 'pepnrt': origpepnrt}))
        h.close()
    
pepseqcount = 0
pepseqindex = dict()
rowcount = 0
maxmcnt = defaultdict(int)
minmcnt = defaultdict(lambda: 1e+20)
monos = ("S","Ox")
for tgid,pepid,pepseq,gsym,mcnt,nrt in rows:
    if pepseq not in pepseqindex:
        pepseqindex[pepseq] = pepseqcount
        pepseqcount += 1
    for mono in monos:
	if mcnt[mono] > maxmcnt[mono]:
	    maxmcnt[mono] = mcnt[mono]    
	if mcnt[mono] < minmcnt[mono]:
	    minmcnt[mono] = mcnt[mono]    

# We skip zero so that the peptide coefficient is the "base" NRT...
monos = ("S","Ox")
monoindex = {}
j = pepseqcount
for mono in monos:
    for i in range(minmcnt[mono]+1,maxmcnt[mono]+1):
	monoindex[(mono,i)] = j
	j += 1
totalindex = j
j += 1
totalcols = j

print max(pepseqindex.values())
print monoindex, totalindex, totalcols

def getvalues(pepseq,mcnt):
    l = [ (pepseqindex[pepseq],1) ]
    for m in monos:
	if (m,mcnt[m]) in monoindex:
	    l.append((monoindex[(m,mcnt[m])],1))
    l.append((totalindex,mcnt['Total']))
    return l

def stop(goodresid,badresid):
    if max(map(abs,goodresid.values())) > opts.thresh:
	return False
    if min(map(abs,badresid.values())) <= opts.thresh:
	return False
    return True

def choosebad(goodresid):
    tgids = []
    values = []
    for tgid,resid in goodresid.items():
	if abs(resid) > opts.thresh:
	    tgids.append(tgid)
	    values.append(abs(resid))
    values = array(values)/sum(values)
    if len(tgids) == 0:
	return set()
    if opts.heuristic == "uniform":
        inds = random.choice(len(tgids),min(opts.choose,len(values)),replace=False)
    elif opts.heuristic == "weighted":
        inds = random.choice(len(tgids),min(opts.choose,len(values)),replace=False,p=values)
    elif opts.heuristic == "worst":
	sortedinds = sorted(range(len(values)),key=lambda i: values[i],reverse=True)
        inds = sortedinds[:opts.choose]
	# figure out if there are ties at the value of the last index chosen...
        criticalvalue = values[inds[-1]]
	criticalvalueinds = set()
	for i in range(len(tgids)):
	    if abs(values[i]-criticalvalue) < 1e-5:
		criticalvalueinds.add(i)
	if len(criticalvalueinds) > len(criticalvalueinds&set(inds)):
	    # print "inds",", ".join(map(str,sorted(inds)))
	    # print "crit",", ".join(map(str,sorted(criticalvalueinds)))
	    toremove = (set(inds)&criticalvalueinds)
	    # print "torm",", ".join(map(str,sorted(toremove)))
	    ntoreplace = len(toremove)
	    inds = set(inds)-toremove
	    # print "inds",", ".join(map(str,sorted(inds)))
	    replaceinds = random.choice(len(criticalvalueinds),ntoreplace,replace=False)
	    criticalvalueinds = list(criticalvalueinds)
	    inds = inds|set(map(lambda i: criticalvalueinds[i],replaceinds))
	    # print "inds",", ".join(map(str,sorted(inds)))
    elif opts.heuristic == "all":
	inds = list(range(len(values)))
    else:
	raise RuntimeError("bad heuristic value")
    return set(map(lambda i: tgids[i],inds))

def choosegood(badresid):
    goods = []
    for tgid in badresid:
	if abs(badresid[tgid]) < opts.thresh:
	    goods.append(tgid)
    return set(goods)

def status(iteration,badtg,newgoodtg,newbadtg,goodresid,beta):
    goodtg = set(goodresid)
    print "%3d"%(iteration,),"GTG:%4d(%3d)"%(len(goodtg),len(newgoodtg)), "BTG:%3d(%3d)"%(len(badtg),len(newbadtg)),
    for mono in monos:
        for i in range(1,maxmcnt[mono]+1):
            print "%s%d:%.2f"%(mono,i,beta[monoindex[(mono,i)]],),
    rms = np.sqrt(np.mean(array(goodresid.values())**2))
    print "#:%+.2f"%(beta[totalindex],),"RMS:%.2f"%(rms),
   
    pepseqcounts = {}
    pepidcounts = {}
    for tgid in tgs:
        pepseq = tgs[tgid]['pepseq']
        if pepseq not in pepseqcounts:
	    pepseqcounts[pepseq] = [0,0]
        pepid = tgs[tgid]['pepid']
        if pepid not in pepidcounts:
	    pepidcounts[pepid] = [0,0]
        if tgid in goodtg:
	    pepseqcounts[pepseq][0] += 1
	    pepidcounts[pepid][0] += 1
        if tgid in badtg:
	    pepseqcounts[pepseq][1] += 1
	    pepidcounts[pepid][1] += 1
   
    allgood = 0;allbad = 0;mixed = 0
    for pepseq in pepseqcounts:
        if pepseqcounts[pepseq][0] == 0:
	    allbad += 1
        elif pepseqcounts[pepseq][1] == 0:
	    allgood += 1
        else:
	    mixed += 1
    print "SEQ:%s,%s,%s"%(allgood,mixed,allbad),
    allgood1= 0;allbad1= 0;mixed1= 0
    for pepid in pepidcounts:
        if pepidcounts[pepid][0] == 0:
	    allbad1 += 1
        elif pepidcounts[pepid][1] == 0:
	    allgood1 += 1
        else:
	    mixed1 += 1
    print "PEP:%s,%s,%s"%(allgood1,mixed1,allbad1)
    return {'rms': rms, 'badtg': set(badtg), 'allbadpep': allbad1, 'beta': copy.copy(beta)}

def tgregress(badtg,beta=None):
    goodtg = alltg-badtg;
    A0 = zeros(shape=(len(goodtg),totalcols))
    A1 = zeros(shape=(len(badtg),totalcols))
    nrtvec0 = zeros(shape=(len(goodtg),))
    nrtvec1 = zeros(shape=(len(badtg),))
    goodtgids = []
    badtgids  = []
    a0rowindex = 0
    a1rowindex = 0
    for tgid,pepid,pepseq,gsym,mcnt,nrt in rows:
	if tgid in goodtg:
	    for j,v in getvalues(pepseq,mcnt):
	        A0[a0rowindex,j] = v
            nrtvec0[a0rowindex] = float(nrt)
	    goodtgids.append(tgid)
	    a0rowindex += 1
        else:
	    for j,v in getvalues(pepseq,mcnt):
	        A1[a1rowindex,j] = v
            nrtvec1[a1rowindex] = float(nrt)
	    badtgids.append(tgid)
	    a1rowindex += 1
    
    returnbeta=False
    if isinstance(beta,type(None)):
        beta = linalg.lstsq(A0,nrtvec0)[0]
	returnbeta = True
    nrtfit0 = A0.dot(beta)
    goodresid = dict(zip(goodtgids,(nrtvec0-nrtfit0)))
    nrtfit1 = A1.dot(beta)
    badresid = dict(zip(badtgids,(nrtvec1-nrtfit1)))
    if returnbeta:
        return beta,goodresid,badresid
    else:
        return goodresid,badresid

def pepregress(pepids,beta,nrts):
    A = zeros(shape=(len(pepids),totalcols))
    nrtvec = zeros(shape=(len(pepids),))
    thepepids = []
    rowindex = 0
    seen = set()
    for tgid,pepid,pepseq,gsym,mcnt,nrt in rows:
	if pepid in seen or pepid not in pepids or pepid not in nrts:
	    continue
	seen.add(pepid)
	for j,v in getvalues(pepseq,mcnt):
	    A[rowindex,j] = v
        nrtvec[rowindex] = float(nrts[pepid])
	thepepids.append(pepid)
	rowindex += 1
    
    nrtfit = A.dot(beta)
    resid = dict(zip(thepepids,(nrtvec-nrtfit)))
    return resid

from numpy import zeros, linalg, array, random
import numpy as np

random.seed()

alltg = set(map(itemgetter(0),rows))
badtgstats = []

for iters in range(opts.iterations):

    iteration = 1
    badtg = set()
    while True:

	beta,goodresid,badresid = tgregress(badtg)
	# print beta
    
        done = stop(goodresid,badresid)
        # print done

        newbadtg = choosebad(goodresid)
        newgoodtg = choosegood(badresid)

	# print newbadtg, newgoodtg
	
        if len(newbadtg) == 0 and len(newgoodtg) == 0:
	    done = True

	if done:
	    break
    
        if opts.progress:
	    status(iteration,badtg,newgoodtg,newbadtg,goodresid,beta)
    
        iteration += 1
        badtg |= newbadtg
        badtg -= newgoodtg

    # print "here"
    stats = status(iteration,badtg,newgoodtg,newbadtg,goodresid,beta)
    badtgstats.append(stats)

bestallbadstats = 1e+20
beststats = None
for stats in badtgstats:
    if stats['allbadpep'] < bestallbadstats:
	beststats = stats
        bestallbadstats = stats['allbadpep']
beta = beststats['beta']
badtg = beststats['badtg']

goodresid,badresid = tgregress(badtg,beta)
newbadtg = map(lambda t: t[0],filter(lambda t: abs(t[1]) > opts.thresh, goodresid.items()))
newgoodtg = map(lambda t: t[0],filter(lambda t: abs(t[1]) <= opts.thresh, badresid.items()))
status(-1,badtg,set(),newbadtg,goodresid,beta)

newpepnrt = defaultdict(list)
for tgid,pepid,pepseq,gsym,mcnt,nrt in rows:
    if tgid in badtg:
	continue
    newpepnrt[pepid].append(nrt)
for pepid in list(newpepnrt):
    newpepnrt[pepid] = np.median(array(newpepnrt[pepid]))
oldresid = pepregress(set(newpepnrt),beta,origpepnrt)
newresid = pepregress(set(newpepnrt),beta,newpepnrt)
for pepid in oldresid:
    if (abs(newresid[pepid] - oldresid[pepid]) > 0.1) or (abs(newpepnrt[pepid]-origpepnrt[pepid]) > 0.1): 
        print pepid,origpepnrt[pepid],oldresid[pepid],newpepnrt[pepid],newresid[pepid]

if not opts.upload:
   sys.exit(0)

w = GPTWiki(quiet=True)
for tg in w.itertransgroups():
    tgid = tg.get('id')
    if tgid in badtg:
        tg.set('pepnrt','drop')
    elif tgid in alltg:
        tg.set('pepnrt','use')
    else:
	tg.delete('pepnrt')
    if w.put(tg):
	print tgid

for pep in w.iterpeptides():
    pepid = pep.get('id')
    if pepid in origpepnrt:
	if pepid in newpepnrt:
	    pep.set('nrt',newpepnrt[pepid])
	else:
	    pep.delete('nrt')
    else:
	pep.delete('nrt')
    if w.put(pep):
	print pepid
