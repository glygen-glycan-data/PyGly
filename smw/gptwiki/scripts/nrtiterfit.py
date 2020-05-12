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
    data = json.loads(open(opts.cachefile).read())
    tgrows = data['tgrows']
    peprows = data['peprows']
    tgs = data['tgs']
    origpepnrt = data['pepnrt']
else:
    monos = "NHFS"
    tgs = dict()
    origpepnrt = dict()
    tgrows = []
    peprows = []
    pepseen = set()
    w = GPTWiki(quiet=True)
    for tg in w.itertransgroups():
        specname = tg.get('spectra')
        spec = w.get(specname)
	if spec.get('type',"DDA") != "DDA":
	    continue
        pepid = tg.get('peptide')
        p = w.get(pepid)
        pepacc = p.get('id')
        pepseq = p.get('sequence')
        pepname = p.get('name')
        pepnrt = p.get('nrt')
	pepnrtobs = int(p.get('nrtobs',0))
	if pepnrtobs == 0:
	    pepnrt = None
        nrt = tg.get('nrt')
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
	if pepid not in pepseen:
	    peprows.append((pepid,pepseq,gsym,mcnt,pepnrt))
	    pepseen.add(pepid)
	    if pepnrt != None:
                origpepnrt[pepid] = float(pepnrt)
        if nrt != None:
            tgrows.append((tg.get('id'),pepid,pepseq,gsym,mcnt,float(nrt)))
            tgs[tg.get('id')] = dict(pepseq=pepseq,pepid=pepid)

    if opts.cachefile:
        h = open(opts.cachefile,'w')
        h.write(json.dumps({'tgrows': tgrows, 'peprows': peprows, 'tgs': tgs, 'pepnrt': origpepnrt}))
        h.close()
    
pepseqcount = 0
pepseqindex = dict()
rowcount = 0
maxmcnt = defaultdict(int)
minmcnt = defaultdict(lambda: 1e+20)
monos = ("S","Ox")
for tgid,pepid,pepseq,gsym,mcnt,nrt in tgrows:
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
    values = np.array(values)/sum(values)
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
    rms = np.sqrt(np.mean(np.array(goodresid.values())**2))
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
    for tgid,pepid,pepseq,gsym,mcnt,nrt in tgrows:
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
        beta = linalg.lstsq(A0,nrtvec0,rcond=1e-15)[0]
	returnbeta = True
    nrtfit0 = A0.dot(beta)
    goodresid = dict(zip(goodtgids,(nrtvec0-nrtfit0)))
    nrtfit1 = A1.dot(beta)
    badresid = dict(zip(badtgids,(nrtvec1-nrtfit1)))
    if returnbeta:
        return beta,goodresid,badresid
    else:
        return goodresid,badresid

def pepregress(pepids,beta,nrts=None):
    A = zeros(shape=(len(pepids),totalcols))
    nrtvec = zeros(shape=(len(pepids),))
    thepepids = []
    rowindex = 0
    for pepid,pepseq,gsym,mcnt,nrt in peprows:
	if pepid not in pepids or (nrts != None and pepid not in nrts) or pepseq not in pepseqindex:
	    continue
	for j,v in getvalues(pepseq,mcnt):
	    A[rowindex,j] = v
	if nrts != None and pepid in nrts:
            nrtvec[rowindex] = float(nrts[pepid])
	else:
	    nrtvec[rowindex] = 0.0
	thepepids.append(pepid)
	rowindex += 1
    
    nrtfit = A.dot(beta)
    fits = dict(zip(thepepids,nrtfit))
    if nrts != None:
        resid = dict(zip(thepepids,(nrtvec-nrtfit)))
        return resid
    return fits

from numpy import zeros, linalg, random
import numpy as np

random.seed()

alltg = set(map(itemgetter(0),tgrows))
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
status(-1,badtg,newgoodtg,newbadtg,goodresid,beta)

newpepnrtlst = defaultdict(list)
newpepnrtobs = dict()
newpepnrt = dict()
goodpeps = set()
for tgid,pepid,pepseq,gsym,mcnt,nrt in tgrows:
    if tgid in badtg:
	continue
    newpepnrtlst[pepid].append(nrt)
    goodpeps.add(pepseq)

goodpepids = set()
for pepid,pepseq,gsym,mcnt,nrt in peprows:
    if pepseq in goodpeps:
	goodpepids.add(pepid)

for pepid in list(newpepnrtlst.keys()):
    newpepnrtobs[pepid] = len(newpepnrtlst[pepid])
    newpepnrt[pepid] = float(np.median(np.array(newpepnrtlst[pepid])))

oldresid = pepregress(set(origpepnrt),beta,origpepnrt)
newresid = pepregress(set(newpepnrt),beta,newpepnrt)
pepfits = pepregress(goodpepids,beta)
for pepid in sorted((set(origpepnrt)&set(newpepnrt)&set(pepseqindex))):
    if (abs(newresid[pepid] - oldresid[pepid]) > 0.01) or (abs(newpepnrt[pepid]-origpepnrt[pepid]) > 0.01): 
        print pepid,round(origpepnrt[pepid],3),round(oldresid[pepid]),round(newpepnrt[pepid],3),round(newresid[pepid],3)
for pepid in sorted((set(origpepnrt)-set(newpepnrt))):
    print pepid,round(origpepnrt[pepid],3),round(oldresid.get(pepid,1e+20),3),"lost"
for pepid in sorted((set(newpepnrt)-set(origpepnrt))):
    print pepid,round(newpepnrt[pepid],3),round(newresid[pepid],3),"gained"

if not opts.upload:
    sys.exit(0)

if raw_input("Upload? ") not in ("Y","y"):
    sys.exit(0)

w = GPTWiki(quiet=True)
for tg in w.itertransgroups():
    tgid = tg.get('id')
    if tgid in badtg:
        tg.set('pepnrt','drop (%+.3f)'%(badresid[tgid]))
    elif tgid in alltg:
        tg.set('pepnrt','use (%+.3f)'%(goodresid[tgid]))
    else:
	tg.delete('pepnrt')
    if opts.upload and w.put(tg):
	print tgid

for pep in w.iterpeptides():
    pepid = pep.get('id')
    if newpepnrt.get(pepid) != None:
	pep.set('nrt',newpepnrt[pepid])
	pep.set('nrtobs',newpepnrtobs[pepid])
    elif False and pepid in goodpepids and pepfits.get(pepid) != None:
	pep.set('nrt',pepfits[pepid])
	pep.set('nrtobs',0)
    else:
	pep.delete('nrt')
	pep.delete('nrtobs')
    if opts.upload and w.put(pep):
	print pepid
