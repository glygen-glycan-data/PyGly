#!/bin/env python27
import sys, csv, re
from collections import Counter, defaultdict
import itertools 
import numpy as np
import copy
import math

__all__ = [ 'LatticeFDR' , 'LatticeDim', 'RoundedLatticeDim',
            'AbsRoundedLatticeDim',
            'IntegerLatticeDim', 'NaturalLatticeDim',
            'RoundedLogLatticeDim', 'INCREASING', 'DECREASING',
            'AdaptiveLatticeDim','AdaptiveAbsLatticeDim']

INCREASING=-1
DECREASING=+1
NOCONVERSION=(lambda x: x)

class LatticeDim(object):
    preprocess = False
    def __init__(self,tag,
                 dirn=INCREASING,
                 conv=NOCONVERSION,
                 min=-1e+20,
                 max=+1e+20,
                 default=None,
                 display=None):
        self.tag = tag
        self.dir = dirn
        self.conv = conv
        self.min = min
        self.max = max
        self.default = default
        self.values = set()
	self.origvalues = defaultdict(set)
        self.ordered = False
        if display != None:
            self.display = display
        else:
            self.display = self._display
    def nvalues(self):
        return len(self.values)
    def value(self,it):
        if self.ordered:
            self.values = set(self.values)
            self.ordered = False
        try:
            val = it[self.tag]
            val = min(max(val,self.min),self.max)
        except KeyError:
            val = self.default
        origval = val
        val = self.conv(val)
        if self.dir == DECREASING:
            val = -val
        self.values.add(val)
	self.origvalues[val].add(origval)
        return val
    def binvalues(self,freq,fdrfunc):
        return
    def order(self):
        if not self.ordered:
            self.values = sorted(self.values)
            self.ordered = True
    def _display(self,x):
        if self.dir == DECREASING:
            return -x
        else:
            return x
    def thresh(self,x):
	if self.dir == DECREASING:
	    return min(self.origvalues[x])
	else:
	    return max(self.origvalues[x])

class IntegerLatticeDim(LatticeDim):
    def __init__(self,tag,**kw):
        kw['conv'] = int
        LatticeDim.__init__(self,tag,**kw)

class NaturalLatticeDim(IntegerLatticeDim):
    def __init__(self,tag,**kw):
        kw['min'] = 0
        IntegerLatticeDim.__init__(self,tag,**kw)

class BooleanIntLatticeDim(LatticeDim):
    def __init__(self,tag,**kw):
	kw['conv'] = lambda x: (int(x)!=0)
        LatticeDim.__init__(self,tag,**kw)

class RoundedLatticeDim(LatticeDim):
    def __init__(self,tag,places=0,scale=None,**kw):
	if scale!=None:
	    places = math.log(1.0/scale,10.0)
        kw['conv'] = lambda x: round(x*(10**places))
        kw['display'] = lambda x: (x if self.dir == INCREASING else -x)/(10.0**places)
        LatticeDim.__init__(self,tag,**kw)

class AdaptiveLatticeDim(LatticeDim):
    preprocess = True
    def __init__(self,tag,points=10,**kw):
        kw['conv'] = self.conv
        kw['display'] = self.display
        LatticeDim.__init__(self,tag,**kw)
        self.preprocessdone = False
        self.points = points

    @staticmethod
    def fdrdump(msg,fdr):
        print ">>",msg
	for thr,val in sorted(fdr.items()):
	    print "%.6f\t%.6f"%(thr,val)
	print
	
    def binvalues(self,freq,fdrfunc=None):
	print "===========\n"+self.tag+"\n==========="
        nt=Counter(); nf=Counter(); fdr={}
        lastthr = None
        for thr in sorted(freq):
            nt[thr] = nt[lastthr] + freq[thr][False]
            nf[thr] = nf[lastthr] + freq[thr][True]
            if fdrfunc != None:
                fdr[thr] = fdrfunc(nt[thr],nf[thr])
            else:
                fdr[thr] = float(nf[thr])/float(nt[thr])
            lastthr = thr

        lastfdr = 1.0
        for thr in sorted(freq,reverse=True):
            if fdr[thr] > lastfdr:
                fdr[thr] = lastfdr
            lastfdr = fdr[thr]

	self.fdrdump("Monotonic",fdr)

	if len(fdr) > 0:
          lastthr = min(fdr)
          for thr in sorted(fdr)[1:]:
	    if (fdr[lastthr]==fdr[thr]):
		del fdr[lastthr]
            # elif fdr[lastthr]>(fdr[thr]-0.01):
            #    del fdr[thr]
	    #    continue
            lastthr = thr

	self.fdrdump("Spaced",fdr)

        self.critical_points = sorted(fdr)
        if len(self.critical_points) > self.points:
            self.critical_points = self.critical_points[:(self.points-1)]+[self.critical_points[-1]]
        print "Critical points:",len(self.critical_points)," ".join(map(str,self.critical_points))
        self.preprocessdone = True

    def conv(self,val):
        if not self.preprocessdone:
            return val
        if self.dir == DECREASING:
            for v in self.critical_points:
                if -val <= v:
                    return -v
            return -self.critical_points[-1]
        else:
            for v in self.critical_points:
                if val <= v:
                    return v
            return self.critical_points[-1]
    def display(self,val):
        return val

class AdaptiveAbsLatticeDim(AdaptiveLatticeDim):
    def __init__(self,tag,**kw):
        kw['conv'] = self.conv
        kw['display'] = self.display
        super(AdaptiveAbsLatticeDim,self).__init__(tag,**kw)
    def conv(self,val):
        return super(AdaptiveAbsLatticeDim,self).conv(abs(val))
    def display(self,val):
        return abs(val)

class RoundedLogLatticeDim(LatticeDim):
    def __init__(self,tag,places=0,scale=None,base=10,**kw):
	if scale!=None:
	    places = math.log(1.0/scale,10.0)
        kw['conv'] = lambda x: round(math.log(x,base)*(10**places))
        kw['display'] = lambda x: round(math.exp((x if self.dir == INCREASING else -x)/10.0**places,base),3)
        LatticeDim.__init__(self,tag,**kw)

class AbsRoundedLatticeDim(LatticeDim):
    def __init__(self,tag,places=0,scale=None,**kw):
	if scale!=None:
	    places = math.log(1.0/scale,10.0)
        kw['conv'] = lambda x: round(abs(x)*(10**places))
        kw['display'] = lambda x: (x if self.dir == INCREASING else -x)/(10.0**places)
        LatticeDim.__init__(self,tag,**kw)

class LatticeFDR:
    def __init__(self,fdrfunc=None):
        self.dims = []
        self.fdrfunc = fdrfunc
    def add_dimension(self,dim):
        self.dims.append(dim)
    def set_decoy(self,dim):
        self.decoy = dim
    def compute(self,items,targetfdr=None,matrix=False):
        
        self.ndim = len(self.dims)

        # print repr(items)

        for d in self.dims:
            if not d.preprocess:
                continue
            # print d.tag,
            freq = defaultdict(Counter)
            for it in items:
                freq[d.value(it)][self.decoy.value(it)] += 1
            # print freq
            # print min(freq),
            d.binvalues(freq,self.fdrfunc)
            # print max(freq)

        c = Counter()
        for it in items:
            key = tuple(map(lambda d: d.value(it),self.dims)+[self.decoy.value(it)])
            # print key
            c[key] += 1

        for d in self.dims:
            if d.preprocess:
                d.values = d.critical_points

        self.ndims = []
        for d in self.dims:
            self.ndims.append(d.nvalues())
            d.order()
        self.ndims = tuple(self.ndims)

	print >>sys.stderr, "Lattice dimensions:", self.ndims

	if min(self.ndims) == 0:
	    print >>sys.stderr, "No observations!"
	    raise RuntimeError("No spectral hits")
    
        target = np.zeros(dtype=np.float32,shape=self.ndims)
        decoy = np.zeros(dtype=np.float32,shape=self.ndims)

        for i in self.allinds():
            key = self.ind2key(i)
            target[i] = c[tuple(key+[False])]
            decoy[i] = c[tuple(key+[True])]

        for i in range(self.ndim):
            target = target.cumsum(axis=i)
            decoy = decoy.cumsum(axis=i)

        fdr = np.zeros(dtype=np.float32,shape=self.ndims)
        for i in self.allinds():
            if self.fdrfunc:
                fdr[i] = self.fdrfunc(target[i],decoy[i])
            else:
                fdr[i] = decoy[i]/target[i]

        for i in self.allindsrev():
            ub = 1.0
            for i1 in self.allneighbors(i):
                if fdr[i1] < ub:
                    ub = fdr[i1]
            if fdr[i] > ub:
                fdr[i] = ub

        if matrix:
	    np.savetxt('fdr.txt',fdr)
            np.savetxt('target.txt',target)

        dp = np.zeros(shape=self.ndims,dtype=np.float32)
        dr = np.zeros(shape=self.ndims,dtype=np.uint8)

        dp[i] = 0
        dr[i] = 0
        for i in self.allinds():
            if sum(i) == 0:
                continue
            edges = [(0,None)]
            for i0 in self.allneighborsrev(i):
                # dpi0 = dp[i0] + (fdr[i]-fdr[i0])*target[i0]
                dpi0 = target[i0]
                dri0 = self.encode(i,i0)
                edges.append((dpi0,dri0))
            dp[i],dr[i] = max(edges)

	traceinds = []
        i = tuple([-1]*self.ndim)
        if targetfdr != None:
	    i = np.unravel_index(np.argmax(target[fdr<=targetfdr]),target.shape)
	    i = tuple(map(lambda ii: -(self.ndims[ii]-i[ii]),range(self.ndim)))
        if matrix:
	    traceinds.append(map(lambda ii: self.ndims[ii]+i[ii], range(self.ndim)))
        lastfdr = None
        trace = [self.ind2key(i) + [fdr[i],target[i],decoy[i],dp[i]]]
        while True:
           if fdr[i] != lastfdr and trace[-1][:self.ndim] != self.ind2key(i):
               row = self.ind2key(i) + [fdr[i],target[i],decoy[i],dp[i]]
               trace.append(row)
           lastfdr = fdr[i]
           if dr[i] == 0:
               break
           i = self.decode(i,dr[i])
           if matrix:
	       traceinds.append(map(lambda ii: self.ndims[ii]+i[ii], range(self.ndim)))
        if trace[-1][:self.ndim] != self.ind2key(i):
            if matrix:
	        traceinds.append(map(lambda ii: self.ndims[ii]+i[ii], range(self.ndim)))
            trace.append(self.ind2key(i) + [fdr[i],target[i],decoy[i],dp[i]])

        if matrix:
	    np.savetxt('fdr.txt',fdr,fmt="%.2f")
            np.savetxt('target.txt',target,fmt="%3d")
	    wh = open('trace.txt','w')
	    for t in reversed(traceinds):
           	print >>wh, " ".join(map(str,t)) 
	    wh.close()

        self.trace = list(reversed(trace))
        return self.trace

    def fdr(self,it):
        key = map(lambda d: d.value(it),self.dims)
        for t in self.trace:
            if key <= t[:self.ndim]:
                return float(t[self.ndim])
        return 1.0

    def tracekeys(self):
        return self.keys() + ["FDR","Target","Decoy","AUC"]

    def keys(self):
        return [d.tag for d in self.dims]

    def dump(self):
        tkeys = self.tracekeys()
        width = map(len,tkeys)
        rows = []
        for t in self.trace:
            l = []
            for i in range(self.ndim):
                l.append(("%s"%self.dims[i].thresh(t[i]))[:5])
	    l.append("%.5f"%t[-4])
            for i in range(-3,0):
                l.append("%s"%int(t[i]))
            for i,v in enumerate(l):
                if len(v) > width[i]:
                    width[i] = len(v)
            rows.append(l)
        print " ".join(map(lambda t: "%*s"%t,zip(width,tkeys)))
        for r in rows:
            print " ".join(map(lambda t: "%*s"%t,zip(width,r)))

    def allinds(self):
        return itertools.product(*map(range,self.ndims))

    def allindsrev(self):
        return itertools.product(*map(lambda n: range(n-1,-1,-1),self.ndims))

    def allneighbors(self,i):
        for k in range(1,len(i)+1):
            for d in itertools.combinations(range(len(i)),k):
                bad = False
                for di in d:
                    if i[di]+1 >= self.ndims[di]:
                        bad = True
                if bad:
                    continue
                i1 = list(i)
                for di in d:
                    i1[di] += 1
                yield tuple(i1)

    def allneighborsrev(self,i):
        for k in range(1,len(i)+1):
            for d in itertools.combinations(range(len(i)),k):
                bad = False
                for di in d:
                    if i[di]-1 < 0:
                        bad = True
                if bad:
                    continue
                i1 = list(i)
                for di in d:
                    i1[di] -= 1
                yield tuple(i1)

    def ind2key(self,i):
        key = []
        for j,d in enumerate(self.dims):
            key.append(d.values[i[j]])
        return key

    @staticmethod
    def encode(i1,i0):
        retval = 0
        for j in range(len(i1)):
            if i1[j] > i0[j]:
                retval += 2**j
        return retval

    @staticmethod
    def decode(i,d0):
        i0 = list(i)
        for j in range(len(i)):
            if d0&(2**j):
               i0[j] -= 1
        return tuple(i0)

def get_rows(fn):
    groupingkey = "peptide_group_label"
    rankingkey = "main_var_xx_swath_prelim_score"
    rows = csv.DictReader(open(fn),dialect="excel-tab")
    last_group = None                                                                                                   
    for row in sorted(rows,key=lambda r: (r[groupingkey],-1*float(r[rankingkey]))):
        group = row[groupingkey]
        if group == last_group:
            continue
        # print row
        # if not row.get('delta_rt'):
	#     row['delta_rt'] = float(row['RT'])-float(row['assay_rt'])
	for k in row:
	    try:
                row[k] = float(row[k])
            except ValueError:
		pass
        row['decoy'] = (int(row['decoy']) == 1)
        yield row
        last_group = group

if __name__ == "__main__":

    import sys
    ndecoys = int(sys.argv[2])

    for l in open(sys.argv[1]):
	sl = l.split()
	scores = filter(lambda x: re.search(r'score',x),sl)
	break

    badscores = """
	var_bseries_score
	var_manhatt_score
	var_intensity_score
	var_isotope_overlap_score
	var_massdev_score
	var_massdev_score_weighted
	var_norm_rt_score
	var_yseries_score
	var_elution_model_fit_score
	rt_score
	xx_swath_prelim_score
    """.split()

    for bsc in badscores:
	scores.remove(bsc)

    # scores = """
    #	xx_lda_prelim_score
    #	main_var_xx_swath_prelim_score
    # 	delta_rt
    # """.split()
    
    fdrfunc = (lambda t,d: (float(d)/ndecoys)/float(t) if t > 0 else 1e+20)

    fdr = LatticeFDR(fdrfunc=fdrfunc)
    for s in scores:
        fdr.add_dimension(AdaptiveLatticeDim(s,dirn=DECREASING))
    fdr.set_decoy(LatticeDim("decoy"))
    fdr.compute(list(get_rows(sys.argv[1])),targetfdr=0.10)
    fdr.dump()
