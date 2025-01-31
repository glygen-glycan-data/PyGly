
import math, copy
from . import odometer
from operator import itemgetter
from collections import defaultdict

class IsoShape:

    def __init__(self,isotopes,composition,maxpos=5):

        self.table = []
        o = odometer.composite_odometer()
        elements = list(composition.keys())
        nElements = len(elements)
        o.set_size(nElements)
        for i,e in enumerate(elements):
            nIso = len(isotopes[e])
            o.get_value(i).set_size(nIso-1)
            o.get_value(i).set_max(composition[e])
            o.get_value(i).set_total_max(composition[e])
            for j in range(nIso-1):
                o.get_value(i).set_weighted_total_weights(j+1,j)
        o.set_min(0)
        o.set_weighted_total_max(maxpos)

        probvec = {}
        massvec = {}
        for i,e in enumerate(elements):
            probvec[e] = list(map(itemgetter(1),isotopes[e]))
            for ii,p in enumerate(probvec[e]):
                if ii == 0:
                    continue
                if p == 0.0:
                    o.get_value(i).set_max(0,ii-1)
            massvec[e] = list(map(itemgetter(0),isotopes[e]))

        o.init()
        self.mono = None
        while o.inrange():
            # o.write(sys.stdout)
            # print
            lnp = 0
            m = 0
            isos = {}
            for i,e in enumerate(elements):
                iso = copy.deepcopy(o.get_value(i).values())
                iso.insert(0,composition[e]-o.get_value(i).sum())
                isos[e] = iso
                lnp += self.lnmultinomial(probvec[e],iso,composition[e])
                m += self.mass(massvec[e],iso)
            if lnp > math.log(1e-5):
                self.table.append((math.exp(lnp),m,isos))
            if not self.mono:
                self.mono = m
            o.inc()

    def write(self,h,ch=1):
        t=defaultdict(float); am=defaultdict(float);aam=0.0;totp=0.0
        for (p,m,i) in self.table:
            # print >>h, p,m,i
            t[int(m*ch-self.mono+.5)] += p
            am[int(m*ch-self.mono+.5)] += m*ch*p
            aam += m*ch*p
            totp += p
            # print >>h, int(m*ch-self.mono+.5), m*ch, p, i
        # print >>h, "    %.4f"%(self.mono/ch)
        maxprob = max(t.values())
        totit = 0
        whtmz = 0
        for i,(m,p) in enumerate(sorted(t.items())):
            am[i] /= p
            if i == 0:
                print >>h, "%3d %10.4f %8.4f %8.4f %7.3f%%"%(i,am[i]/ch,(am[i]-self.mono)/ch,0.0,100*p/maxprob)
            else:
                print >>h, "%3d %10.4f %8.4f %8.4f %7.3f%%"%(i,am[i]/ch,(am[i]-self.mono)/ch,(am[i]-am[i-1])/ch,100*p/maxprob)
        print >>h, "  A %10.4f"%((aam/totp)/ch,)

    def clusterIntensities(self):
        t=defaultdict(float)
        maxNC13 = 0
        totalprob = 0.0
        for (p,m,i) in self.table:
            nC13 = int(round(m-self.mono))
            t[nC13] += p
            totalprob += p
            if nC13 > maxNC13:
                maxNC13 = nC13
        for k in t.keys():
            t[k] /= totalprob
        return [t[k] for k in range(maxNC13+1)]
        
    def mass(self,m,n):
        return sum([ mi*ni for (mi,ni) in zip(m,n) ])
    
    def multinomial(self,p,n,N):
        return math.exp(self.gammaln(N+1) - sum([self.gammaln(n[i]+1) for i in xrange(0,len(n))]) + \
                        sum([n[i]*math.log(p[i]) for i in xrange(0,len(n)) if p[i] > 0 and n[i]>0]))

    def lnmultinomial(self,p,n,N):
        # print p,n,N
        return self.gammaln(N+1) - sum([self.gammaln(n[i]+1) for i in xrange(0,len(n))]) + \
               sum([n[i]*math.log(p[i]) for i in xrange(0,len(n)) if p[i] > 0])

    def gammaln(self,xx):
        cof = (76.18009172947146,-86.50532032941677,24.01409824083091,\
               -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5); 
        y=xx; x=xx;
        tmp = x+5.5;
        tmp -= (x+0.5)*math.log(tmp);
        series=1.000000000190015;
        for j in xrange(0,6):
            y += 1
            series += cof[j]/y
        return -tmp+math.log(2.5066282746310005*series/x);

