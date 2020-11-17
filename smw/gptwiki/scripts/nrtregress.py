#!/bin/env python27

import matplotlib
matplotlib.use('Agg')

import getwiki
from analysis.regression import SimpleLinearRegression as SLR
from analysis.regression import IteratedRobustLinearRegression as IRLR
from analysis.chromatogram import GaussianPeakFit as GPF
from operator import itemgetter
import csv, sys, os, os.path
from collections import defaultdict
from pylab import *

w = getwiki.GPTWiki()

slr = SLR()
tgreg = IRLR(max_rvalue=0.99,max_removed=50)

spec = w.get(sys.argv[1])
params0 = {}
try:
    params0['slope'] = float(spec.get('nrtslope'))
    params0['intercept'] = float(spec.get('nrtintercept'))
except TypeError:
    pass

tgpoints = []
for tg in w.iterspectgs(sys.argv[1]):
    tgprt = tg.get('prt')
    if not tgprt:
	continue
    pepid = tg.get('peptide')
    pep = w.get(pepid)
    pepnrt = pep.get('nrt')
    if not pepnrt:
        continue
    tgprt = float(tgprt)
    pepnrt = float(pepnrt)
    tgpoints.append((pepnrt,tgprt))

tgparams0 = slr.fit(tgpoints)
tgparams = tgreg.fit(tgpoints)
print "ORIGLR:",round(params0['slope'],3),round(params0['intercept'],3),
print round(1.0/params0['slope'],3),round(-params0['intercept']/params0['slope'],3)
print "TG SLR:",round(tgparams0['slope'],3),round(tgparams0['intercept'],3),
print round(1.0/tgparams0['slope'],3),round(-tgparams0['intercept']/tgparams0['slope'],3)
print "TG RLR:",round(tgparams['slope'],3),round(tgparams['intercept'],3),
print round(1.0/tgparams['slope'],3),round(-tgparams['intercept']/tgparams['slope'],3),
print tgparams['removed'],round(tgparams['r'],3)

minx = 1e+20
minx = min(minx,*map(itemgetter(0),tgpoints))
maxx = -1e+20
maxx = max(maxx,*map(itemgetter(0),tgpoints))
extremes = [minx,maxx]
plot([minx,maxx],slr.y(params0,extremes),'-.',color="cornflowerblue")
plot([minx,maxx],slr.y(tgparams0,extremes),':',color="cornflowerblue")
plot([minx,maxx],tgreg.y(tgparams,extremes),'--',color="cornflowerblue")
plot(map(itemgetter(0),tgpoints),map(itemgetter(1),tgpoints),'b.')
xlabel("NRT")
ylabel("Obs RT")
savefig(sys.argv[1]+".nrt.png")
# show()
