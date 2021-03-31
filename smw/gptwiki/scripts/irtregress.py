#!/bin/env python27

import getwiki
from analysis.regression import SimpleLinearRegression as SLR
from analysis.regression import IteratedRobustLinearRegression as IRLR
from analysis.chromatogram import GaussianPeakFit as GPF
from operator import itemgetter
import csv, sys, os, os.path
from collections import defaultdict
from pylab import *

# w = getwiki.GPTWiki()

irtpeps = dict()
for r in csv.DictReader(open(sys.argv[1]),dialect='excel-tab'):
    if r["Name"] not in ("irt-a",):
        irtpeps[r["Name"]] = float(r["NRT"])

pf = GPF()
slr = SLR()
irtreg = IRLR(max_rvalue=0.99)
tgreg = IRLR(max_rvalue=0.99,max_removed=50)

xics = defaultdict(list)
for r in csv.DictReader(open(sys.argv[2]),dialect='excel-tab'):
    for k in r:
        if k == "rt":
            continue
        name = k.split('(')[1][:-1]
        xics[name].append((float(r["rt"]),float(r[k])))

if len(xics) == 0:
    sys.exit(1)

irtpoints = []
for k in xics:
    if k not in irtpeps:
	continue
    try:
        xs,xe,pa = pf.fit(xics[k])
    except (RuntimeError,FloatingPointError,ValueError):
	continue
    irtpoints.append((irtpeps[k],pa[0]))

spec = os.path.split(sys.argv[2])[1].split('.')[0]
try:
    irtparams0 = slr.fit(irtpoints)
    irtparams = irtreg.fit(irtpoints)
except RuntimeError:
    # print "Can't fit regression lines"
    sys.exit(1)

print "XICS:",os.path.split(sys.argv[2])[1]
print "iRT Peaks",len(irtpoints)
print "iRT SLR:",round(irtparams0['slope'],3),round(irtparams0['intercept'],3),
print round(1.0/irtparams0['slope'],3),round(-irtparams0['intercept']/irtparams0['slope'],3)
print "iRT RLR:",round(irtparams['slope'],3),round(irtparams['intercept'],3),
print round(1.0/irtparams['slope'],3),round(-irtparams['intercept']/irtparams['slope'],3),
print irtparams['removed'],round(irtparams['r'],3)
print "|nrtintercept=%(intercept)f\n|nrtslope=%(slope)f"%irtparams
print

sys.exit(0)

tgpoints = []

for tg in w.itertgs(spectra=spec):
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

if len(tgpoints) == 0:
    sys.exit(0)

tgparams0 = slr.fit(tgpoints)
tgparams = tgreg.fit(tgpoints)
print "TG SLR:",round(tgparams0['slope'],3),round(tgparams0['intercept'],3),
print round(1.0/tgparams0['slope'],3),round(-tgparams0['intercept']/tgparams0['slope'],3)
print "TG RLR:",round(tgparams['slope'],3),round(tgparams['intercept'],3),
print round(1.0/tgparams['slope'],3),round(-tgparams['intercept']/tgparams['slope'],3),
print tgparams['removed'],round(tgparams['r'],3)

minx = 1e+20
minx = min(minx,*map(itemgetter(0),irtpoints))
minx = min(minx,*map(itemgetter(0),tgpoints))
maxx = -1e+20
maxx = max(maxx,*map(itemgetter(0),irtpoints))
maxx = max(maxx,*map(itemgetter(0),tgpoints))
extremes = [minx,maxx]
plot([minx,maxx],slr.y(tgparams0,extremes),':',color="cornflowerblue")
plot([minx,maxx],slr.y(irtparams0,extremes),':',color="salmon")
plot([minx,maxx],tgreg.y(tgparams,extremes),'--',color="cornflowerblue")
plot([minx,maxx],irtreg.y(irtparams,extremes),'--',color="salmon")
plot(map(itemgetter(0),tgpoints),map(itemgetter(1),tgpoints),'b.')
plot(map(itemgetter(0),irtpoints),map(itemgetter(1),irtpoints),'r.')
xlabel("NRT")
ylabel("Obs RT")
savefig(sys.argv[2]+".irt.png")
# show()
