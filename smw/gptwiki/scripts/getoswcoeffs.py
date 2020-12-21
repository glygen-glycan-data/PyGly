#!/bin/env python27

import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import getwiki

import csv,sys
import matplotlib.pyplot as plt
from analysis.regression import IteratedRobustLinearRegression as IRLR
from analysis.regression import SimpleLinearRegression as SLR
from analysis.fdr import CombinedAnalysisFDR

from operator import itemgetter

if len(sys.argv) < 4:
    print >>sys.stderr, 'please enter the table.tsv, library file, #decoys !'
    sys.exit(1)

fdr = CombinedAnalysisFDR(ndecoys=int(sys.argv[3]),
                          rankingkey="main_var_xx_swath_prelim_score")
fdr.add_ids(sys.argv[1])
tgs = None
for qv,sc,trg in sorted(map(itemgetter(1,0,2),fdr.allqvalues())):
    if qv > 0.05:
        break
    # print >>sys.stderr, (len(tgs) if tgs != None else 0), sc, qv, trg
    if tgs != None and len(tgs) >= 10 and qv > 0.01:
	break
    tgs = dict()
    for row in fdr.filter_ids(sys.argv[1],score=sc):
        tgid = row['peptide_group_label']
        rt = float(row['RT'])/60.0
        tgs[tgid] = rt

if not tgs or len(tgs) < 10:
    print >>sys.stderr, "Not enough good quality transition groups"
    sys.exit(1)

# print >>sys.stderr, (len(tgs) if tgs != None else 0), sc, qv, trg

points = []
seentg = set()
for row in csv.DictReader(open(sys.argv[2]), dialect='excel-tab'):
    tgid = row['transition_group_id']
    if tgid in seentg:
        continue
    seentg.add(tgid)
    if tgid in tgs:
        assay_rt = float(row['Tr_recalibrated'])
        points.append((assay_rt,tgs[tgid]))

def print_params(params,msg=None):
    if msg:
        print >>sys.stderr, msg,
    print >>sys.stderr, round(params['slope'],3),round(params['intercept'],3),
    print >>sys.stderr, round(1.0/params['slope'],3),round(-params['intercept']/params['slope'],3),
   
    if 'removed' in params:
        print >>sys.stderr, params['retained'],params['removed'],round(params['r'],3)
    else:
	print >>sys.stderr, params['npoints'],0,round(params['r'],3)

try:
    reg0 = SLR()
    params0 = reg0.fit(points)
    print_params(params0," SLR:")
except RuntimeError:
    pass

try:
    reg1 = IRLR(min_points=5,max_rvalue=0.99,max_removed=1e+20)
    params1 = reg1.fit(points)
    print_params(params1,"IRLR:")
except RuntimeError:
    print >>sys.stderr, "LC Calibration failed"
    sys.exit(1)

params1["y0"] = -params1['intercept']/params1['slope']
params1["y1"] = 1.0/(params1['slope']*60.0)-params1['intercept']/params1['slope']

# pairs = sorted([(60*p[1],p[0]) for p in params1['retained_points']])
pairs = [(0,params1["y0"]),(1,params1["y1"])]

npairs = len(pairs)
pairs = '\n'.join(['      <Pair from="%s" to="%s"/>'%(p[0],p[1]) for p in pairs])

print (("""
<?xml version="1.0" encoding="UTF-8"?>
<TrafoXML version="1.0" xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/TrafoXML_1_0.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <Transformation name="linear" nrtslope="%s" nrtintercept="%s">
    <Pairs count="%s">
%s
    </Pairs>
  </Transformation>
</TrafoXML>
"""%(params1['slope'],params1['intercept'],npairs,pairs)).strip())

