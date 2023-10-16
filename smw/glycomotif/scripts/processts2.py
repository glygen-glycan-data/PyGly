#!/bin/env python2

import csv, sys
from collections import defaultdict

from getwiki import GlycoMotifWiki

inputfile = sys.argv[1]

def assigncat(vals,total,specgrp,highlow):

    detvals = [ x for x in vals if x > 1 ]
    numdet = len(detvals)
        
    if numdet == 0:
        return "Not detected",None
    
    mindet = min(detvals)
    maxdet = max(detvals)
    
    if vals[0] >= 50*max(vals[1],0.1):
        return "Highly enriched",1
    if vals[0] >= 5*max(vals[1],0.1):
        return "Moderately enriched",1
    for i in range(2,min(specgrp+1,numdet+1)):
        grpmean = sum(vals[:i])/i
        if grpmean > 5*max(vals[i],0.1):
        # if vals[i-1] >= 5*max(vals[i],0.1):
            return "Group enriched",i
    
    if numdet == total:
        if mindet >= highlow:
            return "Ubiquitous high",(numdet,len([ x for x in vals if x > highlow ]))
        else:
            return "Ubiquitous low",(numdet,len([ x for x in vals if x > highlow ]))
    if 0 < numdet < total:
        if maxdet > highlow:
            return "Mixed high",(numdet,len([ x for x in vals if x > highlow ]))
        else:
            return "Mixed low",(numdet,len([ x for x in vals if x > highlow ]))
    
    raise "No assignment made"

def process_pairs(g,gn,pairs,total,grpspec,highlow):
    pairs = sorted(pairs,reverse=True)
    values = [ p[0] for p in pairs ]
    groups = [ p[1] for p in pairs ]
    cat,ind = assigncat(values,total,grpspec,highlow)
    spec = "-"
    if 'enriched' in cat:
        spec = ",".join(groups[:ind])
    elif 'high' in cat or 'low' in cat:
        spec = "detected=%.2f%%,high=%.2f%%"%(100*ind[0]/total,100*ind[1]/total)
    print("\t".join(map(str,[g,gn,cat,spec])))

fieldnames = ['Tissue','Cell type']
rows = csv.DictReader(open(inputfile),dialect='excel-tab')
allgroups = set()
fieldname = None
for r in rows:
    gn = r['Gene name']
    if not fieldname:
        for fn in fieldnames:
            if fn in r:
                fieldname = fn
                break
    allgroups.add(r[fieldname])

total = len(allgroups)
grpspec = 5
highlow = 10

rows = csv.DictReader(open(inputfile),dialect='excel-tab')
lastgene = None
for r in rows:
    g = r['Gene']
    gn = r['Gene name']
    if gn.startswith("ENS"):
        continue
    if g != lastgene:
        if lastgene != None:
            values = [ (pairs[k],k) for k in allgroups ]
            process_pairs(lastgene,lastgn,values,len(allgroups),grpspec,highlow)
        pairs = defaultdict(float)
        lastgene = g
        lastgn = gn
    pairs[r[fieldname]] = float(r['nTPM'])
if lastgene != None:
    values = [ (pairs[k],k) for k in allgroups ]
    process_pairs(lastgene,lastgn,values,len(allgroups),grpspec,highlow)
