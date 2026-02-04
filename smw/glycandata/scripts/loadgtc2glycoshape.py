#!/bin/env python3.12

import sys, time
from collections import defaultdict

from getwiki import GlycanData, Glycan
w = GlycanData()

import findpygly
from pygly.GlycanResource import GlycoShapeSourceFile

gsdb = GlycoShapeSourceFile()

glycoshape = defaultdict(set)

for gsdbid,gtc,dummy,dummy1 in gsdb.allgtc():
    glycoshape[gtc].add(gsdbid)

for m in w.iterglycan():
    start=time.time()
    acc = m.get('accession')
    if len(glycoshape[acc]) > 0:
        m.set_annotation(value=list(glycoshape[acc]),property="GlycoShape",source="GlycoShape",type="CrossReference")
    else:
        m.delete_annotations(source="GlycoShape",property="GlycoShape",type="CrossReference")
    if w.put(m):                                                                                                                     
        print("%s updated in %.2f sec"%(acc,time.time()-start,),file=sys.stderr)                                                      
    else:                                                                                                                             
        print("%s checked in %.2f sec"%(acc,time.time()-start,),file=sys.stderr)
