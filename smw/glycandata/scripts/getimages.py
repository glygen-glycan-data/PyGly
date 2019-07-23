#!/bin/env python27

import sys, os, os.path
import findpygly
from pygly.GlyTouCan import GlyTouCan

def accessions(args):
    if len(args) == 0:
        for it in sys.stdin:
            yield it.strip()
    else:
        for fn in args:
            for it in open(fn):
                yield it.strip()

notation = sys.argv[1]
style = sys.argv[2]

gtc = GlyTouCan(usecache=True)
for i,gtcacc in enumerate(accessions(sys.argv[3:])):
    imgfn = "%s.png"%(gtcacc,)
    if os.path.exists(imgfn):
	continue
    imgstr,width,height = gtc.getimage(gtcacc,style=style,notation=notation,avoidcache=True,trials=3)
    if imgstr and width and height:
	print "writing:",imgfn
        wh = open(imgfn,'w')
        wh.write(imgstr)
        wh.close()
