#!/bin/env python27

import sys, os, os.path
import findpygly
from pygly.GlycanResource import GlyTouCan
from hashlib import md5
import urllib

def accessions(args):
    if len(args) == 0:
        for it in sys.stdin:
            yield it.strip()
    else:
        for fn in args:
            for it in open(fn):
                yield it.strip()

notation = 'snfg'
style = 'extended'
format = 'svg'

if len(sys.argv) > 1:
    notation = sys.argv[1]
if len(sys.argv) > 2:
    style = sys.argv[2]
if len(sys.argv) > 3:
    format = sys.argv[3]

gtc = GlyTouCan(usecache=False)
for gtcacc in accessions(sys.argv[4:]):
    imgfn = "%s.%s"%(gtcacc,format)
    if os.path.exists(imgfn):
	continue
    imgstr = gtc.getimage(gtcacc,style=style,notation=notation,format=format)
    if not imgstr:
	if style == "extended":
	    try:
	        imgstr = urllib.urlopen("https://image.glycosmos.org/%s/%s/%s"%(notation,format,gtcacc,)).read()
	    except IOError:
		pass
	continue
    if md5(imgstr).hexdigest().lower() == "e7183de88ac19ecc75544e939a2d056e":
	continue
    print "writing:",imgfn
    wh = open(imgfn,'w')
    wh.write(imgstr)
    wh.close()
