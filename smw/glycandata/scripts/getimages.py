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

badmd5sum = """
66ceaa814dd7355a3a8d111af7c3075a
e7183de88ac19ecc75544e939a2d056e
""".split()

def goodimg(format,imgstr):
    if not imgstr:
        return False
    if format == "png" and not imgstr.startswith('\x89PNG\x0d\x0a\x1a\x0a'):
	return False
    if format == "svg" and '<svg ' not in imgstr:
	return False
    if md5(imgstr).hexdigest().lower() in badmd5sum:
	return False
    return True

gtc = GlyTouCan(usecache=False)
for gtcacc in accessions(sys.argv[4:]):
    imgfn = "%s.%s"%(gtcacc,format)
    if os.path.exists(imgfn):
	continue
    imgstr = None
    if not imgstr:
        imgstr = gtc.getimage(gtcacc,style=style,notation=notation,format=format)
	if not goodimg(format,imgstr):
	    imgstr = None
    if not imgstr:
	if style == "extended":
	    try:
	        imgstr = urllib.urlopen("https://image.glycosmos.org/%s/%s/%s"%(notation,format,gtcacc,)).read()
	    except IOError:
		pass
	    if not goodimg(format,imgstr):
	        imgstr = None
    if not imgstr:
	if style == "extended" and notation == "snfg" and format == "png":
	    try:
	        imgstr = urllib.urlopen("https://api.glycosmos.org/wurcs2image/experimental/%s/binary/%s"%(format,gtcacc)).read()
	    except IOError:
		pass
	    if not goodimg(format,imgstr):
	        imgstr = None
    if not imgstr:
	if style == "extended" and notation == "snfg" and format == "png":
            wurcs = gtc.getseq(gtcacc,'wurcs')
	    try:
	        imgstr = urllib.urlopen("https://api.glycosmos.org/wurcs2image/experimental/%s/binary/%s"%(format,wurcs)).read()
	    except IOError:
		pass
	    if not goodimg(format,imgstr):
	        imgstr = None
    if not imgstr:
	continue
    print "writing:",imgfn
    wh = open(imgfn,'w')
    wh.write(imgstr)
    wh.close()
