#!/bin/env python3.12

import sys, os, os.path, time
import findpygly
from pygly.GlycanResource import GlyTouCan
from hashlib import md5
import urllib

import findpygly
from getwiki import GlycanData
w = GlycanData()

notation = 'snfg'
style = 'extended'
format = 'svg'

if len(sys.argv) > 1:
    notation = sys.argv[1]
if len(sys.argv) > 2:
    style = sys.argv[2]
if len(sys.argv) > 3:
    format = sys.argv[3]

thepath = 'images/%s/%s/%s'%(notation,style,format)
if len(sys.argv) > 4:
    thepath = sys.argv[4]

badmd5sum = """
f6ad129737d637b63264ae827edd9887
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

for g in w.iterglycan():
    gtcacc = g.get('accession')
    imgfn = "%s/%s.%s"%(thepath,gtcacc,format)
    if os.path.exists(imgfn):
	continue
    imgstr = urllib.urlopen("https://glymage.glyomics.org/image/%s/extended/%s.%s"%(notation,gtcacc,format)).read()
    if not goodimg(format,imgstr):
	continue
    print("writing:",imgfn)
    wh = open(imgfn,'w')
    wh.write(imgstr)
    wh.close()
    time.sleep(1)
