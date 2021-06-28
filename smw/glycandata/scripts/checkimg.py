#!/bin/env python2

import sys, os, os.path
import findpygly
from pygly.GlyTouCan import GlyTouCan
from pygly.alignment import GlycanImageEqual

import csv
from collections import defaultdict

imagecrcfile = sys.argv[1]
imagecrc = defaultdict(lambda:defaultdict(set))
gtc2crc = dict()
for r in csv.DictReader(open(imagecrcfile),dialect='excel-tab'):
    key = (int(r['Image-Size']),r['Image-CRC'])
    style = r['Image-Style']
    gtcacc = r['GlyTouCanAccession']
    imagecrc[key][style].add(gtcacc)
    gtc2crc[(style,gtcacc)] = key

tocheck = list()
for key in imagecrc:
    for style in ('extended','normal','compact'):
        if len(imagecrc[key][style]) > 1:
	    tocheck.append((style,sorted(imagecrc[key][style])))

for style,gtcacc in gtc2crc:
    key = gtc2crc[(style,gtcacc)]
    if style == 'extended':
	key1 = gtc2crc[('normal',gtcacc)]
	if not imagecrc[key]['extended'] <= imagecrc[key1]['normal']:
	    print "%s images for %s are identical, but their normal images are not"%(style,", ".join(sorted(imagecrc[key]['extended'])))
	key2 = gtc2crc[('compact',gtcacc)]
	if not imagecrc[key]['extended'] <= imagecrc[key2]['compact']:
	    print "%s images for %s are identical, but their compact images are not"%(style,", ".join(sorted(imagecrc[key]['extended'])))
    elif style == 'normal':
	key2 = gtc2crc[('compact',gtcacc)]
	if not imagecrc[key]['normal'] == imagecrc[key2]['compact']:
	    print "%s images for %s are identical, but their compact images are not"%(style,", ".join(sorted(imagecrc[key]['normal'])))
    elif style == 'compact':
	key1 = gtc2crc[('normal',gtcacc)]
	if not imagecrc[key]['compact'] == imagecrc[key1]['normal']:
	    print "%s images for %s are identical, but their normal images are not"%(style,", ".join(sorted(imagecrc[key]['compact'])))

# imagetopocmp = GlycanImageEqual()
# 
# gtc = GlyTouCan(usecache=True)
# 
# for style,gtcaccs in sorted(tocheck):
#     glys = map(gtc.getGlycan,gtcaccs)
#     for i in range(0,len(glys)):
# 	for j in range(i+1,len(glys)):
# 	    if not imagetopocmp.eq(glys[i],glys[j]):
#		size,crc = gtc2crc[(style,gtcaccs[i])]
# 	        print "%s images for %s and %s are identical (%s,%s), but their (image-relevant) topologies are not equal"%(style,gtcaccs[i],gtcaccs[j],size,crc)

