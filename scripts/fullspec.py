#!/usr/bin/env python27

import sys, os, os.path
import findpygly
from pygly.GlyTouCan import GlyTouCan
from pygly.GlycanFormatter import GlycoCTFormat

glycoct_format = GlycoCTFormat()
gtc = GlyTouCan()

for l in sys.stdin:
    acc = l.strip()
    g =  gtc.getGlycan(acc)
    if g and g.fully_determined():
	print acc,True
	if not os.path.exists('%s.txt'%(acc,)):
	    seq = gtc.getseq(acc,'glycoct')
	    if not seq:
		try:
		    seq = glycoct_format.toStr(g)
		except:
		    pass
	    if seq:
	        wh = open('%s.txt'%(acc,),'w')
	        wh.write(seq.strip()+'\n')
	        wh.close()
    else:
	print acc,(False if g else None)
