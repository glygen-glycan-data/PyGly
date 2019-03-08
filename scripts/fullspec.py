#!/usr/bin/env python27

import sys
import findpygly
from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

for l in sys.stdin:
    acc = l.strip()
    if gtc.getGlycan(acc).fully_determined():
	seq = gtc.getseq(acc,'glycoct')
	if seq and not os.path.exists('%s.txt'%(acc,)):
	    wh = open('%s.txt'%(acc,),'wb')
	    wh.write(seq)
	    wh.close()
