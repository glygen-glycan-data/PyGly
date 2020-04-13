#!/bin/env python27
import sys, re
from collections import defaultdict
for l in sys.stdin:
    sl = re.split(r'(\d+)',l.strip())
    comp = defaultdict(int)
    for i in range(0,len(sl)-1,2):
	comp[sl[i]] = int(sl[i+1])
    outstr = ""
    for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P','KDN','HexA']:
	if comp[k] > 0:
	    outstr += k + str(comp[k])
    print outstr
