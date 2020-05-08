#!/bin/env python27
import sys, re
from collections import defaultdict

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f','--format',type='choice',choices=['UCKBCOMP','SHORTUCKB'],default='UCKBCOMP',
                  dest='format')
parser.add_option('-p','--position',type='int',default=None,dest='pos')
parser.add_option('-d','--delim',type='string',default=None,dest='delim')
parser.add_option('-r','--require',type='string',default=None,dest='require')
parser.add_option('-a','--all',action='store_true',default=False,dest='all')

opt,arg = parser.parse_args()
if opt.pos != None:
    opt.pos -= 1
if opt.require:
    opt.require = re.compile(opt.require)

if opt.delim == None:
    outdelim = "\t"
    indelim = None
else:
    outdelim = opt.delim
    indelim = opt.delim

keys = ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P','KDN','HexA']

for l in sys.stdin:
    l = l.strip()
    if opt.pos != None:
	sl = l.split(indelim)
	if len(sl) <= opt.pos:
	    if opt.all:
	        print outdelim.join(sl)
	    continue
	elif opt.require and not opt.require.search(sl[opt.pos]):
	    if opt.all:
	        print outdelim.join(sl)
	    continue
	compstr = sl[opt.pos]
    else:
	if opt.require and not opt.require.search(l):
	    if opt.all:
		print l
	    continue
	compstr = l
    prefix = ""
    if compstr.startswith('comp_'):
	prefix = 'comp_'
	compstr = compstr[5:]
    sl = re.split(r'(\d+)',compstr)
    comp = defaultdict(int)
    for i in range(0,len(sl)-1,2):
	comp[sl[i]] = int(sl[i+1])
    for key in comp:
        if comp[key] > 0 and key not in keys:
	    print >>sys.stderr, "Can't normalize",l
	    sys.exit(1)
    outstr = ""
    for k in ['HexNAc','Hex','dHex','NeuAc','NeuGc','Pent','S','P','KDN','HexA']:
	if comp[k] > 0 or opt.format == 'UCKBCOMP':
	    outstr += (k + str(comp[k]))
    if not opt.all and outstr == compstr:
	continue
    if prefix:
	outstr = prefix + outstr;
    if opt.pos != None:
	sl = l.split(indelim)
	sl[opt.pos] = outstr
	print outdelim.join(sl)
    else:
        print outstr
