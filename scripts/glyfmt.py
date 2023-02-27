#!/bin/env python2
import sys, os.path
import findpygly
from pygly.GlycanFormatter import *
from pygly.GlycanBuilderSVGParser import *

clsinst = eval("%s()"%sys.argv[1])
if len(sys.argv) > 2:
    files = sys.argv[2:]
else:
    files = sys.stdin.read().splitlines()
bad = 0
for f in files:
    seq = open(f).read()
    try:
        g = clsinst.toGlycan(seq)
        print("+++", os.path.split(f)[1])
        print(GlycoCTFormat().toXML(g))
        for m in g.all_nodes(undet_subst=True):
            print(m)
        if not g.repeated():
            print(g.underivitized_molecular_weight())
        else:
            print(g.underivitized_molecular_weight(repeat_times=1))
    except GlycanParseError as e:
        print("!!!", os.path.split(f)[1], e)
        bad += 1
print("Failed: %d/%d"%(bad,len(files)))
