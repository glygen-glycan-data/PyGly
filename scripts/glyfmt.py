#!/bin/env python2
import sys, os.path
import findpygly
from pygly.GlycanFormatter import *
from pygly.GlycanBuilderSVGParser import *

if sys.argv[1] not in ("glycoct","wurcs","svg"):
    print("Parser should be one of: glycoct, wurcs, svg.")
    exit(1)
if sys.argv[1] == "glycoct":
    clsname = "GlycoCTFormat"
elif sys.argv[1] == "wurcs":
    clsname = "WURCS20Format"
elif sys.argv[1] == "svg":
    clsname = "GlycanBuilderSVG"

clsinst = eval("%s()"%(clsname,))
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
        # print(GlycoCTFormat().toXML(g))
        print(GlycoCTFormat().toStr(g))
        for m in g.all_nodes(undet_subst=True):
            print(m)
        if not g.repeated():
            print(g.underivitized_molecular_weight())
        else:
            print(g.underivitized_molecular_weight(repeat_times=1))
    except GlycanParseError as e:
        print("!!!", os.path.split(f)[1], e)
        print(e)
        bad += 1
    except KeyError:
        pass
print("Failed: %d/%d"%(bad,len(files)))