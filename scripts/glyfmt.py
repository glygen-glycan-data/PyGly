#!/bin/env python2
import sys, os.path
import findpygly
from pygly.GlycanFormatter import *
from pygly.GlycanBuilderSVGParser import *
from pygly.CompositionFormatter import *

if sys.argv[1] not in ("glycoct","wurcs","svg","comp","iupac","iupac2"):
    print("Parser should be one of: glycoct, wurcs, svg, comp, iupac, iupac2.")
    exit(1)
if sys.argv[1] == "glycoct":
    clsname = "GlycoCTFormat"
elif sys.argv[1] == "wurcs":
    clsname = "WURCS20Format"
elif sys.argv[1] == "svg":
    clsname = "GlycanBuilderSVG"
elif sys.argv[1] == "comp":
    clsname = "CompositionFormat"
elif sys.argv[1] == "iupac":
    clsname = "IUPACLinearFormat"
elif sys.argv[1] == "iupac2":
    clsname = "IUPACParserExtended1"

clsinst = eval("%s()"%(clsname,))
if len(sys.argv) > 2:
    files = sys.argv[2:]
else:
    files = sys.stdin.read().splitlines()
bad = 0
for f in files:
    if os.path.exists(f):
        seq = open(f).read()
    else:
        seq = f.strip()
    try:
        g = clsinst.toGlycan(seq)
        print("+++", os.path.split(f)[1])
        # print(GlycoCTFormat().toXML(g))
        # print(GlycoCTFormat().toStr(g))
        for m in g.all_nodes(undet_subst=True):
            print(m)
        if not g.repeated():
            print(g.underivitized_molecular_weight())
        else:
            print(g.underivitized_molecular_weight(repeat_times=1))
        print(g.glycoct())
    except GlycanParseError as e:
        print("!!!", os.path.split(f)[1], str(e))
        bad += 1
    except Exception as e:
        print("!!!", os.path.split(f)[1], str(e))
        bad += 1
print("Failed: %d/%d"%(bad,len(files)))
