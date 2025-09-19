#!/bin/env python3
import sys, os.path
import findpygly
from pygly.GlycanFormatter import *
from pygly.GlycanBuilderSVGParser import *
from pygly.CompositionFormatter import *
from pygly.GlycanResource import GlyTouCanNoPrefetch as GlyTouCan
from pygly.GlycanMultiParser import GlycanMultiParser

if sys.argv[1] not in ("glycoct","wurcs","svg","comp","iupac","iupac2","iupac3","gtc","auto"):
    print("Parser should be one of: glycoct, wurcs, svg, comp, iupac, iupac2, iupac3, auto, gtc.")
    exit(1)
gtc = None
clsname = None
if sys.argv[1] == "glycoct":
    clsname = "GlycoCTFormat"
elif sys.argv[1] in ("wurcs","gtc"):
    clsname = "WURCS20Format"
    if sys.argv[1] == 'gtc':
        gtc = GlyTouCan()
elif sys.argv[1] == "svg":
    clsname = "GlycanBuilderSVG"
elif sys.argv[1] == "comp":
    clsname = "CompositionFormat"
elif sys.argv[1] == "iupac":
    clsname = "IUPACLinearFormat"
elif sys.argv[1] == "iupac2":
    clsname = "IUPACParserExtended1"
elif sys.argv[1] == "iupac3":
    clsname = "IUPACParserGlyTouCanExtended"
elif sys.argv[1] == "auto":
    clsname = "GlycanMultiParser"

iupacseq = IUPACLinearFormat()

clsinst = eval("%s()"%(clsname,))
if len(sys.argv) > 2:
    files = sys.argv[2:]
else:
    files = sys.stdin.read().splitlines()
bad = 0
for f in files:
    if os.path.exists(f):
        seq = open(f).read()
    elif re.search(r'^G[0-9]{5}[A-Z]{2}$',f):
        seq = gtc.getseq(f,'wurcs')
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
        print(iupacseq.toStr(g))
        # print(clsinst.toStr(g))
        if clsname != "GlycanMultiParser":
            print("Parser:",clsname)
        else:
            print("Parser:",clsinst.parser_name())
    except GlycanParseError as e:
        # traceback.print_exc()
        print("!!!", os.path.split(f)[1], str(e))
        bad += 1
    except Exception as e:
        print("!!!", os.path.split(f)[1], str(e))
        bad += 1
print("Failed: %d/%d"%(bad,len(files)))
