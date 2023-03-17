#!/bin/env python

import sys, re
import findpygly
from pygly.GlycanFormatter import IUPACLinearFormat, GlycanParseError, GlycoCTFormat

ind = int(sys.argv[1])
fmt = IUPACLinearFormat()
gct = GlycoCTFormat()
for l in sys.stdin:
    sl = l.split()
    seq = sl[ind]
    seq = re.sub(r'\(([?ab][?\d]-[?\d])\)','\g<1>',seq)
    seq = seq.replace('[','(')
    seq = seq.replace(']',')')
    try:
        g = fmt.toGlycan(seq)
	# 	l = l.strip() + "\n"+gct.toStr(g) + "\n"
	l = l.strip() + "\t"+fmt.toStr(g) + "\n"
    except (GlycanParseError,KeyError):
        # sys.stdout.write(l)
        continue #raise
    sys.stdout.write(l)
