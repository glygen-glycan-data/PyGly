#!/bin/env python3

import sys
import findpygly
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoPrefetch, GlyTouCanNoCache, GlyTouCan
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from pygly.manipulation import WURCSArchetype, Archetype
from pygly.alignment import GlycanEqual

gtc = GlyTouCanNoPrefetch(verbose=False)
wp = WURCS20Format()
warch = WURCSArchetype()
arch = Archetype()
geq = GlycanEqual()

# badacc = """
#     G63937EY G69594ET G04765HR G07417TQ G11392MR G27289CV G31198WT G49548QF G52742DF G62779TI G70119FR G72112TE
#     G79266BG G84928GN G86111EQ
# """

# badacc = """
#     G19924OM G72462KE
# """

badacc = ""

for line in sys.stdin:
    acc,aacc = line.split()
    if acc in badacc.split():
        continue
    print acc,aacc 
    wurcs = gtc.getseq(acc,format='wurcs')
    awurcs = gtc.getseq(aacc,format='wurcs')
    try:
        gly = wp.toGlycan(wurcs)
    except GlycanParseError:
        gly = None
    try:
        awgly = wp.toGlycan(awurcs)
    except GlycanParseError:
        awgly = None
    if not gly or not awgly:
        continue
    if gly.repeated():
        continue
    agly = arch(gly)
    # if not geq.eq(agly,awgly):
    #     continue
    archwurcs = warch(wurcs)
    print(wurcs)
    print(archwurcs)
    print(awurcs)
    mono0 = wurcs.split(']')[0].split('[')[1]
    mono1 = archwurcs.split(']')[0].split('[')[1]
    mono2 = awurcs.split(']')[0].split('[')[1]
    assert mono1 == mono2, " ".join([mono0, mono1, mono2])
    wagly = wp.toGlycan(archwurcs)
    if not geq.eq(awgly,agly) or not geq.eq(wagly,agly):
        print("archetype:")
        print(awgly.glycoct())
        print("manip arch:")
        print(agly.glycoct())
        print("wurcs manip arch:")
        print(wagly.glycoct())
    assert(geq.eq(awgly,agly))
    assert(geq.eq(wagly,agly))
