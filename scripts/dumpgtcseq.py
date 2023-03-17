#!/bin/env python2

import sys, glob, os
import findpygly
from pygly.GlycanResource import GlyTouCan, GlyCosmos
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError

seqtype = sys.argv[1].split('/')[-1]
assert(seqtype in ('wurcs','glycoct','genglycoct'))

allfn = set(glob.glob(sys.argv[1]+"/G*.txt"))

gco = GlyCosmos(usecache=False)
archived = set(map(lambda d: d['accession'],gco.archived()))
gtc = GlyTouCan(verbose=False,usecache=False)

for acc in gtc.allaccessions():

    if acc in archived:
	continue

    filename = sys.argv[1] + "/" + acc + ".txt"

    seq = None                                                                                                               
    if seqtype == 'wurcs':
        seq = gtc.getseq(acc,'wurcs')
    elif seqtype == 'glycoct':
        seq = gtc.getseq(acc,'glycoct')
    elif seqtype == 'genglycoct':
        seq = gtc.glycoct(acc)

    if not seq:                                                                                                              
        continue                                                                                                             

    if filename not in allfn:
        print >>sys.stderr, "Write:",acc
        wh = open(filename,'w')
        wh.write(seq) 
        wh.close()
    else:
        allfn.remove(filename)

for fn in allfn:
    print >>sys.stderr, "Remove:",fn
    os.unlink(fn)
