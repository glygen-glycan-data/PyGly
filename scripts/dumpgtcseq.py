#!/bin/env python2

import sys, glob
import findpygly
from pygly.GlycanResource import GlyTouCan, GlyCosmos

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

    if filename in allfn:
        allfn.remove(filename)  
        continue

    seq = None                                                                                                               
    if seqtype == 'wurcs':
        seq = gtc.getseq(acc,'wurcs')

    if seqtype == 'glycoct':
        seq = gtc.getseq(acc,'glycoct')
                                                                                                                             
    if seqtype == 'genglycoct':
        seq = gtc.glycoct(acc)

    if not seq:                                                                                                              
        continue                                                                                                             
                                                                                                                             
    print >>sys.stderr, "Write:",acc                                                                                         
    wh = open(filename,'w')                                                                                                  
    wh.write(seq)                                                                                                            
    wh.close()

for fn in allfn:
    print >>sys.stderr, "Remove:",acc
    os.unlink(fn)
