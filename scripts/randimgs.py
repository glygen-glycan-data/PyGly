#!/bin/env python2
from __future__ import print_function

import sys, os, random
import findpygly
from pygly.GlycanImage import GlycanImage
from pygly.GlycanResource.GlyTouCan import GlyTouCanNoCache, GlyTouCan


accs = list(open(sys.argv[1]).read().split())

batch = 10
iterations = 100
scale_options = [ 0.5, 1.0, 2.0, 4.0, ]
redend_options = [ True, False ]
orient_options = [ "RL", "LR", "TB", "BT" ]
notation_options = [ "snfg", "cfg" ]
display_options = [ "normal", "normalinfo", "compact" ]
opaque_options = [ True, False ]

gtc = GlyTouCan()
# gtc = GlyTouCanNoCache()

seen = set()

for j in range(iterations):
    imageWriter = GlycanImage()
    imageWriter.set('scale',random.choice(scale_options))
    imageWriter.set('reducing_end',random.choice(redend_options))
    imageWriter.set('orientation',random.choice(orient_options))
    imageWriter.set('notation',random.choice(notation_options))
    imageWriter.set('display',random.choice(display_options))
    imageWriter.set('opaque',random.choice(opaque_options))
    imageWriter.force(True)
    # imageWriter.verbose(True)

    for acc in random.sample(accs,k=batch):
        if acc in seen:
            continue
        outfile = acc + ".png"
        if os.path.exists(outfile):
            continue
        gly = gtc.getGlycan(acc,format='wurcs')
        if not gly:
            continue
        seq = gtc.getseq(acc,format='wurcs')
        if not seq:
            continue
        comp = gly.iupac_composition(floating_substituents=False,
                                     aggregate_basecomposition=False)
        bad = False
        for k,v in comp.items():
            if k in ('Glc','Gal','Man','NeuAc','NeuGc','Fuc','GlcNAc','GalNAc','Count'):
                continue
            if v <= 0:
                continue
            bad = True
            break
        if bad:
            continue
        imageWriter.writeImage(seq,outfile)
        wh = open(acc + ".log",'w')
        for k in ('scale','reducing_end','orientation','notation','display','opaque'):
            print(k+":",imageWriter.get(k),file=wh)
        print("composition:",comp,file=wh)
        wh.close()
