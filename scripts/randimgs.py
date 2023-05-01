#!/bin/env python2
from __future__ import print_function

import sys, os, random, time, atexit
import findpygly
from pygly.GlycanImage import GlycanImage
from pygly.GlycanResource import GlyTouCan, GlyCosmos
from pygly.GNOme import GNOme

imagenum = int(sys.argv[1]) if len(sys.argv) > 1 else 100
mode = sys.argv[2] if len(sys.argv) > 2 else 'png'
badaccfile = sys.argv[3] if len(sys.argv) > 3 else None

batch = 10
iterations = imagenum//batch
scale_options = [ 0.5, 1.0, 2.0, 4.0, ]
redend_options = [ True, False ]
orient_options = [ "RL", "LR", "TB", "BT" ]
notation_options = [ "snfg", "cfg" ]
display_options = [ "normal", "normalinfo", "compact" ]
opaque_options = [ True, False ]

print("GlyCosmos archived...",file=sys.stderr)
start = time.time()
gco = GlyCosmos(usecache=False)
archived = set(map(lambda d: d['accession'],gco.archived()))
print("GlyCosmos archived complete. (%s secs.)"%(time.time()-start,),file=sys.stderr)

print("GlyTouCan accessions...",file=sys.stderr)
start = time.time()
gtc = GlyTouCan(verbose=False,usecache=False,prefetch=True)
accs = list(filter(lambda acc: acc not in archived,gtc.allaccessions()))
dummy = gtc.getseq('G00912UN','wurcs')
print("GlyTouCan accessions complete. (%s secs.)"%(time.time()-start,),file=sys.stderr)

print("GNOme setup...",file=sys.stderr)
start = time.time()
gnome = GNOme()
print("GNOme setup complete. (%s secs.)"%(time.time()-start,),file=sys.stderr)

import subprocess, atexit
xvfbproc = subprocess.Popen(["Xvfb",":1"])
os.environ["DISPLAY"] = "localhost:1.0"
def killxvfb(proc):
    try:
        proc.terminate()
    except OSError:
        pass
atexit.register(killxvfb,xvfbproc)

if badaccfile is not None:
#uses the accessions your model was trained on, to avoid testing on them
    trained_accessions = set()
    with open(badaccfile) as f:
        for l in f:
            trained_accessions.add(l.rstrip())

seen = set()
for j in range(iterations):
    imageWriter = GlycanImage()
    imageWriter.set('scale',random.choice(scale_options))
    imageWriter.set('reducing_end',random.choice(redend_options))
    imageWriter.set('orientation',random.choice(orient_options))
    imageWriter.set('notation',random.choice(notation_options))
    imageWriter.set('display',random.choice(display_options))
    #imageWriter.set('opaque',random.choice(opaque_options))
    imageWriter.set('format',mode)
    imageWriter.force(True)
    # imageWriter.verbose(True)

    count = 0
    while count < batch:
        acc = random.choice(accs)
        if acc in seen:
            continue
        seen.add(acc)
        outfile = acc + "." + mode
        if os.path.exists(outfile):
            continue
        gly = gtc.getGlycan(acc,format='wurcs')
        if not gly:
            continue
        seq = gtc.getseq(acc,format='wurcs')
        if not seq:
            continue
        if gly.undetermined():
            continue
        if not gly.has_root():
            continue
        if gly.repeated():
            continue
        topoacc = gnome.get_topology(acc)
        if not topoacc:
            continue
        comp = gly.iupac_composition(floating_substituents=False,
                                     aggregate_basecomposition=False)
        if comp['Count'] < 3:
            continue
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
        print(acc,file=sys.stderr)
        imageWriter.writeImage(seq,outfile)
        wh = open(acc + ".log",'w')
        for k in ('scale','reducing_end','orientation','notation','display','opaque'):
            print(k+":",imageWriter.get(k),file=wh)
        print("composition:",comp,file=wh)
        print("topology:",topoacc,file=wh)
        wh.close()
        count += 1
