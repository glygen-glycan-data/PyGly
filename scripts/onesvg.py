#!/bin/env python2

# This is just a terrible way to make one specific svg without waiting on the web service.
from __future__ import print_function

import sys, os, random, time, atexit
import findpygly
from pygly.GlycanImage import GlycanImage
from pygly.GlycanResource import GlyTouCan, GlyCosmos
from pygly.GNOme import GNOme

batch = 10
iterations = 100
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


imageWriter = GlycanImage()
imageWriter.set('scale',4.0)
imageWriter.set('reducing_end',False)
imageWriter.set('orientation',"TB")
imageWriter.set('notation',"snfg")
imageWriter.set('display',"compact")
# imageWriter.set('opaque',random.choice(opaque_options))
imageWriter.set('format','svg')
imageWriter.force(True)
# imageWriter.verbose(True)

#here is where I hardcoded the specific svg
acc = "G26091BM"

outfile = acc + ".svg"
if os.path.exists(outfile):
    sys.exit(1)
gly = gtc.getGlycan(acc,format='wurcs')
if not gly:
    sys.exit(1)
seq = gtc.getseq(acc,format='wurcs')
if not seq:
    sys.exit(1)
if gly.undetermined():
    sys.exit(1)
if not gly.has_root():
    sys.exit(1)
if gly.repeated():
    sys.exit(1)
topoacc = gnome.get_topology(acc)
if not topoacc:
    sys.exit(1)
comp = gly.iupac_composition(floating_substituents=False,
                             aggregate_basecomposition=False)
if comp['Count'] < 3:
    sys.exit(1)
bad = False
for k,v in comp.items():
    if k in ('Glc','Gal','Man','NeuAc','NeuGc','Fuc','GlcNAc','GalNAc','Count'):
        continue
    if v <= 0:
        continue
    bad = True
    break
if bad:
    sys.exit(1)
print(acc,file=sys.stderr)
imageWriter.writeImage(seq,outfile)
wh = open(acc + ".log",'w')
for k in ('scale','reducing_end','orientation','notation','display','opaque'):
    print(k+":",imageWriter.get(k),file=wh)
print("composition:",comp,file=wh)
print("topology:",topoacc,file=wh)
wh.close()
