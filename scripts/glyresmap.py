#!/bin/env python2

import sys, glob, hashlib, os, os.path, traceback
import findpygly
from pygly.alignment import GlycanEqual, GlycanTopoEqual, GlycanImageEqual
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCan, GlyTouCanNoPrefetch
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError
from pygly.GlycanBuilderSVGParser import GlycanBuilderSVG

wp = WURCS20Format()
gp = GlycoCTFormat()
sp = GlycanBuilderSVG()

glyeq = GlycanImageEqual()

gtc = GlyTouCan()

for fn in sys.argv[1:]:
    acc = os.path.split(fn)[1].rsplit('.',1)[0]
    seqstr = open(fn).read().strip()
    seqhash = hashlib.md5(seqstr).hexdigest().lower()
    if seqstr.startswith('WURCS'):
        try:
            gly = wp.toGlycan(seqstr)
            fmt = 'WURCS'
        except GlycanParseError:
            continue
    elif seqstr.startswith('RES'):
        try:
            gly = gp.toGlycan(seqstr)
            fmt = 'GlycoCT'
        except GlycanParseError:
            continue
    elif fn.endswith('.svg'):
        try:
            gly = sp.toGlycan(seqstr)
            fmt = 'SVG'
        except GlycanParseError:
            print(acc)
            traceback.print_exc(file=sys.stdout)
            continue
    cannongly = gtc.getGlycan(acc,'wurcs')
    if not cannongly:
        continue
    idmap = []
    if glyeq.eq(gly,cannongly,idmap=idmap):
        bad = False
        glyids = [ m.external_descriptor_id() for m in gly.all_nodes() ]
        idmapglyids = [ t[0].external_descriptor_id() for t in idmap ]
        cannonglyids = [ m.id() for m in cannongly.all_nodes() ]
        cannonidmapglyids = [ t[1].id() for t in idmap ]
        idmapids = [ (t[0].external_descriptor_id(),t[1].id()) for t in idmap ]
        if len(glyids) != len(idmap) or len(cannonglyids) != len(idmap):
            print("Sequence %s:%s has an incorrect number of nodes in idmap: %d, %d, %d"%(acc,seqhash,len(glyids),len(idmap),len(cannonglyids)))
            bad = True
        if len(glyids) != len(set(glyids)) or len(cannonglyids) != len(set(cannonglyids)):
            print("Sequence %s:%s has non-unique monosaccharide ids: %s,%s"%(acc,seqhash,glyids,cannonglyids))
            bad = True
        if len(idmapglyids) != len(set(idmapglyids)) or len(cannonidmapglyids) != len(set(cannonidmapglyids)):
            print("Sequence %s:%s idmap has non-unique monosaccharide ids: %s"%(acc,seqhash,idmapids))
            bad = True
        if len(idmapids) != len(set(idmapids)):
            print("Sequence %s:%s idmap has non-unique monosaccharide id pairs: %s"%(acc,seqhash,idmapids))
            bad = True
        for mi,mj in idmap:
            if not glyeq.monoeq(mi,mj):
                print("Sequence %s:%s:%s has a mismatched monosaccharide"%(acc,seqhash,fmt))
                print(mi)
                print(mj)
                bad = True
        if not bad:
            print("Sequence %s:%s:%s: %s"%(acc,seqhash,fmt,", ".join(map(lambda t: "%s->%d"%t,sorted(idmapids)))))
        else:
            sys.exit(1)
    else:
        print("Sequence %s:%s:%s does not equal cannonical WURCS sequence"%(acc,seqhash,fmt))
        print(gly.glycoct())
        print(cannongly.glycoct())
        sys.exit(1)
