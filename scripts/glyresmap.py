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
    canongly = gtc.getGlycan(acc,'wurcs')
    if not canongly:
        continue
    idmap = []
    if glyeq.eq(gly,canongly,idmap=idmap):
        bad = False
        glyids = [ m.external_descriptor_id() for m in gly.all_nodes(undet_subst=True) ]
        idmapglyids = [ t[0].external_descriptor_id() for t in idmap ]
        canonglyids = [ m.id() for m in canongly.all_nodes(undet_subst=True) ]
        canonidmapglyids = [ t[1].id() for t in idmap ]
        idmapids = [ ti for t in idmap for ti in glyeq.monoidmap(*t) ]
        if len(glyids) != len(idmap) or len(canonglyids) != len(idmap):
            print("Sequence %s:%s has an incorrect number of nodes in idmap: %d, %d, %d"%(acc,seqhash,len(glyids),len(idmap),len(canonglyids)))
            bad = True
        if (not gly.undetermined() and len(glyids) != len(set(glyids))) or len(canonglyids) != len(set(canonglyids)):
            print("Sequence %s:%s has non-unique monosaccharide ids: %s,%s"%(acc,seqhash,glyids,canonglyids))
            bad = True
        if (not gly.undetermined() and len(idmapglyids) != len(set(idmapglyids))) or len(canonidmapglyids) != len(set(canonidmapglyids)):
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
            # print(gly.external_descriptor_ids())
            # print(canongly.external_descriptor_ids())
            # for mj in canongly.all_nodes(subst=True,undet_subst=True):
            #     print(mj.external_descriptor_id())
            # for mi,mj in idmap:
            #     print(mi)
            #     print(mj)
            #     print(glyeq.monoidmap(mi,mj))
            print("Sequence %s:%s:%s: %s"%(acc,seqhash,fmt,", ".join(map(lambda t: "%s->%s"%t,sorted(idmapids)))))
        else:
            sys.exit(1)
    else:
        print("Sequence %s:%s:%s does not equal canonical WURCS sequence"%(acc,seqhash,fmt))
        print(gly.glycoct())
        print(canongly.glycoct())
        sys.exit(1)
