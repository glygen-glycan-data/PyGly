#!/bin/env python27

from getwiki import GlycoMotifWiki, AllMotif

w = GlycoMotifWiki()
AllMotifpageid = AllMotif.id

from glycancomparison import GlycanTopologySameAs
gtsa = GlycanTopologySameAs()

from GlycanFormatter import GlycoCTFormat, WURCS20Format
wurcs_parser = WURCS20Format()
glycoct_parser = GlycoCTFormat()


glycanobj = {}
for m in w.itermotif():

    gtcid = m.get("glytoucan")
    if gtcid in glycanobj:
        continue

    gseq = m.get("glycoct")
    objloaded = False
    if gseq:
        try:
            g = glycoct_parser.toGlycan(gseq)
            objloaded = True
        except:
            pass

    if not objloaded:
        wseq = m.get("wurcs")
        if wseq:
            try:
                g = wurcs_parser.toGlycan(wseq)
            except:
                continue
        else:
            continue

    glycanobj[gtcid] = g
print "Total readable glycans: %s " % len(glycanobj)


matchpool = []
for acc1, g1 in glycanobj.items():
    for acc2, g2 in glycanobj.items():

        if acc1 == acc2:
            continue

        if gtsa.get(g1, g2):
            res = set([acc1, acc2])
            matchpool.append(res)

# Cluster the comparison results
anotherround = True
while (anotherround):
    anotherround = False
    for i, s1 in enumerate(matchpool):
        for j, s2 in enumerate(matchpool):

            if j <= i:
                continue

            if s1.intersection(s2):
                new = s1.union(s2)
                matchpool.append(new)
                anotherround = True
                matchpool.pop(j)
                matchpool.pop(i)
                break

        if anotherround:
            break

glycaninthematchpool = set()
for s in matchpool:
    glycaninthematchpool = glycaninthematchpool.union(s)


for m in w.itermotif():
    
    if m.get("collection") != AllMotifpageid:
        continue

    gtcid = m.get("glytoucan")
    if not gtcid in glycaninthematchpool:
        equivalent = AllMotifpageid + "." + gtcid
        equivalents = [ equivalent ]
        m.set("topology", equivalents)
        w.put(m)
        print "Singleton %s" % m.get("id")
        continue

    equivalents = []
    for matchres in matchpool:
        if gtcid in matchres:
            equivalents = list(matchres)
            break

    equivalents = map(lambda x: str(AllMotifpageid + "." + x), equivalents)

    print m.get("id")
    m.set("topology", equivalents)
    w.put(m)
