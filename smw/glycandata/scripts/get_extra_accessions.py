import sys

from getwiki import GlycanDataWiki, Glycan

w = GlycanDataWiki()

from pygly.GlyTouCan import GlyTouCan

gtc = GlyTouCan()

accs = set(open('../data/glytoucan_accessions.txt').read().split())
newaccs = set()

for acc in accs:
    topo = gtc.gettopo(acc)
    comp = gtc.getcomp(acc)
    if topo:
        newaccs.add(topo)
    if comp:
        for newacc in gtc.hascomp(comp):
            newaccs.add(newacc)    
        newaccs.add(comp)
newaccs = (newaccs-accs)

wh = open('../data/extra_accessions.txt','w')
for acc in sorted(newaccs):
    print >> wh, acc
wh.close()


