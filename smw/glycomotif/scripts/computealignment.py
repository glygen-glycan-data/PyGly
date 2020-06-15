import os
import sys
import csv
from pygly.GlycanResource import GlyTouCan
from pygly.GlycanFormatter import WURCS20Format, GlycoCTFormat
from pygly.alignment import GlyTouCanMotif


wp = WURCS20Format()
gp = GlycoCTFormat()

gtcm = GlyTouCanMotif()
gtc = GlyTouCan()


from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()
motifs = {}
for m in w.itermotif():

    acc = m.get("glytoucan")
    if acc in motifs:
        continue

    try:
        motifs[acc] = wp.toGlycan(str(m.get("wurcs")))
    except:

        try:
            motifs[acc] = gp.toGlycan(m.get("glycoct"))
        except:
            continue


accs = set()
res = []
for acc, f, s in gtc.allseq(format='wurcs'):
    try:
        g = wp.toGlycan(s)
    except:
        continue

    for macc, motif in motifs.items():

        a = gtcm.leq(motif, g, rootOnly=True, anywhereExceptRoot=False)
        b = gtcm.leq(motif, g, rootOnly=False, anywhereExceptRoot=True)

        if a or b:
            l = "%s\t%s\t%s" % (macc, acc, a)
            res.append(l)


fp = sys.argv[1]
# fp = "../data/motif_alignment.tsv"
output_file = open(fp, 'w')
for l in sorted(res):
    output_file.write(l + "\n")
output_file.close()







