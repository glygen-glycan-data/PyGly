#!/bin/env python3
import sys, os.path, csv, copy, traceback
import findpygly
from pygly.Monosaccharide import Substituent, Mod, constantString, Anomer
from pygly.GlycanFormatter import WURCS20Format, GlycanParseError
from pygly.WURCS20MonoFormatter import  WURCS20MonoFormat
from pygly.manipulation import WURCSManipulation, WURCSArchetype
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCan, GlyTouCanNoPrefetch, GlyCosmos, GlyCosmosNoCache
from pygly.alignment import GlycanEqual, GlycanTopoEqual, GlycanImageEqual, GlycanSubsumption, MotifInclusive

snfgfile = sys.argv.pop(1)

wm = WURCSManipulation()
wp = WURCS20Format()
wmp = WURCS20MonoFormat()
subsumption = GlycanSubsumption()
motifincl = MotifInclusive()
glyeq = GlycanEqual()

# AUd1122h AUd21122h AUd21122h_5*NCC/3=O AUd21122h_5*NCCO/3=O AUdxxxxh
# AUdxxxxxh AUdxxxxxh_5*NCC/3=O AUdxxxxxh_5*NCCO/3=O h1112h h1112h_2*NCC/3=O
# h1122h h1122h_2*NCC/3=O h2112h h2112h_2*NCC/3=O h2121h h2122h h2122h_2*N
# h2122h_2*NCC/3=O h2212h h2212h_2*NCC/3=O h222h hU122h hxxxxh_2*N
# hxxxxh_2*NCC/3=O o1112h o1112h_2*NCC/3=O o1122h o1122h_2*NCC/3=O
# o2112A o2112h o2112h_2*NCC/3=O o2121h o2122A o2122h o2122h_2*N
# o2122h_2*NCC/3=O o2122m_2*NCC/3=O o212h o2212h o2212h_2*NCC/3=O o222h
# oxxxh oxxxxh oxxxxh_2*NCC/3=O u1112h u1112h_2*NCC/3=O u1112m u11221h
# u1122h u1122h_2*N u1122h_2*NCC/3=O u1221m u2111h_2*NCC/3=O u2112A u2112h
# u2112h_2*N u2112h_2*NCC/3=O u211h u2121A u2121h u2121h_2*NCC/3=O u2122A
# u2122h u2122h_2*N u2122h_2*NCC/3=O u212h u2211m u2212h u2212h_2*NCC/3=O
# u2222h u2222h_2*NCC/3=O u222h

extras = list("""
a2122h-1x_1-5_2*NSO/3=O/3=O
a1122h-1x_1-5_2*NSO/3=O/3=O
a2112h-1x_1-5_2*NSO/3=O/3=O
a2212h-1x_1-5_2*NSO/3=O/3=O
a2111h-1x_1-5_2*NSO/3=O/3=O
a1112h-1x_1-5_2*NSO/3=O/3=O
a2121h-1x_1-5_2*NSO/3=O/3=O
a2222h-1x_1-5_2*NSO/3=O/3=O
""".split())

if len(sys.argv) > 1:
    gtc = GlyTouCanNoPrefetch()
    iterable = sys.argv[1:]
else:
    gtc = GlyTouCan()
    iterable = gtc.allaccessions()

gco = GlyCosmos()
archived = set(map(lambda d: d['accession'],gco.archived()))

snfgmonos = []
for r in csv.DictReader(open(snfgfile),dialect='excel-tab'):
    if r['Color'] == 'White':
        continue
    mono = wp.toGlycan(r['WURCS'].strip()).root()
    snfgmonos.append(mono)   

for mstr in extras:
    mono = wmp.get(mstr)
    snfgmonos.append(mono)

# for m in snfgmonos:
#     print(m)

wurcsmonostr = dict()

for acc in iterable:
    if acc in archived:
        continue
    seq = gtc.getseq(acc,format='wurcs')
    if not seq:
        print(acc,False,"No WURCS sequence")
        continue

    comp = wm.deconstruct(seq)
    for md in comp['monos']:
        wms = wm.monoseq(md)
        if wms in wurcsmonostr:
            if wurcsmonostr[wms] == False:
                print(acc,False,"Bad monosaccharide: ["+wms+"]")
            elif wurcsmonostr[wms] == True:
                print(acc,True,"Good monosaccharide: ["+wms+"]")
            elif wurcsmonostr[wms] == None:
                print(acc,None,"Unknown monosaccharide: ["+wms+"]")
            continue

        try:
            m = wmp.get(wms)
        except GlycanParseError:
            wurcsmonostr[wms] = None
            good = False
            print(acc,None,"Unknown monosaccharide: "+wms+"]")
            continue

        tocheck = [ m ]
        m1 = m.clone()
        aldi = False
        for pos,mod in m1.mods():
            if mod == Mod.aldi and pos == (1,):
                aldi = True
                break
        subtoremove = []
        for sublink in m1.substituent_links():
            if sublink.child().name() in (Substituent.phosphate,Substituent.sulfate,Substituent.methyl):
                subtoremove.append(sublink)
        if aldi or len(subtoremove) > 0:
            if aldi:
                m1.remove_mod(Mod.aldi,(1,))
                m1.set_ring_start(None)
                m1.set_ring_end(None)
                m1.set_anomer(Anomer.missing)
            for sublink in subtoremove:
                m1.remove_substituent_link(sublink)
                sublink.child().del_parent_link(sublink)
                sublink.child().set_connected(False)
                sublink.child().clear_links()
            tocheck.append(m1)

        monogood = False
        for mi in tocheck:
            for snfgmono in snfgmonos:
                if subsumption.monoleq(snfgmono,mi) or motifincl.monoleq(snfgmono,mi):
                    monogood = True
                    break
            if monogood:
                break

        wurcsmonostr[wms] = monogood;
        if not monogood:
            good = False
            print(acc,False,"Bad monosaccharide: ["+wms+"]")
        else:
            print(acc,True,"Good monosaccharide: ["+wms+"]")
