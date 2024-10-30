#!/bin/env python3.12
import sys, os.path, csv, copy, traceback, re
import findpygly
from pygly.Monosaccharide import Substituent, Mod, constantString, Anomer
from pygly.GlycanFormatter import GlycanParseError
from pygly.WURCS20MonoFormatter import  WURCS20MonoFormat
from pygly.manipulation import WURCSManipulation
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCan, GlyTouCanNoPrefetch, GlyCosmos, GlyCosmosNoCache
from pygly.alignment import GlycanSubsumption, MotifInclusive
from pygly.combinatorics import itersubmatchings
from operator import itemgetter

wm = WURCSManipulation()
wmp = WURCS20MonoFormat()
subsumption = GlycanSubsumption()
motifincl = MotifInclusive()

snfgfile = sys.argv[1]
sys.argv.pop(1)

snfgsubstwurcs = list(filter(None,"""
*OCC/3=O
*OCC^RC/4N/3=O
*OCC^RNCC/6=O/4C/3=O
*N=^XCC/3N
*N=^XCC/3N/2C
*N=^XCNC/4C/3C
*OC=O
*OCCO/3=O
*OCC^XNCC/6=O/4CCCN/11=O/3=O
*OCCCC^XCN/7=O/6NC/3=O
*OCCN/3=O
*OCC^XCO/4O/3=O
*OCC^XCOC/4OC/3=O
*OCCCCO/3=O
*OCCC^XCO/5O/3=O
*OCCC^RC/5O/3=O
*OCCC^SC/5O/3=O
*OCC^XC/4O/3=O
*OC
*N
*NCC/3=O
*NCCO/3=O
*NSO/3=O/3=O
*OPO/3O/3=O
*OCCC/4=O/3=O
*OC^XO*/3CO/6=O/3C
*OSO/3=O/3=O
*NCCCSO/6=O/6=O
""".splitlines()))

snfgmonowurcs = list(filter(None,"""
#axxxxh-1x_1-5
a2122h-1x_1-5
a1122h-1x_1-5
a2112h-1x_1-5
a2212h-1x_1-5
a2111h-1x_1-5
a2222h-1x_1-5
a1112h-1x_1-5
a2121h-1x_1-5
#axxxxh-1x_1-5_2*NCC/3=O
a2122h-1x_1-5_2*NCC/3=O
a1122h-1x_1-5_2*NCC/3=O
a2112h-1x_1-5_2*NCC/3=O
a2112h-1x_1-5_2*NCC/3=O
a2212h-1x_1-5_2*NCC/3=O
a2111h-1x_1-5_2*NCC/3=O
a2222h-1x_1-5_2*NCC/3=O
a1112h-1x_1-5_2*NCC/3=O
a2121h-1x_1-5_2*NCC/3=O
#axxxxh-1x_1-5_2*N
a2122h-1x_1-5_2*N
a1122h-1x_1-5_2*N
a2112h-1x_1-5_2*N
a2212h-1x_1-5_2*N
a2212h-1x_1-5_2*N
a2111h-1x_1-5_2*N
a2222h-1x_1-5_2*N
a1112h-1x_1-5_2*N
a2121h-1x_1-5_2*N
#axxxxA-1x_1-5
a2122A-1x_1-5
a1122A-1x_1-5
a2112A-1x_1-5
a2212A-1x_1-5
a2111A-1x_1-5
a2222A-1x_1-5
a1112A-1x_1-5
a2121A-1x_1-5
#axxxxm-1x_1-5
a2122m-1x_1-5
a2211m-1x_1-5
a2212m-1x_1-5
a2111m-1x_1-5
a1112m-1x_1-5
a1221m-1x_1-5
#axxxxm-1x_1-5_2*NCC/3=O
a2122m-1x_1-5_2*NCC/3=O
a2211m-1x_1-5_2*NCC/3=O
a2111m-1x_1-5_2*NCC/3=O
a1112m-1x_1-5_2*NCC/3=O
a1221m-1x_1-5_2*NCC/3=O
#adxxxm-1x_1-5
ad122m-1x_1-5
a1d22m-1x_1-5
a2d12m-1x_1-5
a2d22m-1x_1-5
ad222m-1x_1-5
a1d21m-1x_1-5
#axxxh-1x_1-5
a211h-1x_1-5
a112h-1x_1-5
a212h-1x_1-5
a222h-1x_1-5
#Aadxxxxxh-2x_2-6
Aad21122h-2x_2-6
Aad21122h-2x_2-6_5*NCC/3=O
Aad21122h-2x_2-6_5*NCCO/3=O
Aad21122h-2x_2-6_5*N
#Aadxxxxxm-2x_2-6_5*N_7*N
Aad22111m-2x_2-6_5*N_7*N
Aad21122m-2x_2-6_5*N_7*N
Aad21111m-2x_2-6_5*N_7*N
Aad11122m-2x_2-6_5*N_7*N
a2122m-1x_1-5_2*N_4*N
a11221h-1x_1-5
Aad1122h-2x_2-6
Aad112A-2x_2-6
a11222h-1x_1-5
a2122h-1x_1-5_2*NCC/3=O_3*OC^RCO/4=O/3C
a2122h-1x_1-5_2*NCCO/3=O_3*OC^RCO/4=O/3C
a2122h-1x_1-5_2*N_3*OC^RCO/4=O/3C
a26h-1x_1-4_3*CO
ha122h-2x_2-6
ha112h-2x_2-6
ha121h-2x_2-6
ha222h-2x_2-6
""".splitlines()))

def lines(filename):
    for i,l in enumerate(open(filename)):
        sl = l.split(None,1)
        if re.search(r'^G\d{5}[A-Z]{2}$',sl[0]):
            yield "%s:%s"%(sl[0],sl[1].strip())
        else:
            yield "L%06d:%s"%(i+1,l.strip())

if len(sys.argv) > 1:
    if os.path.exists(sys.argv[1]):
        iterable = lines(sys.argv[1])
    else:
        gtc = GlyTouCanNoPrefetch()
        iterable = sys.argv[1:]
else:
    gtc = GlyTouCanNoCache()
    iterable = sorted(gtc.allaccessions())

gco = GlyCosmosNoCache(verbose=False)
archived = set(map(lambda d: d['accession'],gco.archived()))

snfgmonos = []

for r in csv.DictReader(open(snfgfile),dialect='excel-tab'):
    if r['Color'] == 'White':
        continue
    md = wm.deconstruct(r['WURCS'])['monos'][0]
    mono = wmp.get(md['wurcs'])
    snfgmonos.append((md,mono))

for mono in snfgmonowurcs:
    mono = mono.strip()
    if mono.startswith('#'):
        continue
    if mono in [ t[0]['wurcs'] for t in snfgmonos]:
        continue
    assert False, "Unexpected monosaccharide WURCS: "+mono+"."
    m = wmp.get(mono)
    md = wm.monodeconstruct(mono)
    snfgmonos.append((md,m))

snfgmonos.sort(key=lambda t: -len(t[0].get('substituents',[])))

snfgsubsts = set()

for subst in snfgsubstwurcs:
    subst = subst.strip()
    if subst.startswith('#'):
        continue
    s = wmp.getsubst(subst[1:])
    snfgsubsts.add(subst[1:])

extrasubsts = set()

for monos in snfgmonos:
    for sd in monos[0]['substituents']:
        if sd['subst'] not in snfgsubsts:
            extrasubsts.add(sd['subst'])

def testmono(snfgmono,mono):
    mono.set_anomer(Anomer.missing)
    if subsumption.monoleq(snfgmono,mono):
         return True
    aldi = False
    for pos,mod in mono.mods():
        if mod == Mod.aldi and pos == (1,):
            aldi = True
            break
    if aldi:
        mono1 = mono.clone()
        mono1.remove_mod(Mod.aldi,(1,))
        mono1.set_ring_start(None)
        mono1.set_ring_end(None)
        if subsumption.monoleq(snfgmono,mono1):
            return True
    return False

def checkmono(md):
    subst = md['substituents']
    for sub in subst:
        if not checksubst(sub['subst'],extra=True)[0]:
            return False,"No matching SNFG substituent:",sub['subst']
    for snfgmd,snfgmono in snfgmonos:
        snfgsubst = snfgmd['substituents']
        for snfgsubst1,subst1,rest_subst1 in itersubmatchings(snfgsubst,subst,lambda a,b: a['subst'] == b['subst']):
            md1 = copy.copy(md)
            del md1['wurcs']
            md1['substituents'] = subst1
            wms = wm.monoseq(md1)
            try: 
                mono = wmp.get(wms)
            except GlycanParseError:
                return False,"No matching SNFG monosaccharde:",None
            result = testmono(snfgmono,mono)
            result1 = True
            for sub in rest_subst1:
                if not checksubst(sub['subst'])[0]:
                    result1 = False
                    break
            if result and result1:
                return True,"Matched SNFG monosaccharide:",snfgmd['wurcs']
    for sub in subst:
        if not checksubst(sub['subst'])[0]:
            return False,"No matching SNFG substituent:",sub['subst']
    return False,"No matching SNFG monosaccharide:",None

def checksubst(wurcsseq,extra=False):
    for snfgwurcsseq in snfgsubsts:
        if wurcsseq == snfgwurcsseq:
            return True,"Matched SNFG substituent",wurcsseq
    if extra:
        for extrawurcsseq in extrasubsts:
            if wurcsseq == extrawurcsseq:
                return True,"Matched extra substituent",wurcsseq
    return False,"No matching SNFG substituent",None

resultcache = dict()

for acc in iterable:
    seq = None
    if re.search(r"^G\d{5}[A-Z]{2}$",acc):
        if acc in archived:
            continue
        seq = gtc.getseq(acc,format='wurcs')
    else:
        acc,seq = acc.split(':',1)
    if not seq:
        print(acc,False,"No WURCS sequence")
        continue

    # print(acc,seq)

    comp = wm.deconstruct(seq)
    # print(comp)
    for wmd in comp['monos']:
        wms = wmd['wurcs']
        if wms in resultcache:
            print(acc,"Monosaccharide",resultcache[wms][0],wms,*resultcache[wms][1:])
            continue

        result = checkmono(wmd)
        resultcache[wms] = result
        print(acc,"Monosaccharide",resultcache[wms][0],wms,*resultcache[wms][1:])

    for ld in comp['links']:
        wss = ld.get('subst')
        if not wss:
            continue

        if wss in resultcache:
            print(acc,"Substituent",resultcache[wss][0],wss,*resultcache[wss][1:])
            continue

        result = checksubst(wss)
        resultcache[wss] = result
        print(acc,"Substituent",resultcache[wss][0],wss,*resultcache[wss][1:])
