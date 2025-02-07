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

snfgfiles = []
for p in sys.argv[1:]:
    f = os.path.split(p)[1]
    if f.startswith('snfg') and f.endswith('.tsv'):
        snfgfiles.append(p)
    else:
        break
for i in range(len(snfgfiles)):
    sys.argv.pop(1)

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

for f in snfgfiles:
  if not os.path.split(f)[1].startswith('snfgmono'):
    continue
  for r in csv.DictReader(open(f),dialect='excel-tab'):
    if r['Color'] == 'White':
        continue
    md = wm.deconstruct(r['WURCS'])['monos'][0]
    if md['wurcs'] in [ t[0]['wurcs'] for t in snfgmonos ]:
        assert False, "Repeated monosaccharide WURCS: "+md['wurcs']+"."
    mono = wmp.get(md['wurcs'])
    snfgmonos.append((md,mono))

snfgmonos.sort(key=lambda t: -len(t[0].get('substituents',[])))

snfgsubsts = set()

for f in snfgfiles:
  if not os.path.split(f)[1].startswith('snfgsubst'):
    continue
  for r in csv.DictReader(open(f),dialect='excel-tab'):
    subst = r['WURCS']
    s = wmp.getsubst(subst[1:])
    if subst[1:] in snfgsubsts:
        assert False, "Repeated substituent WURCS: "+subst+"."
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
    if aldi or mono.noring():
        mono1 = mono.clone()
        mono1.remove_mod(Mod.aldi,(1,))
        mono1.set_ring_start(None)
        mono1.set_ring_end(None)
        if subsumption.monoleq(snfgmono,mono1):
            return True
    return False

def checkmono(md):
    subst = md['substituents']
    nsubst = len(subst)
    ndsubst= len(set(sub['subst'] for sub in subst))
    for sub in subst:
        if not checksubst(sub['subst'],extra=True)[0]:
            return False,nsubst,ndsubst,0,0,"No matching substituent:",sub['subst']
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
                return False,nsubst,ndsubst,0,0,"No matching monosaccharide."
            result = testmono(snfgmono,mono)
            result1 = True
            for sub in rest_subst1:
                if not checksubst(sub['subst'])[0]:
                    result1 = False
                    break
            if result and result1:
                return True,nsubst,ndsubst,len(rest_subst1),len(set(sub['subst'] for sub in rest_subst1)),"Matched monosaccharide:",snfgmd['wurcs']
    for sub in subst:
        if not checksubst(sub['subst'])[0]:
            return False,nsubst,ndsubst,0,0,"No matching substituent:",sub['subst']
    return False,nsubst,ndsubst,0,0,"No matching monosaccharide."

def checksubst(wurcsseq,extra=False):
    for snfgwurcsseq in snfgsubsts:
        if wurcsseq == snfgwurcsseq:
            return True,"Matched substituent",wurcsseq
    if extra:
        for extrawurcsseq in extrasubsts:
            if wurcsseq == extrawurcsseq:
                return True,"Matched extra substituent",wurcsseq
    return False,"No matching substituent",wurcsseq

resultcache = dict()

print("\t".join(["Accession","#Mono","#FlSubst","Type","Matched?","WURCS","#Subst","#Distinct","#ExtraSubst","#DistinctExtraSubst","Explanation","MatchedWURCS"]))

for acc in iterable:
    seq = None
    if re.search(r"^G\d{5}[A-Z]{2}$",acc):
        if acc in archived:
            continue
        seq = gtc.getseq(acc,format='wurcs')
    else:
        acc,seq = acc.split(':',1)
    if not seq:
        # print(acc,False,"No WURCS sequence")
        continue

    # print(acc,seq)

    comp = wm.deconstruct(seq)
    # print(comp)
    nmono = len(comp['monos'])
    nsubst = len([ ld.get('subst') for ld in comp['links'] if ld.get('subst') ])
    for wmd in comp['monos']:
        wms = wmd['wurcs']
        if wms in resultcache:
            print("\t".join(map(str,[acc,nmono,nsubst,"Monosaccharide",resultcache[wms][0],wms,*resultcache[wms][1:]])))
            continue

        result = checkmono(wmd)
        resultcache[wms] = result
        print("\t".join(map(str,[acc,nmono,nsubst,"Monosaccharide",resultcache[wms][0],wms,*resultcache[wms][1:]])))

    for ld in comp['links']:
        wss = ld.get('subst')
        if not wss:
            continue

        if wss in resultcache:
            print("\t".join(map(str,[acc,nmono,nsubst,"Substituent",resultcache[wss][0],wss,"","","","",*resultcache[wss][1:]])))
            continue

        result = checksubst(wss)
        resultcache[wss] = result
        print("\t".join(map(str,[acc,nmono,nsubst,"Substituent",resultcache[wss][0],wss,"","","","",*resultcache[wss][1:]])))
