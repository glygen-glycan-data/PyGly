#!/bin/env python3
import sys, os.path, csv, copy, traceback
import findpygly
from pygly.Monosaccharide import Substituent, Mod, constantString, Anomer
from pygly.GlycanFormatter import WURCS20Format, GlycanParseError
from pygly.WURCS20MonoFormatter import  WURCS20MonoFormat
from pygly.manipulation import WURCSManipulation, WURCSArchetype
from pygly.GlycanResource import GlyTouCanNoCache, GlyTouCan, GlyTouCanNoPrefetch, GlyCosmos, GlyCosmosNoCache
from pygly.alignment import GlycanEqual, GlycanTopoEqual, GlycanImageEqual, GlycanSubsumption, MotifInclusive

snfgskel = sys.argv.pop(1)
snfgskel = set(open(snfgskel).read().split())
notsnfgskel = set()
while len(sys.argv) > 1 and os.path.exists(sys.argv[1]):
    notsnfgskel.update(set(open(sys.argv[1]).read().split()))
    sys.argv.pop(1)

if len(sys.argv) > 1:
    gtc = GlyTouCanNoPrefetch()
    iterable = sys.argv[1:]
else:
    gtc = GlyTouCan()
    iterable = gtc.allaccessions()

gco = GlyCosmos()
archived = set(map(lambda d: d['accession'],gco.archived()))

wm = WURCSManipulation()

for acc in iterable:
    if acc in archived:
        continue
    seq = gtc.getseq(acc,format='wurcs')
    if not seq:
        print(acc,False,"No WURCS sequence")
        continue
    allgood = True
    if '~' in seq:
        print(acc,False,"Repeating glycan")
        allgood = False
    if ',0+/' in seq:
        print(acc,False,"Zero-plus link count glycan")
        allgood = False
    elif seq.split('/',2)[1].endswith('+'):
        print(acc,False,"Undetermined link count glycan")
        allgood = False
    comp = wm.deconstruct(seq)
    floating = 0
    for l in comp['links']:
        if 'seq' in l and '-' in l['seq'] and '*' in l['seq']:
            print(acc,False,"Substituent link")
            allgood = False
            break
        if 'seq' in l and '-' not in l['seq'] and '*' in l['seq']:
            if l['seq'].endswith('*OSO/3=O/3=O') or l['seq'].endswith('*OPO/3O/3=O') or l['seq'].endswith('*OC'):
                pass
            else:
                print(acc,False,"Bad floating substituent")
                allgood = False
            floating += 1
    if comp['counts'][2] not in (0,comp['counts'][1]+floating-1):
        print(acc,False,"Bad link count glycan")
        allgood = False
    for m in comp['monos']:
        mseq = '['+wm.monoseq(m)+']'
        if mseq in snfgskel:
            print(acc,True,"Good monosaccharide: "+mseq)
        elif mseq in notsnfgskel:
            print(acc,False,"Bad monosaccharide: "+mseq)
            allgood = False
        else:
            print(acc,None,"Unknown monosaccharide: "+mseq)
            allgood = False
    if allgood:
        print(acc,True,"SNFG glycan")
