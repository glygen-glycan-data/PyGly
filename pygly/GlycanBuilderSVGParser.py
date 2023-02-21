from __future__ import print_function

from . Monosaccharide import *
from . Glycan import Glycan
from . MonoFactory import MonoFactory

import xml.etree.ElementTree as ET

import re, sys, traceback
import copy
import string
from collections import defaultdict

import StringIO

from . GlycanFormatter import GlycanFormatter
from . GlycanFormatterExceptions import *

class GlycanBuilderSVGMonosaccharideLookupError(GlycanBuilderSVGParseError):
    def __init__(self,res):
        self.message = "GlycanBuilderSVG parser: Bad residue name %s."%(res,)

class GlycanBuilderSVGSubstituentLookupError(GlycanBuilderSVGParseError):
    def __init__(self,res):
        self.message = "GlycanBuilderSVG parser: Bad substituent residue name %s."%(res,)

class GlycanBuilderSVGRingLookupError(GlycanBuilderSVGParseError):
    def __init__(self,ringstr):
        self.message = "GlycanBuilderSVG parser: Bad ring value %s."%(ringstr,)

class GlycanBuilderSVGMissingRootError(GlycanBuilderSVGParseError):
    def __init__(self):
        self.message = "GlycanBuilderSVG parser: Can't determine the root monosaccharide."

class GlycanBuilderSVG(GlycanFormatter):
    substmap = {'S': {'type': Substituent.sulfate, 
                      'parent_type': Linkage.oxygenPreserved, 
                      'child_type': Linkage.nitrogenAdded},
                'P': {'type': Substituent.phosphate, 
                      'parent_type': Linkage.oxygenPreserved, 
                      'child_type': Linkage.nitrogenAdded},
                'N': {'type': Substituent.amino, 
                      'parent_type': Linkage.oxygenLost, 
                      'child_type': Linkage.nitrogenAdded},
                'Ac': {'type': Substituent.acetyl,
                       'parent_type': Linkage.oxygenPreserved, 
                       'child_type': Linkage.nitrogenAdded},
                'Gc': {'type': Substituent.glycolyl,
                       'parent_type': Linkage.oxygenPreserved, 
                       'child_type': Linkage.nitrogenAdded},
                'Me': {'type': Substituent.methyl,
                       'parent_type': Linkage.oxygenPreserved, 
                       'child_type': Linkage.nitrogenAdded},
                'NAc': {'type': Substituent.nAcetyl,
                        'parent_type': Linkage.oxygenLost, 
                        'child_type': Linkage.nitrogenAdded},
               }
    ringmap = {'p': {'ring_start': 1, 'ring_end': 5},
               'f': {'ring_start': 1, 'ring_end': 4},
               'o': {'ring_start': 0, 'ring_end': 0},
              }
    def __init__(self):
        self.monofactory = MonoFactory()

    def toGlycan(self,s):
        s = s.replace('data:','data.')
        doc = ET.parse(StringIO.StringIO(s))
        monos = dict()
        rootmonoindex = 1
        root = doc.getroot()
	ns = root.tag.split('}')[0]+'}'
        for g in root.getiterator(ns+'g'):
            if g.attrib.get('data.type') in ('Residue','Monosaccharide'):
                 # handle monosaccharide
                 res = g.attrib['data.residueName']
                 try:
                     m = self.monofactory.new(res)
                 except KeyError:
                     raise GlycanBuilderSVGMonosaccharideLookupError(res)
                 anomer = g.attrib['data.residueAnomericState']
                 if anomer == 'a':
                     m.set_anomer(Anomer.alpha)
                 elif anomer == 'b':
                     m.set_anomer(Anomer.beta)
                 elif anomer == '?':
                     m.set_anomer(Anomer.missing)
                 ringstr = g.attrib['data.residueRingSize']
                 if ringstr != "?":
                     if ringstr in self.ringmap:
                         m.set_ring_start(self.ringmap[ringstr]['ring_start'])
                         m.set_ring_end(self.ringmap[ringstr]['ring_end'])
                         if ringstr == "o":
                             m.set_anomer(Anomer.uncyclized)
                     else:
                         raise GlycanBuilderSVGRingLookupError(ringstr)
                 resid = int(g.attrib['data.residueIndex'])
                 if g.attrib.get('data.residueIsReducingEnd') == "true":
                     rootmonoindex = resid
                 if g.attrib.get('data.residueIsAlditol') == "true":
                     m.add_mod(1,Mod.aldi)
                 m.set_id(resid)
                 m.set_external_descriptor_id(g.attrib['ID'])
                 monos[resid] = m
            elif g.attrib.get('data.type') == 'Substituent':
                # handle substituent
                subres = g.attrib['data.residueName']
                if subres not in self.substmap:
                    raise GlycanBuilderSVGSubstituentLookupError(subres)
                sub = Substituent(self.substmap[subres]['type'])
                if g.attrib['data.parentPositions'] == "?":
                    parent_pos = None
                else:
                    parent_pos = set(map(int,g.attrib['data.parentPositions'].split('/')))
                if g.attrib['data.childPositions'] == "?":
                    child_pos = None
                else:
                    child_pos = set(map(int,g.attrib['data.childPositions'].split('/')))
                parent = monos[int(g.attrib['data.parentResidueIndex'])]
                parent.add_substituent(sub, 
                                       parent_pos=parent_pos, parent_type=self.substmap[subres]['parent_type'], 
                                       child_pos=child_pos, child_type=self.substmap[subres]['child_type'])
                parent_eid = parent.external_descriptor_id()
                parent_eid += ';'+g.attrib['ID']
                parent.set_external_descriptor_id(parent_eid)

        for g in root.getiterator(ns+'g'):
            if g.attrib.get('data.type') == 'Linkage':
                if not g.attrib.get('data.parentResidueIndex') or not g.attrib.get('data.childResidueIndex'):
                    continue
                parent = monos[int(g.attrib['data.parentResidueIndex'])]
                if g.attrib['data.parentPositions'] == "?":
                    parentpos = None
                else:
                    parentpos = set(map(int,g.attrib['data.parentPositions'].split('/')))
                child = monos[int(g.attrib['data.childResidueIndex'])]
                if g.attrib['data.childPositions'] == "?":
                    childpos = None
                else:
                    childpos = set(map(int,g.attrib['data.childPositions'].split('/')))
                parent.add_child(child, parent_pos=parentpos, child_pos=childpos,
                                        parent_type=Linkage.oxygenPreserved,
                                        child_type=Linkage.oxygenLost)

        if rootmonoindex not in monos:
            raise GlycanBuilderSVGMissingRootError()
        rootmono = monos[rootmonoindex]
        gly = Glycan(rootmono)
        return gly
