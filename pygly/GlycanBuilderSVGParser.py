from __future__ import print_function

from . Monosaccharide import *
from . Glycan import Glycan
from . MonoFactory import MonoFactory

import xml.etree.ElementTree as ET

import re, sys, traceback
import copy
import string
from collections import defaultdict

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
try:
    from io import BytesIO
except ImportError:
    pass

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

class GlycanBuilderSVGUnsupportedStructureError(GlycanBuilderSVGParseError):
    def __init__(self):
        self.message = "GlycanBuilderSVG parser: Unsupported structure."

class GlycanBuilderSVGUnexpectedLinkConnectionError(GlycanBuilderSVGParseError):
    def __init__(self):
        self.message = "GlycanBuilderSVG parser: Unexpected link connection."

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
                'NS': {'type': Substituent.nsulfate, 
                       'parent_type': Linkage.oxygenPreserved, 
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
                'F': {'type': Substituent.fluoro,
                      'parent_type': Linkage.oxygenLost, 
                      'child_type': Linkage.nitrogenAdded},
                'I': {'type': Substituent.iodo,
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
        s = s.replace(b'data:',b'data.')
        doc = ET.parse(BytesIO(s))
        undetroot = dict()
        undetrootlinkpos = dict()
        monos = dict()
        rootmonoindex = 1
        root = doc.getroot()
        ns = root.tag.split('}')[0]+'}'
        for g in root.iter(ns+'g'):
            if g.attrib.get('data.type') in ('Residue','Monosaccharide'):
                 # handle monosaccharide
                 res = g.attrib['data.residueName']
                 try:
                     m = self.monofactory.new(res)
                 except KeyError as e:
                     m = None
                 if not m:
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
                 if g.attrib.get('data.residueMultiplicity') != None:
                     undetroot[resid] = int(g.attrib.get('data.residueMultiplicity'))
                 if g.attrib.get('data.residueUndeterminedMultiplicity') != None:
                     undetroot[resid] = int(g.attrib.get('data.residueUndeterminedMultiplicity'))
                     linkpos = [g.attrib.get('data.residueUndeterminedChildPos'),
                                g.attrib.get('data.residueUndeterminedParentPos')]
                     try:
                         linkpos[0] = int(linkpos[0])
                     except (ValueError,TypeError):
                         if linkpos[0] == '?':
                             linkpos[0] = -1
                     try:
                         linkpos[1] = int(linkpos[1])
                     except (ValueError,TypeError):
                         if linkpos[1] == '?':
                             linkpos[1] = -1
                     undetrootlinkpos[resid] = tuple(linkpos)
                  
            elif g.attrib.get('data.type') == 'Substituent':
                # handle substituent
                subres = g.attrib['data.residueName']
                if subres not in self.substmap:
                    raise GlycanBuilderSVGSubstituentLookupError(subres)
                sub = Substituent(self.substmap[subres]['type'])
                if 'data.residueIndex' in g.attrib:
                    # floating substituent
                    resid = int(g.attrib['data.residueIndex'])
                    sub.set_id(resid)
                    sub.set_external_descriptor_id(g.attrib['ID'])
                    monos[resid] = sub
                    
                    assert g.attrib.get('data.residueUndeterminedMultiplicity') != None
                    undetroot[resid] = int(g.attrib.get('data.residueUndeterminedMultiplicity'))
                    linkpos = [g.attrib.get('data.residueUndeterminedChildPos'),
                               g.attrib.get('data.residueUndeterminedParentPos')]
                    try:
                        linkpos[0] = int(linkpos[0])
                    except (ValueError,TypeError):
                        if linkpos[0] == '?':
                            linkpos[0] = -1
                    try:
                        linkpos[1] = int(linkpos[1])
                    except (ValueError,TypeError):
                        if linkpos[1] == '?':
                            linkpos[1] = -1
                    undetrootlinkpos[resid] = tuple(linkpos)

                else:
                    if g.attrib['data.parentPositions'] == "?":
                        parent_pos = None
                    else:
                        try:
                            parent_pos = set(map(int,g.attrib['data.parentPositions'].split('/')))
                        except (KeyError,ValueError):
                            raise GlycanBuilderSVGUnexpectedLinkConnectionError()
                    if g.attrib['data.childPositions'] == "?":
                        child_pos = None
                    else:
                        try:
                            child_pos = set(map(int,g.attrib['data.childPositions'].split('/')))
                        except (KeyError,ValueError):
                            raise GlycanBuilderSVGUnexpectedLinkConnectionError()
                    try:
                        parent = monos[int(g.attrib['data.parentResidueIndex'])]
                    except (KeyError,ValueError):
                        raise GlycanBuilderSVGUnexpectedLinkConnectionError()

                    if sub.name() == Substituent.sulfate and parent_pos == set([2,]):
                        for sl in parent.substituent_links():
                            if sl.child().name() == Substituent.amino and sl.parent_pos() == set([2,]):
                                subres = "NS"
                                parent.remove_substituent_link(sl)
                                sub = Substituent(self.substmap[subres]['type'])
                                break

                    parent.add_substituent(sub, 
                                           parent_pos=parent_pos, parent_type=self.substmap[subres]['parent_type'], 
                                           child_pos=child_pos, child_type=self.substmap[subres]['child_type'])
                    sub.set_external_descriptor_id(g.attrib['ID']);
                    # parent_eid = parent.external_descriptor_id()
                    # parent_eid += ';'+g.attrib['ID']
                    # parent.set_external_descriptor_id(parent_eid)

        for g in root.iter(ns+'g'):
            if g.attrib.get('data.type') == 'Linkage':
                if not g.attrib.get('data.parentResidueIndex') or not g.attrib.get('data.childResidueIndex'):
                    continue
                try:
                    parent = monos[int(g.attrib['data.parentResidueIndex'])]
                except KeyError:
                    raise GlycanBuilderSVGUnexpectedLinkConnectionError()
                if g.attrib['data.parentPositions'] == "?":
                    parentpos = None
                else:
                    try:
                        parentpos = set(map(int,g.attrib['data.parentPositions'].split('/')))
                    except ValueError:
                        raise GlycanBuilderSVGUnexpectedLinkConnectionError()
                try:
                    child = monos[int(g.attrib['data.childResidueIndex'])]
                except KeyError:
                    raise GlycanBuilderSVGUnexpectedLinkConnectionError()
                if g.attrib['data.childPositions'] == "?":
                    childpos = None
                else:
                    childpos = set(map(int,g.attrib['data.childPositions'].split('/')))
                parent.add_child(child, parent_pos=parentpos, child_pos=childpos,
                                        parent_type=Linkage.oxygenPreserved,
                                        child_type=Linkage.oxygenLost)

        coremonos = set(monos)
        undetroots = set()
        if len(undetroot) > 0:
            toexplore = set(undetroot)
            while len(toexplore) > 0:
                resid = toexplore.pop()
                if resid in coremonos:
                    coremonos.remove(resid)
                for c in monos[resid].children():
                    toexplore.add(c.id())
            for resid in undetroot:
                for i in range(undetroot[resid]):
                    if i > 0:
                        if monos[resid].is_monosaccharide():
                            m = monos[resid].deepclone()
                        else:
                            m = monos[resid].clone()
                    else:
                        m = monos[resid]
                    childpos,parentpos = undetrootlinkpos[resid]
                    for coreresid in coremonos:
                        linkage = UninstantiatedLinkage(child_pos=childpos,
                                                        parent_pos=parentpos,
                                                        parent_type=Linkage.oxygenPreserved,
                                                        child_type=Linkage.oxygenLost)
                        monos[coreresid].add_child_with_special_linkage(m,linkage)
                    m.set_connected(False)
                    undetroots.add(m)

        if len(coremonos) == 0:
            gly = Glycan()
        else:
            if rootmonoindex not in coremonos:
                raise GlycanBuilderSVGMissingRootError()
            rootmono = monos[rootmonoindex]
            gly = Glycan(rootmono)
        gly.set_undetermined(undetroots)
        gly.set_ids()
        return gly
