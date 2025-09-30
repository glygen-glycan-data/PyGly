#!/bin/env python3.12

import sys,inspect,csv
from operator import itemgetter
from collections import defaultdict
import findpygly
from pygly.Monosaccharide import SuperClass, Stem, Config, Mod, Substituent, Anomer, Linkage
from pygly.manipulation import Composition

import logging

class ClassifierEngine(object):

    def __init__(self,glytoucan=None,alignments=None,verbose=False):
        self._classifiers = []
        for n,c in  inspect.getmembers(sys.modules[__name__], inspect.isclass):
            if hasattr(c,'_class'):
                self._classifiers.append(c())
        self.log = logging.getLogger('ClassifierEngine')
        if verbose:
            self.log.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(name)s.%(funcName)s:%(lineno)s: %(message)s')
        handler = logging.StreamHandler()
        handler.setFormatter(formatter)
        self.log.addHandler(handler)        
        self._gtc = glytoucan
        self._al = alignments
        self.init(None)

    def init(self,acc):
        self._acc = acc
        self._motifs = None
        self._monocnt = None
        self._hasmono = None
        self._structure = None
        self._tests = dict()

    def __str__(self):
        retstr = []
        retstr.append("acc: "+str(self._acc))
        retstr.append("motifs: "+str(self._motifs))
        retstr.append("monocnt: "+str(self._monocnt))
        retstr.append("hasmono: "+str(self._hasmono))
        retstr.append("structure: "+str(self._structure))
        retstr.append("tests: "+str(self._tests))
        return "\n".join(retstr)

    def getmotifs(self):
        if self._acc not in self._al:
            self._motifs = set()
        else:
            self._motifs = set(self._al[self._acc]['Loose'])
        self.log.debug("Motifs: %s",self._motifs)

    def gethasmono(self):
        self._hasmono = set()
        if not self._monocnt:
            self.getmonocnt()
        for key in self._monocnt:
            if key != 'Monosaccharide':
                self._hasmono.add(key)
        self.log.debug("HasMono: %s",self._hasmono)

    def getmonocnt(self):
        self._monocnt = defaultdict(int)
        if self._structure == None:
            self.getStructure()
        if self._structure:
            if not self._structure.repeated():
                comp = self._structure.iupac_composition()
                comp1 = self._structure.iupac_composition(floating_substituents=False)
                comp2 = defaultdict(int)
                comp3 = defaultdict(int)
            else:
                comp = self._structure.iupac_composition(repeat_times=1)
                comp1 = self._structure.iupac_composition(repeat_times=1,floating_substituents=False)
                comp2 = self._structure.iupac_composition(repeat_times=2)
                comp3 = self._structure.iupac_composition(repeat_times=2,floating_substituents=False)
            for key in comp2:
                comp2[key] -= comp[key]
            for key in comp3:
                comp3[key] -= comp1[key]            
            # print(comp,comp2)
            for ckey,count in comp.items():
                if count > 0 and not ckey.startswith('_'):
                    count1 = (count if comp2[ckey] == 0 else str(count)+"+")
                    if ckey=='Count':
                        self._monocnt['Monosaccharide'] = count1
                    else:
                        self._monocnt[ckey] = count1
            for ckey,count in comp1.items():
                if count > 0 and not ckey.startswith('_') and ckey != "Count":
                    count1 = (count if comp3[ckey] == 0 else str(count)+"+")
                    if ckey.endswith('+aldi'):
                        self._monocnt[ckey] = count1
        else:
            raise RuntimeError("No source for monosaccharide counts for accession: %s"%(self._acc,))
        self.log.debug("MonoComp: %s",self._monocnt)

    def getStructure(self):
        glycan = None
        self._structure = self._gtc.getGlycan(self._acc,format='wurcs')
        if self._structure == None:
            raise RuntimeError("No source for structure for accession: %s"%(self._acc,))
        try:
            self.log.debug("Structure: %s"%(self._structure,))
        except AssertionError:
            self.log.debug("Structure: %r"%(self._structure,))
        except KeyError:
            self.log.debug("Structure: %r"%(self._structure,))

    def classifiers(self):
        return self._classifiers

    def assign(self,acc):
        self.init(acc)
        cls = set()
        for c in self._classifiers:
            if c.match(self):
                cls.add(c.classification())
        cls = sorted(cls)
        retain = set()
        for i in range(len(cls)):
            suppress = False
            for j in range(i+1,len(cls)):
                assert len(cls[i]) <= len(cls[j])
                if cls[i] == cls[j][:len(cls[i])]:
                    suppress = True
                    break
            if not suppress:
                retain.add(cls[i])
        return retain

    def has_motif(self,m):
        if self._motifs == None:
            self.getmotifs()
        return m in self._motifs

    def has_mono(self,m):
        if self._hasmono == None:
            self.gethasmono()
        return m in self._hasmono

    def all_mono(self,*args):
        total = 0
        for a in list(args) + ['Monosaccharide',]:
            mcnt = self.mono_count(a)
            try:
                val = int(mcnt)
            except (TypeError,ValueError):
                val = int(mcnt[:-1])
            if a == 'Monosaccharide':
                total -= val
            else:
                total += val
        return (total == 0)

    def mono_count(self,*args): 
        if len(args) == 0:
            args = ["Monosaccharide"]
        if self._monocnt == None:
            self.getmonocnt()
        if len(args) == 1:
            return self._monocnt[args[0]]
        return sum(map(lambda m: self._monocnt[m],args))

    def has_repeat_units(self):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        return self._structure.repeated()

    class Decorators:
        @staticmethod
        def maketest(key,fn,*args,**kw):
            def wrapped(self):
                if key not in self._tests:
                    self._tests[key] = fn(self,*args,**kw)
                return self._tests[key]
            return wrapped

        @staticmethod
        def needs_structure(f):
            def wrapped(self,*args,**kw):
                if self._structure == None:
                    self.getStructure()
                assert(self._structure)
                return f(self,*args,**kw)
            return wrapped

    def redend(self,anomer):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        r = self._structure.root()
        if r.anomer() == anomer:
            return True
        return False

    redendalpha = Decorators.maketest('redendalpha',redend,Anomer.alpha)
    redendbeta = Decorators.maketest('redendbeta',redend,Anomer.beta)

    class MonoMethods:

        @staticmethod
        def subst_count(m):
            cnt = defaultdict(int)
            for s in m.substituents():
                cnt[s.name()] += 1
                cnt[None] += 1
            return cnt

        @staticmethod
        def mod_count(m):
            cnt = defaultdict(int)
            for mod in m.mods():
                if len(mod[0]) == 1:
                    cnt[(mod[0][0],mod[1])] += 1
                    cnt[mod[1]] += 1
                else:
                    cnt[mod] +=1
                    cnt[mod[1]] += 1
                cnt[None] += 1
            return cnt

        @staticmethod
        def const_match(m,fn,const):
            cnt = fn(m)
            total = 0
            for key,value in const.items():
                if isinstance(value,int):
                    if cnt[key] != value:
                        return False
                    total += cnt[key]
                elif value.startswith('>='):
                    if cnt[key] < int(value[2:]):
                        return False
                    total += cnt[key]
                elif value.startswith('>'):
                    if cnt[key] <= int(value[1:]):
                        return False
                    total += cnt[key]
            if cnt[None] != total:
                return False
            return True
        
        @staticmethod
        def mod_match(m,const):
            return ClassifierEngine.MonoMethods.const_match(m,ClassifierEngine.MonoMethods.mod_count,const)

        @staticmethod
        def subst_match(m,const):
            return ClassifierEngine.MonoMethods.const_match(m,ClassifierEngine.MonoMethods.subst_count,const)

        @staticmethod
        def floating_subst_count(m):
            cnt = 0
            for s in m.substituents():
                if s in Composition.floating_substs:
                    cnt += 1
            return cnt


        @staticmethod
        def hex(m):
            if m.superclass() == SuperClass.HEX and m.stem() == Stem.missing:
                if not m.has_mods() and not m.has_substituents():
                    return True
            return False

        @staticmethod
        def glc(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if not m.has_mods() and not m.has_substituents():
                    return True
            return False

        @staticmethod
        def glca(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if m.has_substituents():
                    return False
                if ClassifierEngine.MonoMethods.mod_match(m,{(6,Mod.a): 1}):
                    return True
            return False

        @staticmethod
        def sulfglc(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if m.has_mods():
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0"}):
                    return True
            return False

        @staticmethod
        def sulfglca(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(6,Mod.a): 1}):
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0"}):
                    return True
            return False

        @staticmethod
        def methglca(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(6,Mod.a): 1}):
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.methyl: ">0"}):
                    return True
            return False

        @staticmethod
        def nsulfglc(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if m.has_mods():
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">=0",Substituent.methyl: ">=0",Substituent.nsulfate: 1}):
                    return True
            return False    

        @staticmethod
        def aldinsulfglc(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(1,Mod.aldi): 1}):
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">=0",Substituent.nsulfate: 1}):
                    return True
            return False

        @staticmethod
        def aminoglc(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if m.has_mods():
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">=0",Substituent.amino: 1}):
                    return True
            return False

        @staticmethod
        def aldiglcnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(1,Mod.aldi): 1}):
                    return False
                if m.is_nacetylated():
                    return True
            return False

        @staticmethod
        def sulfaldiglcnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(1,Mod.aldi): 1}):
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0",Substituent.nAcetyl: 1}):
                    return True
            return False

        @staticmethod
        def gal(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if not m.has_mods() and not m.has_substituents():
                    return True
            return False

        @staticmethod
        def sulfgal(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if m.has_mods():
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0"}):
                    return True
            return False

        @staticmethod
        def nsulfgal(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if m.has_mods():
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">=0",Substituent.nsulfate: 1}):
                    return True
            return False

        @staticmethod
        def aldigalnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(1,Mod.aldi): 1}):
                    return False
                if m.is_nacetylated():
                    return True
            return False

        @staticmethod
        def sulfaldigalnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(1,Mod.aldi): 1}):
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0",Substituent.nAcetyl: 1}):
                    return True
            return False

        @staticmethod
        def man(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.man,):
                if not m.has_mods() and not m.has_substituents():
                    return True
            return False

        @staticmethod
        def hexnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == Stem.missing:
                if not m.has_mods() and m.is_nacetylated():
                    return True
            return False

        @staticmethod
        def glcnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if not m.has_mods() and m.is_nacetylated():
                    return True
            return False

        @staticmethod
        def sulfglcnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.glc,):
                if m.has_mods():
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0",Substituent.nAcetyl: 1}):
                    return True
            return False

        @staticmethod
        def galnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if not m.has_mods() and m.is_nacetylated():
                    return True
            return False

        @staticmethod
        def sulfgalnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if m.has_mods():
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0",Substituent.nAcetyl: 1}):
                    return True
            return False

        @staticmethod
        def fuc(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,) and m.config() == (Config.l,):
                if m.has_substituents():
                    return False
                if m.count_mod(Mod.d) == 1 and m.count_mod() == 1:
                    return True
            return False

        @staticmethod
        def xyl(m):
            if m.superclass() == SuperClass.PENT and m.stem() == (Stem.xyl,) and m.config() == (Config.d,):
                if m.has_substituents():
                    return False
                if m.has_mods():
                    return False
                return True
            return False

        @staticmethod
        def hexa(m):
            if m.superclass() == SuperClass.HEX and m.stem() == Stem.missing:
                if m.has_substituents():
                    return False
                if ClassifierEngine.MonoMethods.mod_match(m,{(6,Mod.a): 1}):
                    return True
            return False

        @staticmethod
        def idoa(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.ido,) and m.config() == (Config.l,):
                if m.has_substituents():
                    return False
                if ClassifierEngine.MonoMethods.mod_match(m,{(6,Mod.a): 1}):
                    return True
            return False

        @staticmethod
        def sulfidoa(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.ido,) and m.config() == (Config.l,):
                if not ClassifierEngine.MonoMethods.mod_match(m,{(6,Mod.a): 1}):
                    return False
                if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0",Substituent.methyl: ">=0"}):
                    return True
            return False

        @staticmethod
        def sia(m):
            if m.superclass() == SuperClass.NON and m.stem() == (Stem.gro,Stem.gal):
                if set(map(itemgetter(1),m.mods())) == set([Mod.a,Mod.keto,Mod.d]):
                    found = False
                    for s in m.substituents():
                        if s.name() in (Substituent.nAcetyl,Substituent.nglycolyl):
                            found = True
                            break
                    if len(list(m.substituents())) == 1 and found:
                        return True
            return False

        @staticmethod
        def sulfsia(m):
            if m.superclass() == SuperClass.NON and m.stem() == (Stem.gro,Stem.gal):
                if set(map(itemgetter(1),m.mods())) == set([Mod.a,Mod.keto,Mod.d]):
                    found = False
                    for s in m.substituents():
                        if s.name() in (Substituent.nAcetyl,Substituent.nglycolyl):
                            found = True
                            break
                    if ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">0", Substituent.nAcetyl: ">=0", Substituent.nglycolyl: ">=0"}) and found:
                        return True
            return False

        @staticmethod
        def neuac(m):
            if m.superclass() == SuperClass.NON and m.stem() == (Stem.gro,Stem.gal):
                if set(map(itemgetter(1),m.mods())) == set([Mod.a,Mod.keto,Mod.d]):
                    found = False
                    for s in m.substituents():
                        if s.name() in (Substituent.nAcetyl,):
                            found = True
                            break
                    if len(list(m.substituents())) == 1 and found:
                        return True
            return False

        @staticmethod
        def hexenxa(m):
            if m.superclass() == SuperClass.HEX: # and m.stem() == (Stem.thr,):
                if not ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">=0"}):
                    return False
                if ClassifierEngine.MonoMethods.mod_match(m,{((4,5),Mod.enx): 1,
                                                             (6,Mod.a): 1,
                                                             Mod.d: ">=0"}):
                    return True
                elif ClassifierEngine.MonoMethods.mod_match(m,{((4,5),Mod.en): 1,
                                                               (6,Mod.a): 1,
                                                               Mod.d: ">=0"}):
                    return True
            return False

        @staticmethod
        def threnxa(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.thr,):
                if not ClassifierEngine.MonoMethods.subst_match(m,{Substituent.sulfate: ">=0"}):
                    return False
                if ClassifierEngine.MonoMethods.mod_match(m,{((4,5),Mod.enx): 1,
                                                             (6,Mod.a): 1,
                                                             Mod.d: ">=0"}):
                    return True
                elif ClassifierEngine.MonoMethods.mod_match(m,{((4,5),Mod.en): 1,
                                                               (6,Mod.a): 1,
                                                               Mod.d: ">=0"}):
                    return True
            return False

    def acceptcore(self,*acceptedfns):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        r = self._structure.root()
        for c in r.children():
            okseen = False
            for fn in acceptedfns:
                if fn(c):
                    okseen = True
                    break
            if okseen:
                continue
            if self.MonoMethods.fuc(c) or self.MonoMethods.sia(c) or self.MonoMethods.sulfsia(c):
                continue
            return False
        return True

    acceptcore1 = Decorators.maketest('acceptcore1',acceptcore,MonoMethods.gal,MonoMethods.sulfgal)
    acceptcore36 = Decorators.maketest('acceptcore36',acceptcore,MonoMethods.glcnac,MonoMethods.sulfglcnac)
    acceptcore57 = Decorators.maketest('acceptcore57',acceptcore,MonoMethods.galnac,MonoMethods.sulfgalnac)
    acceptcore8 = Decorators.maketest('acceptcore8',acceptcore,MonoMethods.gal,MonoMethods.sulfgal)
    acceptcore9 = Decorators.maketest('acceptcore9',acceptcore,MonoMethods.gal,MonoMethods.sulfgal)

    def acceptbarecore(self,*acceptedfns):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        r = self._structure.root()
        if self.MonoMethods.floating_subst_count(r) > 0:
            return False
        for c in r.children():
            okseen = False
            for fn in acceptedfns:
                if fn(c):
                    okseen = True
                    break
            if okseen:
                continue
            return False
        return True

    acceptoglc = Decorators.maketest('acceptoglc',acceptbarecore,MonoMethods.xyl)
    acceptegfoglcnac = Decorators.maketest('acceptegfoglcnac',acceptbarecore,MonoMethods.gal)

    def testncore(self):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        r = self._structure.root()
        if not self.MonoMethods.glcnac(r) and not self.MonoMethods.hexnac(r):
            return False
        ncores = []
        for l in r.links():
            c = l.child()
            if c.anomer() != Anomer.alpha and (l.parent_pos() == None or 4 in l.parent_pos()) and (self.MonoMethods.glcnac(c) or self.MonoMethods.hexnac(c)):
                ncores.append((r,c))
        if len(ncores) == 0:
            return False
        for i in range(len(ncores)):
            nc = ncores.pop(0)
            for l in nc[-1].links():
                c = l.child()
                if c.anomer() != Anomer.alpha and (l.parent_pos() == None or 4 in l.parent_pos()) and (self.MonoMethods.man(c) or self.MonoMethods.hex(c)):
                    ncores.append(tuple(list(nc)+[c]))
        if len(ncores) == 0:
            return False
        for i in range(len(ncores)):
            nc = ncores.pop(0)
            n3man = 0; n6man = 0; n36man = 0;
            for l in nc[-1].links():
                c = l.child()
                if c.anomer() != Anomer.beta and (l.parent_pos() == None or 3 in l.parent_pos()) and (self.MonoMethods.man(c) or self.MonoMethods.hex(c)):
                    n3man += 1
                if c.anomer() != Anomer.beta and (l.parent_pos() == None or 6 in l.parent_pos()) and (self.MonoMethods.man(c) or self.MonoMethods.hex(c)):
                    n6man += 1
                if c.anomer() != Anomer.beta and (l.parent_pos() == None or 3 in l.parent_pos() or 6 in l.parent_pos()) and (self.MonoMethods.man(c) or self.MonoMethods.hex(c)):
                    n36man += 1
            if n36man >= 2 and n3man >= 1 and n6man >= 1:
                ncores.append(nc)
        if len(ncores) == 0:
            return False
        return True

    def count_mono(self,*fnargs):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        cnt = 0
        for m in self._structure.all_nodes():
            for fn in fnargs:
                if fn(m):
                    cnt += 1
                    break
        return cnt

    def is_leaf(self,fn):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        cnt = 0
        for m in self._structure.all_nodes():
            if fn(m):
                if m.has_links(default=False):
                    return False
        return True

    countsulfgals = Decorators.maketest('countsulfgals',count_mono,MonoMethods.sulfgal)
    countbaregals = Decorators.maketest('countbaregals',count_mono,MonoMethods.gal)

    countglc = Decorators.maketest('countglc',count_mono,MonoMethods.glc)
    countglca = Decorators.maketest('countglca',count_mono,MonoMethods.glca)
    countsulfglcnac = Decorators.maketest('countsulfglcnac',count_mono,MonoMethods.sulfglcnac)
    countnsulfglc = Decorators.maketest('countnsulfglc',count_mono,MonoMethods.nsulfglc)
    countbareglcnac = Decorators.maketest('countbareglcnac',count_mono,MonoMethods.glcnac)
    countoptsulfglcnac = Decorators.maketest('countoptsulfglcnac',count_mono,MonoMethods.glcnac,MonoMethods.sulfglcnac,MonoMethods.nsulfglc)
    
    countsulfgalnac = Decorators.maketest('countsulfgalnac',count_mono,MonoMethods.sulfgalnac)
    countnsulfgal = Decorators.maketest('countnsulfgal',count_mono,MonoMethods.nsulfgal)
    countbaregalnac = Decorators.maketest('countbaregalnac',count_mono,MonoMethods.galnac)
    countoptsulfgalnac = Decorators.maketest('countoptsulfgalnac',count_mono,MonoMethods.galnac,MonoMethods.sulfgalnac,MonoMethods.nsulfgal)

    countthrenxa = Decorators.maketest('countthrenxa',count_mono,MonoMethods.threnxa)
    countneuac = Decorators.maketest('countneuac',count_mono,MonoMethods.neuac)
    
    countvalidxxx = Decorators.maketest('countvalidxxx',count_mono,MonoMethods.threnxa,MonoMethods.nsulfglc,MonoMethods.nsulfgal)

    countvalidgag = Decorators.maketest('countvalidgag',count_mono,
                                        MonoMethods.glca,MonoMethods.sulfglca,MonoMethods.methglca,
                                        MonoMethods.aminoglc,MonoMethods.nsulfglc,MonoMethods.aldinsulfglc,
                                        MonoMethods.gal,
                                        MonoMethods.sulfgal,MonoMethods.nsulfgal,
                                        MonoMethods.glcnac,MonoMethods.sulfglcnac,MonoMethods.aldiglcnac,MonoMethods.sulfaldiglcnac,
                                        MonoMethods.galnac,MonoMethods.sulfgalnac,MonoMethods.aldigalnac,MonoMethods.sulfaldigalnac,
                                        MonoMethods.idoa,MonoMethods.sulfidoa,
                                        MonoMethods.hexa,
                                        MonoMethods.hexenxa,
                                        MonoMethods.neuac,MonoMethods.fuc)

    threnxa_is_leaf = Decorators.maketest('threnxa_is_leaf',is_leaf,MonoMethods.threnxa)
    neuac_is_leaf = Decorators.maketest('neuac_is_leaf',is_leaf,MonoMethods.neuac)

    @Decorators.needs_structure
    def nobranch(self):
        for m in self._structure.all_nodes():
            if m.link_count(default=False,repeat=-Linkage.REPEAT_EXIT) > 1:
                return False
        return True

    @Decorators.needs_structure
    def nononfucbranch(self):
        for m in self._structure.all_nodes():
            if m.link_count(default=False,repeat=-Linkage.REPEAT_EXIT) > 1:
                lcnt = m.link_count(default=False,repeat=-Linkage.REPEAT_EXIT)
                for l in m.links(default=False,repeat=Linkage.NON_REPEAT):
                    if ClassifierEngine.MonoMethods.fuc(l.child()):
                        lcnt -= 1
                if lcnt > 1:
                    return False
        return True

    def repeat_unit_sizes_are_even(self):
        if self._structure == None:
            self.getStructure()
        if not self._structure:
            return False
        for nds in self._structure.repeat_nodes():
            if len(list(nds)) % 2 != 0:
                return False
        return True

class Classifier(object):
    _description = ""

    def match(self,data):
        raise NotImplemented()

    def classification(self):
        return tuple(self._class)

class MotifClassifier(Classifier):
    _exceptions = []

    def classification(self):
        return self._id,tuple(self._class),self._goodmotif

    def match(self,data):
        self._goodmotif = None
        good = False
        for m in self._motifs:
            if data.has_motif(m):
                good = True
                self._goodmotif = m
                break
        if good:
            for m in self._exceptions:
                if data.has_motif(m):
                    good = False
                    break
        if good:
            return self.refine(data)
        return False

    def refine(self,data):
        return True

class NGlycanBase(MotifClassifier):
    _id = "GGC.000001"
    _description = """
       N-linked glycan type determined by core alignment of five residue N-Glycan reducing-end motif.
    """.strip()
    _class = ("N-linked",)
    _motifs = ["GGM.001001"]

class NGlycanSuper(MotifClassifier):
    _id = "GGC.000002"
    _description = """
       N-linked glycan type determined by core alignment of HexNAc(2)Hex(3) reducing-end motif. Structures with a core alignment to the N-Glycan core basic motif (GGM.001001) are excluded, these are matched by GGC.000001. Additional checks ensure carbon bonds, anomeric configurations, and specific monosaccharides are consistent with the N-glycan core motif.
    """.strip()
    _class = ("N-linked",)
    _motifs = ["GGM.001027"]
    _exceptions = ["GGM.001001"]

    def refine(self,data):
        return data.testncore()

class NGlycanHybrid(MotifClassifier):
    _id = "GGC.000003"
    _description = """
       Hybrid N-linked glycan subtype determined by core alignment with two antennae with GlcNAc and Man residues. Only valid carbon bond positions are enumerated: 2 or 4 on the alpha-3 Man arm for GlcNAc; and 3 or 6 on the alpha-6 Man arm for Man. 
    """.strip()
    _class = ("N-linked","N-linked hybrid")
    _motifs = ["GGM.001003"]

class NGlycanComplex(MotifClassifier):
    _id = "GGC.000004"
    _description = """
       Complex N-linked glycan subtype determined by core alignment with two antennae GlcNAc residues. Only valid carbon bond positions are enumerated: 2 or 4 on the alpha-3 Man arm; and 2, 4 or 6 on the alpha-6 Man arm. 
    """.strip()
    _class = ("N-linked","N-linked complex")
    _motifs = ["GGM.001004"]

class NGlycanHighMannose(MotifClassifier):
    _id = "GGC.000005"
    _description = """
       High-mannose N-linked glycan subtype determined by core alignment to the N-Glycan core basic motif (GGM.001001), plus checks to ensure all residues, except for two reducing end GlcNAc, are Man. The additional motif (GGM.001002) is included for historical reasons, but does not result in any additional structures. 
    """.strip()

    _class = ("N-linked","N-linked high mannose")
    _motifs = ["GGM.001002","GGM.001001"]

    def refine(self,data):
        if data.has_repeat_units():
            return False
        if data.mono_count("Man") + 2 != data.mono_count():
            return False
        return True

class NGlycanPaucimannose(MotifClassifier):
    _id = "GGC.000006"
    _class = ("N-linked","N-linked paucimannose")
    _motifs = ["GGM.001001"]
        
    def refine(self,data):
        if data.has_repeat_units():
            return False
        if data.mono_count() != data.mono_count("Man","GlcNAc","Fuc"):
            return False
        if data.mono_count('GlcNAc') != 2:
            return False
        if data.mono_count('Man') != 3:
            return False
        if data.mono_count('Fuc') != (data.mono_count() - 5):
            return False
        return True
    
class NGlycanTruncated1(MotifClassifier):
    _id = "GGC.000007"
    _class = ("N-linked","N-linked truncated")
    _motifs = ["GGM.001022"]

class NGlycanTruncated2(MotifClassifier):
    _id = "GGC.000008"
    _class = ("N-linked","N-linked monoantennary")
    _motifs = ["GGM.001041"]
    _exceptions = ["GGM.001027"]

class NGlycanTruncated3(MotifClassifier):
    _id = "GGC.000009"
    _class = ("N-linked","N-linked truncated")
    _motifs = ["GGM.001046"]

class NGlycanBisected(MotifClassifier):
    _id = "GGC.000010"
    _class = ("N-linked","N-linked bisected")
    _motifs = ["GGM.001023"]

class NGlycanCoreFucosylated(MotifClassifier):
    _id = "GGC.000011"
    _class = ("N-linked","N-linked core-fucosylated")
    _motifs = ["GGM.001024"]

class NGlycanArmFucosylated(MotifClassifier):
    _id = "GGC.000012"
    _class = ("N-linked","N-linked arm-fucosylated")
    _motifs = ["GGM.001025"]

class NGlycanAlditolReduced(MotifClassifier):
    _id = "GGC.000013"
    _class = ("N-linked","N-linked alditol-reduced")
    _motifs = ["GGM.001026"]

class NGlycanBiantennary(MotifClassifier):
    _id = "GGC.000014"
    _class = ("N-linked","N-linked biantennary")
    _motifs = ["GGM.001037"]
    _exceptions = ["GGM.001038"]

    def refine(self,data):
        if not data.has_motif("GGM.001004"):
            return False
        return True

class NGlycanTriantennary(MotifClassifier):
    _id = "GGC.000015"
    _class = ("N-linked","N-linked triantennary")
    _motifs = ["GGM.001038"]
    _exceptions = ["GGM.001039"]

    def refine(self,data):
        if not data.has_motif("GGM.001004"):
            return False
        return True

class NGlycanTetraantennary(MotifClassifier):
    _id = "GGC.000016"
    _class = ("N-linked","N-linked tetraantennary")
    _motifs = ["GGM.001039"]

    def refine(self,data):
        if not data.has_motif("GGM.001004"):
            return False
        return True

class OGlycanOGlcNAc(MotifClassifier):
    _id = "GGC.000101"
    _class = ("O-linked","O-GlcNAc")
    _motifs = ["GGM.001028"]

class OGlycanEGFOGlcNAcCore(MotifClassifier):
    _id = "GGC.000102"
    _class = ("O-linked","EGF type O-GlcNAc core")
    _motifs = ["GGM.001040"]

    def refine(self,data):
        if not data.acceptegfoglcnac():
            return False
        if data.redendalpha():
            return False
        return True

class OGlycanHMOligos(MotifClassifier):
    _id = "GGC.000103"
    _class = ("Human Milk Oligosaccharide",)
    _motifs = ["GGM.001029"]
    _exceptions = ["GGM.001101","GGM.001102","GGM.001108","GGM.001109"]

class OGlycanOGal(MotifClassifier):
    _id = "GGC.000104"
    _class = ("O-linked","O-galactose")
    _motifs = ["GGM.001044"]

class OGlycanOGalNAc(MotifClassifier):
    _id = "GGC.000105"
    _class = ("O-linked","O-GalNAc core, other")
    _motifs = ["GGM.001034"]
    _exceptions = ["GGM.001005","GGM.001006","GGM.001007","GGM.001008","GGM.001009","GGM.001010",
                   "GGM.001011","GGM.001012","GGM.001013","GGM.001014","GGM.001015","GGM.001016",
                   "GGM.001017","GGM.001018","GGM.001032","GGM.001033","GGM.001042","GGM.001043",
                   "GGM.001203","GGM.001207"]

class OGlycanOGlucose(MotifClassifier):
    _id = "GGC.000106"
    _class = ("O-linked","O-glucose core")
    _motifs = ["GGM.001031"]

    def refine(self,data):
        if not data.acceptoglc():
            return False
        return True

class OGlycanOMannose1(MotifClassifier):
    _id = "GGC.000107"
    _class = ("O-linked","O-mannose core")
    _motifs = ["GGM.001019"]

class OGlycanOMannose3(MotifClassifier):
    _id = "GGC.000109"
    _class = ("O-linked","O-mannose")
    _motifs = ["GGM.001036"]

    def refine(self,data):
        if data.has_repeat_units():
            return False
        if data.mono_count("Man") != data.mono_count():
            return False
        return True

class OGlycanOFucose1(MotifClassifier):
    _id = "GGC.000110"
    _class = ("O-linked","O-fucose core")
    _motifs = ["GGM.001020","GGM.001021"]

    def refine(self,data):
        if data.redendbeta():
            return False
        return True

class OGlycanOFucose2(MotifClassifier):
    _id = "GGC.000111"
    _class = ("O-linked","O-fucose core")
    _motifs = ["GGM.001035"]

    def refine(self,data):
        if data.redendbeta():
            return False
        return True

class OGlycanCore1(MotifClassifier):
    _id = "GGC.000112"
    _class = ("O-linked","O-linked core 1")
    _motifs = ["GGM.001005","GGM.001006"]
    # unless core 2 or core 9 motifs
    _exceptions = ["GGM.001007","GGM.001008","GGM.001042","GGM.001043",]

    def refine(self,data):
        if not data.acceptcore1():
            return False
        return True

class OGlycanCore2(MotifClassifier):
    _id = "GGC.000113"
    _class = ("O-linked","O-linked core 2")
    _motifs = ["GGM.001007","GGM.001008"]

class OGlycanCore3(MotifClassifier):
    _id = "GGC.000114"
    _class = ("O-linked","O-linked core 3")
    _motifs = ["GGM.001009","GGM.001010"]
    # unless core 2 core 4 motifs
    _exceptions = ["GGM.001007","GGM.001008","GGM.001011","GGM.001012"]

    def refine(self,data):
        if not data.acceptcore36():
            return False
        return True

class OGlycanCore4(MotifClassifier):
    _id = "GGC.000115"
    _class = ("O-linked","O-linked core 4")
    _motifs = ["GGM.001011","GGM.001012"]

class OGlycanCore5(MotifClassifier):
    _id = "GGC.000116"
    _class = ("O-linked","O-linked core 5")
    _motifs = ["GGM.001013","GGM.001014"]

    def refine(self,data):
        if not data.acceptcore57():
            return False
        return True

class OGlycanCore6(MotifClassifier):
    _id = "GGC.000117"
    _class = ("O-linked","O-linked core 6")
    _motifs = ["GGM.001015","GGM.001016"]
    # unless core 2,4 motifs
    _exceptions = ["GGM.001007","GGM.001008","GGM.001011","GGM.001012"]

    def refine(self,data):
        if not data.acceptcore36():
            return False
        return True

class OGlycanCore7(MotifClassifier):
    _id = "GGC.000118"
    _class = ("O-linked","O-linked core 7")
    _motifs = ["GGM.001017","GGM.001018"]

    def refine(self,data):
        if not data.acceptcore57():
            return False
        return True

class OGlycanCore8(MotifClassifier):
    _id = "GGC.000119"
    _class = ("O-linked","O-linked core 8")
    _motifs = ["GGM.001032","GGM.001033"]
    _exceptions = ["GGM.001042","GGM.001043","GGM.001007","GGM.001008"]

    def refine(self,data):
        if not data.acceptcore8():
            return False
        return True

class OGlycanCore9(MotifClassifier):
    _id = "GGC.000120"
    _class = ("O-linked","O-linked core 9")
    _motifs = ["GGM.001042","GGM.001043"]

    def refine(self,data):
        if not data.acceptcore9():
            return False
        return True

class GlycosphingolipidIsoglobo(MotifClassifier):
    _id = "GGC.000201"
    _class = ("Glycosphingolipid","Glycosphingolipid isoglobo series")
    _motifs = ["GGM.001101"]

class GlycosphingolipidMuco(MotifClassifier):
    _id = "GGC.000202"
    _class = ("Glycosphingolipid","Glycosphingolipid muco series")
    _motifs = ["GGM.001102"]

class GlycosphingolipidArthro(MotifClassifier):
    _id = "GGC.000203"
    _class = ("Glycosphingolipid","Glycosphingolipid arthro series")
    _motifs = ["GGM.001103"]

class GlycosphingolipidMollu(MotifClassifier):
    _id = "GGC.000204"
    _class = ("Glycosphingolipid","Glycosphingolipid mollu series")
    _motifs = ["GGM.001104"]

class GlycosphingolipidGala(MotifClassifier):
    _id = "GGC.000205"
    _class = ("Glycosphingolipid","Glycosphingolipid gala series")
    _motifs = ["GGM.001105"]

class GlycosphingolipidLacto(MotifClassifier):
    _id = "GGC.000206"
    _class = ("Glycosphingolipid","Glycosphingolipid lacto series")
    _motifs = ["GGM.001106"]

class GlycosphingolipidNeolacto(MotifClassifier):
    _id = "GGC.000207"
    _class = ("Glycosphingolipid","Glycosphingolipid neo-lacto series")
    _motifs = ["GGM.001107"]

class GlycosphingolipidGanglio(MotifClassifier):
    _id = "GGC.000208"
    _class = ("Glycosphingolipid","Glycosphingolipid ganglio series")
    _motifs = ["GGM.001108"]

class GlycosphingolipidGlobo(MotifClassifier):
    _id = "GGC.000209"
    _class = ("Glycosphingolipid","Glycosphingolipid globo series")
    _motifs = ["GGM.001109"]

class GPIAnchor(MotifClassifier):
    _id = "GGC.000301"
    _class = ("GPI anchor",)
    _motifs = ["GGM.001030"]

class OGlycanCMannosyl(MotifClassifier):
    _id = "GGC.000501"
    _class = ("C-linked","C-linked Man")
    _motifs = ["GGM.001045"]

class GAGBase(MotifClassifier):
    _id = "GGC.000401"
    _class = ("GAG",)
    _motifs = ["GGM.001201",
               "GGM.001202",
               "GGM.001203",
               "GGM.001204",
               "GGM.001205",
               "GGM.001206",
               "GGM.001207",
               "GGM.001208",
               "GGM.001211",
               "GGM.001212",
               "GGM.001213",
               "GGM.001214"]

    def baserefine(self,data):
        # These tests are true for all GAGs...?
        if data.has_mono('GalA') or data.has_mono("Man"):
            data.log.info("Has GalA or Man")
            return False
        if not data.has_repeat_units():
            if data.mono_count() != data.countvalidgag():
                data.log.info("Unexpected GAG monosaccharide (no repeat)")
                return False
        else:
            if int(data.mono_count()[:-1]) != data.countvalidgag():
                data.log.info("Unexpected GAG monosaccharide (w/ repeat)")
                return False
        if data.countthrenxa() == 1 and not data.threnxa_is_leaf():
            data.log.info("Thr-enx-a not a leaf")
            return False
        if data.countthrenxa() > 1:
            data.log.info("Multiple Thr-enx-a observed")
            return False
        if data.has_repeat_units() and not data.repeat_unit_sizes_are_even():
            data.log.info("Repeat unit sizes are not even")
            return False
        return True

    def refine(self,data):
        # These tests only for those GAGs not with an explicit subtype class
        if not self.baserefine(data):
            return False
        if data.countbaregals() > 0:
            data.log.info("Has bare Gal")
            return False
        if not data.nobranch():
            data.log.info("Branched glycan")
            return False
        if data.has_mono('Fuc') or data.has_mono('Sia'):
            data.log.info("Fuc or Sia monosaccharide")
            return False
        if data.countoptsulfglcnac() > 0 and data.countoptsulfgalnac() > 0:
            data.log.info("Mixture of GlcNAc- and GalNAc-related monosaccharides")
            return False
        return True

class XylGAGCore(GAGBase):
    _id = "GGC.000402"
    _class = ("GAG","Xylose Core GAG")
    _motifs = ["GGM.001215"]

    def refine(self,data):
        # no other checks?
        if not self.baserefine(data):
            return False
        if not data.nobranch():
            data.log.info("Branched glycan")
            return False
        return True

class KeratanSulfate(GAGBase):
    _id = "GGC.000403"
    _class = ("GAG","Keratan sulfate")
    _motifs = ["GGM.001209",
               "GGM.001210"]

    def refine(self,data):
        if not self.baserefine(data):
            return False
        if not data.all_mono('GlcNAc','Gal','NeuAc'):
            data.log.info("Unexpected Keratan sulfate monosaccharide")
            return False
        if data.countneuac() == 1 and not data.neuac_is_leaf():
            data.log.info("Non-terminal NeuAc")
            return False
        if data.countneuac() > 1:
            data.log.info("Too many NeuAc")
            return False
        if not data.nobranch():
            data.log.info("Branched glycan")
            return False
        return True

if __name__ == "__main__":

    from getwiki import GlycanData
    from pygly.GlycanResource import GlycoMotif, GlycoMotifNoPrefetch, GlycoMotifNoCache
    from pygly.GlycanResource import GlyTouCan, GlyTouCanNoPrefetch, GlyTouCanNoCache

    w = GlycanData(); gtc = GlyTouCan(); gm = GlycoMotif()

    motifacc = None
    badbymotif = None
    if sys.argv[1].startswith('GGM.'):
        classifier = ClassifierEngine(verbose=False,
                                      glycandata=w,
                                      glytoucan=gtc,
                                      glycomotif=gm)
        if ':' in sys.argv[1]:
            motifacc,startafter = sys.argv[1].split(':',1)
        else:
            motifacc = sys.argv[1]
            startafter = None
        collection,motifacc = motifacc.split('.',1)
        accs = sorted(set(map(itemgetter(0),classifier._gm.getstruct(collection,motifacc))))
        if startafter != None:
            accs = filter(lambda acc: acc >= startafter, accs)

        if len(sys.argv) >= 3:
            badbymotif = defaultdict(set)
            for row in csv.DictReader(open(sys.argv[2]),dialect='excel-tab'):
                badbymotif[row['motif']].add(row['gtcacc'])
    else:
        classifier = ClassifierEngine(verbose=True,
                                      glycandata=w,
                                      glytoucan=gtc,
                                      glycomotif=gm)
        accs = sys.argv[1:]
    for acc in accs:
        any = False
        for asn in classifier.assign(acc):
            if not badbymotif:
                if not any:
                    print(acc)
                print(">>", asn[0], asn[1])
            else:
                print("\t".join([acc,motifacc,asn[0],asn[1],"FP" if (acc in badbymotif[motifacc]) else "TP"]))
            any = True
        if not any:
            if not badbymotif:
                print(acc)
                print(">>","No match")
            else:
                print("\t".join([acc,motifacc,"-","","FN" if (acc not in badbymotif[motifacc]) else "TN"]))
