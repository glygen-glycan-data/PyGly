#!/bin/env python2

import sys,inspect,csv
from operator import itemgetter
from collections import defaultdict
import findpygly
from pygly.Monosaccharide import SuperClass, Stem, Config, Mod, Substituent, Anomer, Linkage

import logging

class ClassifierEngine(object):

    def __init__(self,glycandata=None,glytoucan=None,glycomotif=None,verbose=False):
        self._classifiers = []
	for n,c in  inspect.getmembers(sys.modules[__name__], inspect.isclass):
	    if hasattr(c,'_class'):
		self._classifiers.append(c())
	self.log = logging.getLogger('ClassifierEngine')
	if verbose:
	    self.log.setLevel(logging.DEBUG)
	formatter = logging.Formatter('%(name)s.%(funcname)s:%(lineno)s: %(message)s')
        handler = logging.StreamHandler()
	handler.setFormatter(formatter)
        self.log.addHandler(handler)        
        self._gd = glycandata
	self._gtc = glytoucan
	self._gm = glycomotif
	self.init(None)

    def init(self,glycan):
	self._glycan = glycan
	self._motifs = None
	self._monocnt = None
	self._hasmono = None
	self._structure = None
	self._tests = dict()

    def getmotifs(self):
	self._motifs = set()
        glycan = None
        if self._gd:
            glycan = self._gd.get(self._glycan)
	if glycan:
	    for ann in glycan.annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
		for value in ann.get('value',[]):
		    self._motifs.add(value)
	elif self._gm:
	    for motifacc,strict in self._gm.getmotif(collection='GGM',accession=self._glycan):
		if int(motifacc.split('.',1)[1]) >= 1000:
		    self._motifs.add(motifacc)
        else:
            raise RuntimeError("No source for motif data for accession: %s"%(self._glycan,))
        self.log.debug("Motifs: %s",self._motifs)

    def gethasmono(self):
	self._hasmono = set()
        glycan = None
        if self._gd:
            glycan = self._gd.get(self._glycan)
	if glycan:
	    for m in glycan.get_annotation_values(property="HasMonosaccharide",source="EdwardsLab"):
		self._hasmono.add(m)
	elif self._gtc:
	    if not self._monocnt:
		self.getmonocnt()
	    for key in self._monocnt:
		if key != 'Monosaccharide':
		    self._hasmono.add(key)
        else:
            raise RuntimeError("No source for motif data for accession: %s"%(self._glycan,))
        self.log.debug("HasMono: %s",self._hasmono)

    def getmonocnt(self):
	self._monocnt = defaultdict(int)
        glycan = None
        if self._gd:
            glycan = self._gd.get(self._glycan)
	if glycan:
	    for ann in glycan.annotations(type="MonosaccharideCount",source="EdwardsLab"):
		prop = ann.get('property')
		try:
		    value = int(ann.get('value'))
		except ValueError:
		    value = ann.get('value')
		if prop.endswith('Count'):
		    prop = prop[:-5]
		self._monocnt[prop] = value
	elif self._gtc:
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
                raise RuntimeError("No source for monosaccharide counts for accession: %s"%(self._glycan,))
        else:
            raise RuntimeError("No source for monosaccharide counts for accession: %s"%(self._glycan,))
        self.log.debug("MonoComp: %s",self._monocnt)

    def getStructure(self):
        glycan = None
        if self._gd:
            glycan = self._gd.get(self._glycan)
	if glycan:
	    self._structure = glycan.getGlycan()
	elif self._gtc:
	    self._structure = self._gtc.getGlycan(self._glycan)
        else:
            raise RuntimeError("No source for structure for accession: %s"%(self._glycan,))
        if self._structure == None:
            raise RuntimeError("No source for structure for accession: %s"%(self._glycan,))
        # self.log.debug("Structure: %s"%(self._structure,))

    def assign(self,glycan):
	self.init(glycan)
	if not self._glycan:
	    return set()
	cls = set()
	for c in self._classifiers:
	    if c.match(self):
		cls.add(c.classification())
	# allow blank subtype only if no other classification with type
	bytype = defaultdict(int)
	for cl in cls:
	    if cl[1]:
		bytype[cl[0]] += 1
	retain = set()
	for cl in cls:
	    if cl[1]:
		retain.add(cl)
	    elif bytype[cl[0]] == 0:
		retain.add(cl)
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

    def acceptcore(self,acceptedfn):
	if self._structure == None:
	    self.getStructure()
	if not self._structure:
	    return False
	r = self._structure.root()
	okseen = False
	for c in r.children():
	    if acceptedfn(c) and not okseen:
		okseen = True
		continue
	    if self.MonoMethods.fuc(c) or self.MonoMethods.sia(c):
		continue
	    return False
	return True

    acceptcore1 = Decorators.maketest('acceptcore1',acceptcore,MonoMethods.gal)
    acceptcore36 = Decorators.maketest('acceptcore36',acceptcore,MonoMethods.glcnac)
    acceptcore57 = Decorators.maketest('acceptcore57',acceptcore,MonoMethods.galnac)

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
    countneuac = Decorators.maketest('countthrenxa',count_mono,MonoMethods.neuac)
    
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
	    if len(nds) % 2 != 0:
		return False
	return True

class Classifier(object):

    def match(self,data):
	raise NotImplemented()

	def classification(self):
	    return self._class

class MotifClassifier(Classifier):
    _exceptions = []

    def classification(self):
	return tuple(list(self._class)+[self._goodmotif])

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
    _class = ("N-linked","")
    _motifs = ["GGM.001001"]

class NGlycanSuper(MotifClassifier):
    _class = ("N-linked","")
    _motifs = ["GGM.001027"]
    _exceptions = ["GGM.001001"]

    def refine(self,data):
	return data.testncore()

class NGlycanHybrid(MotifClassifier):
    _class = ("N-linked","N-linked hybrid")
    _motifs = ["GGM.001003"]

class NGlycanComplex(MotifClassifier):
    _class = ("N-linked","N-linked complex")
    _motifs = ["GGM.001004"]

class NGlycanHighMannose(MotifClassifier):
    _class = ("N-linked","N-linked high mannose")
    _motifs = ["GGM.001002","GGM.001001"]

    def refine(self,data):
	if data.has_repeat_units():
	    return False
	if data.mono_count("Man") + 2 != data.mono_count():
	    return False
	return True

class NGlycanPaucimannose(MotifClassifier):
    _class = ("N-linked","N-linked paucimannose")
    _motifs = ["GGM.001001"]
	
    def refine(self,data):
	if data.has_repeat_units():
	    return False
	if data.mono_count() not in (5,6):
	    return False
	if data.mono_count() != data.mono_count("Man","GlcNAc","Fuc"):
	    return False
	if data.mono_count('GlcNAc') != 2:
	    return False
	if data.mono_count('Fuc') != (data.mono_count() - 5):
	    return False
        return True
    
class NGlycanTruncated(MotifClassifier):
    _class = ("N-linked","N-linked truncated")
    _motifs = ["GGM.001022"]

class NGlycanBisected(MotifClassifier):
    _class = ("N-linked","N-linked bisected")
    _motifs = ["GGM.001023"]

class NGlycanCoreFucosylated(MotifClassifier):
    _class = ("N-linked","N-linked core-fucosylated")
    _motifs = ["GGM.001024"]

class NGlycanArmFucosylated(MotifClassifier):
    _class = ("N-linked","N-linked arm-fucosylated")
    _motifs = ["GGM.001025"]

class NGlycanAlditolReduced(MotifClassifier):
    _class = ("N-linked","N-linked alditol-reduced")
    _motifs = ["GGM.001026"]

class NGlycanBiantennary(MotifClassifier):
    _class = ("N-linked","N-linked biantennary")
    _motifs = ["GGM.001037"]
    _exceptions = ["GGM.001038"]

    def refine(self,data):
        if not data.has_motif("GGM.001004"):
            return False
        return True

class NGlycanTriantennary(MotifClassifier):
    _class = ("N-linked","N-linked triantennary")
    _motifs = ["GGM.001038"]
    _exceptions = ["GGM.001039"]

    def refine(self,data):
        if not data.has_motif("GGM.001004"):
            return False
        return True

class NGlycanTetraantennary(MotifClassifier):
    _class = ("N-linked","N-linked tetraantennary")
    _motifs = ["GGM.001039"]

    def refine(self,data):
        if not data.has_motif("GGM.001004"):
            return False
        return True

class OGlycanOGlcNAc(MotifClassifier):
    _class = ("O-linked","O-GlcNAc")
    _motifs = ["GGM.001028"]

class OGlycanOMannose(MotifClassifier):
    _class = ("O-linked","O-mannose")
    _motifs = ["GGM.001019"]

class OGlycanOFucose(MotifClassifier):
    _class = ("O-linked","O-fucose")
    _motifs = ["GGM.001020","GGM.001021"]

    def refine(self,data):
	if data.redendbeta():
	    return False
	return True

class OGlycanCore1(MotifClassifier):
    _class = ("O-linked","O-linked core 1")
    _motifs = ["GGM.001005","GGM.001006"]
    # unless core 2 motifs
    _exceptions = ["GGM.001007","GGM.001008"]

    def refine(self,data):
	if not data.acceptcore1():
	    return False
        return True

class OGlycanCore2(MotifClassifier):
    _class = ("O-linked","O-linked core 2")
    _motifs = ["GGM.001007","GGM.001008"]

class OGlycanCore3(MotifClassifier):
    _class = ("O-linked","O-linked core 3")
    _motifs = ["GGM.001009","GGM.001010"]
    # unless core 2, core 4 motifs
    _exceptions = ["GGM.001007","GGM.001008","GGM.001011","GGM.001012"]

    def refine(self,data):
	if not data.acceptcore36():
	    return False
	return True

class OGlycanCore4(MotifClassifier):
    _class = ("O-linked","O-linked core 4")
    _motifs = ["GGM.001011","GGM.001012"]

class OGlycanCore5(MotifClassifier):
    _class = ("O-linked","O-linked core 5")
    _motifs = ["GGM.001013","GGM.001014"]

    def refine(self,data):
	if not data.acceptcore57():
	    return False
	return True

class OGlycanCore6(MotifClassifier):
    _class = ("O-linked","O-linked core 6")
    _motifs = ["GGM.001015","GGM.001016"]
    # unless core 2 motifs
    _exceptions = ["GGM.001007","GGM.001008"]

    def refine(self,data):
	if not data.acceptcore36():
	    return False
	return True

class OGlycanCore7(MotifClassifier):
    _class = ("O-linked","O-linked core 7")
    _motifs = ["GGM.001017","GGM.001018"]

    def refine(self,data):
	if not data.acceptcore57():
	    return False
	return True

class GlycosphingolipidIsoglobo(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid isoglobo series")
    _motifs = ["GGM.001101"]

class GlycosphingolipidMuco(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid muco series")
    _motifs = ["GGM.001102"]

class GlycosphingolipidArthro(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid arthro series")
    _motifs = ["GGM.001103"]

class GlycosphingolipidMollu(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid mollu series")
    _motifs = ["GGM.001104"]

class GlycosphingolipidGala(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid gala series")
    _motifs = ["GGM.001105"]

class GlycosphingolipidLacto(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid lacto series")
    _motifs = ["GGM.001106"]

class GlycosphingolipidNeolacto(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid neo-lacto series")
    _motifs = ["GGM.001107"]

class GlycosphingolipidGanglio(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid ganglio series")
    _motifs = ["GGM.001108"]

class GlycosphingolipidGlobo(MotifClassifier):
    _class = ("Glycosphingolipid","Glycosphingolipid globo series")
    _motifs = ["GGM.001109"]

class GPIAnchor(MotifClassifier):
    _class = ("GPI anchor","")
    _motifs = ["GGM.001030"]

class GAGBase(MotifClassifier):
    _class = ("GAG","")
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
	    data.log.info("Mixture of GlcNAc- and GalNac-related monosaccharides")
	    return False
        return True

class KeratanSulfate(GAGBase):
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
        classifier = ClassifierEngine(verbose=True)
	accs = sys.argv[1:]
    for acc in accs:
        any = False
	for asn in classifier.assign(acc):
            if not badbymotif:
                if not any:
                    print acc
                print ">>", asn[0], asn[1]
            else:
                print "\t".join([acc,motifacc,asn[0],asn[1],"FP" if (acc in badbymotif[motifacc]) else "TP"])
            any = True
        if not any:
            if not badbymotif:
                print acc
                print ">>","No match"
            else:
                print "\t".join([acc,motifacc,"-","","FN" if (acc not in badbymotif[motifacc]) else "TN"])
