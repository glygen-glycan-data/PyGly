#!/bin/env python2

import sys,inspect
from operator import itemgetter
from collections import defaultdict
import findpygly
from pygly.Monosaccharide import SuperClass, Stem, Config, Mod, Substituent, Anomer

class ClassifierEngine(object):

    def __init__(self):
        self._classifiers = []
	for n,c in  inspect.getmembers(sys.modules[__name__], inspect.isclass):
	    if hasattr(c,'_class'):
		self._classifiers.append(c())
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
	for ann in self._glycan.annotations(property='ClassMotif',type='Motif',source='GlycoMotif'):
	    for value in ann.get('value',[]):
		self._motifs.add(value)

    def gethasmono(self):
	self._hasmono = set()
	for m in self._glycan.get_annotation_values(type="HasMonosaccharide",source="EdwardsLab"):
	    self._hasmono.add(m)

    def getmonocnt(self):
	self._monocnt = defaultdict(int)
	for ann in self._glycan.annotations(type="MonosaccharideCount",source="EdwardsLab"):
	    prop = ann.get('property')
	    value = int(ann.get('value'))
	    if prop.endswith('Count'):
		prop = prop[:-5]
	    self._monocnt[prop] = value

    def getStructure(self):
	self._structure = self._glycan.getGlycan()
	if not self._structure:
	    self._structure = False

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

    def mono_count(self,*args): 
	if len(args) == 0:
	    args = ["Monosaccharide"]
	if self._monocnt == None:
	    self.getmonocnt()
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
        def hex(m):
            if m.superclass() == SuperClass.HEX and m.stem() == Stem.missing:
                if not m.has_mods() and not m.has_substituents():
                    return True
            return False
    
        @staticmethod
        def gal(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if not m.has_mods() and not m.has_substituents():
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
        def galnac(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,):
                if not m.has_mods() and m.is_nacetylated():
                    return True
            return False
    
        @staticmethod
        def fuc(m):
            if m.superclass() == SuperClass.HEX and m.stem() == (Stem.gal,) and m.config() == (Config.l,):
                if m.count_mod(Mod.d) == 1 and m.count_mod() == 1:
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
    _class = ("N-linked","hybrid")
    _motifs = ["GGM.001003"]

class NGlycanComplex(MotifClassifier):
    _class = ("N-linked","complex")
    _motifs = ["GGM.001004"]

class NGlycanHighMannose(MotifClassifier):
    _class = ("N-linked","high mannose")
    _motifs = ["GGM.001002","GGM.001001"]

    def refine(self,data):
	if data.has_repeat_units():
            return False
	if data.mono_count("Man") + 2 != data.mono_count():
	    return False
	return True

class NGlycanPaucimannose(MotifClassifier):
    _class = ("N-linked","paucimannose")
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
    _class = ("N-linked","truncated")
    _motifs = ["GGM.001022"]

class NGlycanBisected(MotifClassifier):
    _class = ("N-linked","bisected")
    _motifs = ["GGM.001023"]

class NGlycanCoreFucosylated(MotifClassifier):
    _class = ("N-linked","core-fucosylated")
    _motifs = ["GGM.001024"]

class NGlycanArmFucosylated(MotifClassifier):
    _class = ("N-linked","arm-fucosylated")
    _motifs = ["GGM.001025"]

class NGlycanAlditolReduced(MotifClassifier):
    _class = ("N-linked","alditol-reduced")
    _motifs = ["GGM.001026"]

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
    _class = ("O-linked","core 1")
    _motifs = ["GGM.001005","GGM.001006"]
    # unless core 2 motifs
    _exceptions = ["GGM.001007","GGM.001008"]

    def refine(self,data):
	if not data.acceptcore1():
	    return False
	return True

class OGlycanCore2(MotifClassifier):
    _class = ("O-linked","core 2")
    _motifs = ["GGM.001007","GGM.001008"]

class OGlycanCore3(MotifClassifier):
    _class = ("O-linked","core 3")
    _motifs = ["GGM.001009","GGM.001010"]
    # unless core 2, core 4 motifs
    _exceptions = ["GGM.001007","GGM.001008","GGM.001011","GGM.001012"]

    def refine(self,data):
	if not data.acceptcore36():
	    return False
	return True

class OGlycanCore4(MotifClassifier):
    _class = ("O-linked","core 4")
    _motifs = ["GGM.001011","GGM.001012"]

class OGlycanCore5(MotifClassifier):
    _class = ("O-linked","core 5")
    _motifs = ["GGM.001013","GGM.001014"]

    def refine(self,data):
	if not data.acceptcore57():
	    return False
	return True

class OGlycanCore6(MotifClassifier):
    _class = ("O-linked","core 6")
    _motifs = ["GGM.001015","GGM.001016"]
    # unless core 2 motifs
    _exceptions = ["GGM.001007","GGM.001008"]

    def refine(self,data):
	if not data.acceptcore36():
	    return False
	return True

class OGlycanCore7(MotifClassifier):
    _class = ("O-linked","core 7")
    _motifs = ["GGM.001017","GGM.001018"]

    def refine(self,data):
	if not data.acceptcore57():
	    return False
	return True

class GlycosphingolipidIsoglobo(MotifClassifier):
    _class = ("Glycosphingolipid","isoglobo series")
    _motifs = ["GGM.001101"]

class GlycosphingolipidMuco(MotifClassifier):
    _class = ("Glycosphingolipid","muco series")
    _motifs = ["GGM.001102"]

class GlycosphingolipidArthro(MotifClassifier):
    _class = ("Glycosphingolipid","arthro series")
    _motifs = ["GGM.001103"]

class GlycosphingolipidMollu(MotifClassifier):
    _class = ("Glycosphingolipid","mollu series")
    _motifs = ["GGM.001104"]

class GlycosphingolipidGala(MotifClassifier):
    _class = ("Glycosphingolipid","gala series")
    _motifs = ["GGM.001105"]

class GlycosphingolipidLacto(MotifClassifier):
    _class = ("Glycosphingolipid","lacto series")
    _motifs = ["GGM.001106"]

class GlycosphingolipidNeolacto(MotifClassifier):
    _class = ("Glycosphingolipid","neo-lacto series")
    _motifs = ["GGM.001107"]

class GlycosphingolipidGanglio(MotifClassifier):
    _class = ("Glycosphingolipid","ganglio series")
    _motifs = ["GGM.001108"]

class GlycosphingolipidGlobo(MotifClassifier):
    _class = ("Glycosphingolipid","globo series")
    _motifs = ["GGM.001109"]

if __name__ == "__main__":

    from getwiki import GlycanData, Glycan
    w = GlycanData()

    classifier = ClassifierEngine()
    for acc in sys.argv[1:]:
	g = w.get(acc)
	for asn in classifier.assign(g):
	    print g.get('accession'),asn[0],asn[1]
