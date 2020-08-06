
from CompositionTable import Composition
from MonoFormatter import MassSym
from Monosaccharide import Monosaccharide, Linkage, Anomer, Substituent, Mod

from collections import defaultdict

tochar = defaultdict(lambda: '?',filter(None,map(lambda s: s.strip(),"""
Gal H
Hex H
GlcNAc N
HexNAc N
NeuAc S
Fuc F
Man H
Glc H
GalNAc N
Xyl X
NeuGc J
""".splitlines())))

def comp2str1(comp):
    comp1 = defaultdict(int)
    for k,v in comp.items():
        comp1[tochar[k]] += v
    s = ""
    for sym in "HNFSXJ?":
        c = comp1.get(sym,0)
        if c == 1:
            s += sym
        elif c > 1:
            s += "%s%d"%(sym,c)
    return s

class GlyCompFilter:
    def __init__(self,glydb):
	self.glydb = glydb
	self.syms = MassSym()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
	    comp = Composition()
	    for m in gr.glycan.all_nodes():
		try:
		    sym = self.syms.toStr(m)
		except KeyError:
                    if m.noring():
                        if m.count_mod() == 1 and \
                               m.count_mod(Mod.aldi) == 1:
                            aglycon = 'aldi?'
                            m1 = m.clone()
                            m1.clear_mods()                    
                            try:
                                sym=self.syms.toStr(m1)
                            except KeyError:
                                sym = 'Xxx'
                        else:
                            sym = 'Xxx'
                    else:
                        sym = 'Xxx'
		comp[sym] += 1
	    gr['composition'] = comp
	    gr['compstr'] = comp2str1(comp)
	    yield gr
