
from Glycan import Glycan
from MonoFactory import MonoFactory
import re

mf = MonoFactory()

class Glycosidase:
    pass

class LeafCleaver(Glycosidase):
    def cleave(self,gly):
	ncleaved = 0
	newcleaved = True
	while newcleaved:
	  newcleaved = False
          for l in gly.all_links():
            m = l.child()
            if m.compatible(self.monosaccharide) and not m.has_children():
                l.parent().del_link(l)
		ncleaved += 1
		newcleaved = True
	return ncleaved

class Serial(Glycosidase):
    def __init__(self,*glycosidases):
	self._glycosidases = glycosidases

    def number(self):
	return len(self._glycosidases)

    def cleave(self,gly):
        ncleaved = []
	for gcd in self._glycosidases:
	    ncleaved.append(gdc.cleave(gly))
 	return tuple(ncleaved)

class Simultaneous(Glycosidase):
    def __init__(self,*glycosidases):
	self._serial = Serial(*glycosidases)

    def number(self):
	return self._serial.number()

    def cleave(self,gly):
        ncleaved = [0]*self.number()
	while True:
	    ncl = self._serial.cleave(gly)
	    if sum(ncl) == 0:
		break
	    ncleaved = map(sum,zip(ncleaved,ncl))
 	return tuple(ncleaved)
                
class Neuraminidase(LeafCleaver):
    monosaccharide = mf.new('NeuAc')

class Sialidase(Neuraminidase):
    pass
                
class Galactosidase(LeafCleaver):
    monosaccharide = mf.new('Gal')

class Beta14Galactosidase(Galactosidase):
    # this one is special, because Fucose on the the sub-terminal GlcNAc
    # protects the Gal.
    parent = mf.new('GlcNAc')
    protector = mf.new('Fuc')
    def cleave(self,gly):
	ncleaved = 0
        for l in gly.all_links():
            m = l.child()
	    if not m.compatible(self.monosaccharide):
		continue
	    if m.has_children():
		continue
	    p = l.parent()
	    if not p.compatible(self.parent):
		continue
	    anyfuc = False
	    for ch in p.children():
		if ch.compatible(self.protector):
		    anyfuc = True
		    break
	    if anyfuc:
		continue
            l.parent().del_link(l)
	    ncleaved += 1
	return ncleaved

class Fucosidase(LeafCleaver):
    monosaccharide = mf.new('Fuc')
