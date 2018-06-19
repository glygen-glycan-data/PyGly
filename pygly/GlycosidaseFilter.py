
from Glycosidase import *
from GlycanImage import GlycanImage
import tempfile, os

class CleavedGlycanImageFilter:
    def __init__(self,glydb):
	self.glydb = glydb
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
	    if '-' in gr.accession:
		del gr['image']
		imageWriter = GlycanImage()
		imageWriter.set('scale',1.0)
		imageWriter.set('display','compact')
		fmt = imageWriter.format()
		dummy,tmpoutfile = tempfile.mkstemp(suffix=".%s"%fmt)
		imageWriter.writeImage(gr.glycan,tmpoutfile)
		h=open(tmpoutfile); gr['image'] = h.read(); h.close()
		os.unlink(tmpoutfile)
	    yield gr

class NeuraminidaseFilter:
    def __init__(self,glydb):
	self.glydb = glydb
        self.neu = Neuraminidase()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
            changed = self.neu.cleave(gr.glycan)
	    gr['sialic-acids-cleaved'] = changed
	    if changed > 0:
                gr.accession += '-S%d'%changed
	    yield gr
	
class GalactosidaseFilter:
    def __init__(self,glydb):
	self.glydb = glydb
        self.gal = Galactosidase()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
            changed = self.gal.cleave(gr.glycan)
	    gr['galactose-cleaved'] = changed
	    if changed > 0:
                gr.accession += '-G%d'%changed
	    yield gr
	    
class Beta14GalactosidaseFilter(GalactosidaseFilter):
    def __init__(self,glydb):
	self.glydb = glydb
        self.gal = Beta14Galactosidase()

class FucosidaseFilter:
    def __init__(self,glydb):
	self.glydb = glydb
        self.gal = Fucosidase()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
            changed = self.gal.cleave(gr.glycan)
	    gr['fucose-cleaved'] = changed
	    if changed > 0:
                gr.accession += '-F%d'%changed
	    yield gr
	    
