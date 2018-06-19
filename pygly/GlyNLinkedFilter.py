
from subtree import linearcodeSubtreeContainedIn
from GlycanFactory import GlycanFactory

class GlyNLinkedFilter:
    def __init__(self,glydb):
	self.glydb = glydb
	gf = GlycanFactory()
	self.nlc = gf['N-Linked Core']
	self.mnlc = gf['Minimal N-Linked Core']
	self.cmp = linearcodeSubtreeContainedIn()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
	    gr['nlinked'] = self.test(gr.glycan)
	    gr['mini-nlinked'] = self.test1(gr.glycan)
	    yield gr
    def test(self,g):
        return self.cmp.compare(self.nlc,g)
    def test1(self,g):
        return self.cmp.compare(self.mnlc,g)
	
    
