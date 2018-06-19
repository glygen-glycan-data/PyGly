
from subtree import linearcodeSubtreeContainedIn
from GlycanFactory import GlycanFactory

class GlyHighMannoseFilter:
    def __init__(self,glydb):
	self.glydb = glydb
	gf = GlycanFactory()
	self.ngc = gf['N-Linked Core']
	self.cmp = linearcodeSubtreeContainedIn()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
	    gr['high-mannose'] = self.test(gr.glycan)
	    yield gr
    def test(self,g):
        if not self.cmp.compare(self.ngc,g):
	    return False
	try:
            symcomp,eltcomp = g.composition()
	except KeyError:
	    return False
	# this is implied by subtree comparison with N-glycan core.
	# print symcomp
	assert(symcomp['GlcNAc'] >= 2)
	assert(symcomp['Man'] >= 3)
	# Assume we cannot be fucosylated
	return ((symcomp['GlcNAc'] == 2) and (sum(symcomp.values()) == symcomp['Man'] + 2))
	    
