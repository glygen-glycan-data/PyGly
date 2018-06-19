
from subtree import linearcodeSubtreeContainedIn
from GlycanFactory import GlycanFactory
from MonoFactory import MonoFactory
from Monosaccharide import Anomer

class GlyLactosamineFilter:
    def __init__(self,glydb):
	self.glydb = glydb
	gf = GlycanFactory()
        mf = MonoFactory()
        self.GlcNAc = mf['GlcNAc']
        self.GlcNAc.set_anomer(Anomer.missing)
        self.Man = mf['adMan']
        self.Man.set_anomer(Anomer.missing)
	self.lac = gf['Lactosamine N-glycans']
	self.ngc = gf['N-Linked Core']
	self.cmp = linearcodeSubtreeContainedIn()

    def nLinkedAntenae(self, g):
        r = g.root()
        # assert r.compatible(self.GlcNAc)

        # Pick out the GlcNAc - there may be a core-fucosylation
        glcnac2 = filter(lambda m: m.compatible(self.GlcNAc), r.children())[0]
        assert glcnac2.compatible(self.GlcNAc)

        # Only one child
	man1 = glcnac2.first_child()
        assert man1.compatible(self.Man)

        # brach manose in parent pos order (or arbitary if none)...
        man13,man16 = filter(lambda m: m.compatible(self.Man),
		             map(lambda l: l.child(),
                                 sorted(man1.links(),key=lambda l: l.parent_pos())))
        try:
            assert man13.compatible(self.Man)
            assert man16.compatible(self.Man)
        except:
            print man13,"\n",self.Man
            raise
        
        antman13 = map(lambda l: l.child(),
                       sorted(man13.links(),key=lambda l: l.parent_pos()))
        antman16 = map(lambda l: l.child(),
                       sorted(man16.links(),key=lambda l: l.parent_pos()))
        return len(antman13)+len(antman16),antman13,antman16

    def __iter__(self):
	return self.next()

    def next(self):
	for gr in self.glydb:
	    gr['lactosamine'] = self.test(gr.glycan)
	    yield gr
    def test(self,g):
        if self.cmp.compare(self.ngc,g) and \
           self.cmp.compare(g,self.lac):
            # print gr.accession
	    try:
                nant,ant13,ant16 = self.nLinkedAntenae(g)
	    except IndexError:
		return False
            # print nant,ant13,ant16
            if nant == 2:
                if len(ant13) == 1:
                    lac = True
                else:
                    lac = False
            else:
                lac = True
            # print lac
	    return lac
	return False
