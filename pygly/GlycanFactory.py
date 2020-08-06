
from ReferenceTable import ReferenceTable
from GlycanFormatter import GlycoCTFormat
from MonoFactory import MonoFactory
from Monosaccharide import Anomer, Linkage

# Should structures specified by their oxford abbreviation use
# undetermined linkages for their antennae? Currently no, but
# perhaps. Right now we just get one linkage instantiation.

class GlycanFactory(ReferenceTable):
    def __init__(self):
        self.fmt = GlycoCTFormat()
	self.mf = MonoFactory()
        super(GlycanFactory,self).__init__()

    def new(self,key):
        return self[key].clone()

    def parseSection(self,name,kv):
        aliases = [name]
        g = self.fmt.toGlycan('\n'.join(kv['GlycoCT'].split()))
        aliases.extend([ s.strip() for s in kv.get('Aliases','').split(';') ])
        return [(a,g) for a in aliases]
    
    def add_mono(self, parent, name, parent_pos,
                 child_pos=1, anomer=Anomer.beta,
                 parent_type=Linkage.oxygenPreserved,
                 child_type=Linkage.oxygenLost):
        m = self.mf.new(name)
        m.set_anomer(anomer)
        parent.add_child(m,parent_pos=parent_pos,
                         child_pos=child_pos,
                         parent_type=parent_type,
                         child_type=child_type)
        return m

    def oxford2Glycan(self,name):
	if name in self:
	    return self.new(name)
	p = 0
	if name[p] == 'F':
	    g = self.new('FM3')
	    p += 1
	else:
	    g = self.new('M3')
	# print repr(g)
	# print self.fmt.toStr(g)
	# print self
	r = g.root()
	glcnac2 = filter(lambda m: m.compatible(self.mf['GlcNAc']), r.children())[0]
	man1 = glcnac2.children()[0]
	man16 = [l.child() for l in man1.links() if l.parent_pos() == 6][0]
	man13 = [l.child() for l in man1.links() if l.parent_pos() == 3][0]
	assert name[p] == 'A'
	nant = int(name[p+1])
        ant = [None]
	if nant in (1,2,3,4):
	    ant.append(self.add_mono(man13,'GlcNAc',parent_pos=2))
        if nant in (2,3,4):
	    ant.append(self.add_mono(man16,'GlcNAc',parent_pos=2))
	if nant in (3,4):
	    ant.append(self.add_mono(man13,'GlcNAc',parent_pos=4))
        if nant in (4,):
	    ant.append(self.add_mono(man16,'GlcNAc',parent_pos=6))
	p += 2
	if p >= len(name):
	    return g
	if name[p] == 'B':
            b = self.add_mono(man1,'GlcNAc',4)
	    name[p] += 1
	    if p >= len(name):
	        return g
        if name[p] == 'F':
	    nfuc = int(name[p+1])
	    assert (nfuc <= nant)
            for fi in range(1,nfuc+1):
                self.add_mono(ant[fi],'Fuc',parent_pos=6,anomer=Anomer.alpha)
	    p += 2
	    if p >= len(name):
	        return g
	assert(name[p] == 'G')
	ngal = int(name[p+1])
        gal = [None]
	assert (ngal <= nant)
        for gi in range(1,ngal+1):
            gal.append(self.add_mono(ant[gi],'Gal',parent_pos=4))
	p += 2
	if p >= len(name):
	    return g
	assert(name[p] == 'S')
	nsia = int(name[p+1])
        sia = [None]
	assert (nsia <= ngal)
        for si in range(1,nsia+1):
            sia.append(self.add_mono(gal[si],'Neu5Ac',parent_pos=6,child_pos=2,anomer=Anomer.alpha))
        return g
