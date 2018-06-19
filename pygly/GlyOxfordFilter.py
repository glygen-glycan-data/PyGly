
from MonoFactory import MonoFactory
from Monosaccharide import Anomer

class GlyOxfordFilter:
    def __init__(self,glydb):
	self.glydb = glydb
        mf = MonoFactory()
	self.GlcNAc = mf['GlcNAc']
        self.GlcNAc.set_anomer(Anomer.missing)
        self.Man = mf['adMan']
        self.Man.set_anomer(Anomer.missing)
	self.Fuc = mf['Fuc']
	self.Fuc.set_anomer(Anomer.missing)

    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
	    if gr.get('nlinked',False):
		try:
	            symcomp,eltcomp = gr.glycan.composition()	
		except KeyError:
		    continue
	        if gr.get('high-mannose',False):
		    gr['oxford'] = 'M%d'%symcomp['Man']
		elif gr.get('lactosamine',False):
		    nant = symcomp['GlcNAc'] - 2
		    bisec = False
		    ngal = symcomp['Gal']
		    nsial = symcomp['NeuAc']
		    nfuc = symcomp['Fuc']
		    corefuc = False

		    # check for core fuc
		    # check for bisecting GlcNAc
		    r = gr.glycan.root()
		    # assert r.compatible(self.GlcNAc)

		    corefuc = (len(filter(lambda m: m.compatible(self.Fuc), r.children())) > 0)
		    if corefuc:
		        nfuc -= 1

		    # Pick out the GlcNAc - there may be a core-fucosylation
		    glcnac2 = filter(lambda m: m.compatible(self.GlcNAc), r.children())[0]
                    assert glcnac2.compatible(self.GlcNAc)

		    # Only one child
		    man1 = glcnac2.first_child()
		    assert man1.compatible(self.Man)

		    bisec = (len(filter(lambda m: m.compatible(self.GlcNAc), man1.children())) > 0)
		    if bisec:
		        nant -= 1
		    ox = ""
		    if corefuc:
		        ox += "F"
		    ox += "A%d"%nant
		    if bisec:
		        ox += "B"
		    if nfuc > 0:
		        ox += "F%d"%nfuc
		    if ngal > 0:
	                ox += "G%d"%ngal
		    if nsial > 0:
		        ox += "S%d"%nsial
		    if nant > 0:
                        gr['oxford'] = ox
	    yield gr
	            
