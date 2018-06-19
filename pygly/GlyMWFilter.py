
class GlyMWFilter:
    def __init__(self,glydb,comp_table,perm_comp_table,mass_table,**kw):
	self.glydb = glydb
	self.cmp = comp_table
        self.perm = perm_comp_table
	self.mt = mass_table
	self.addh2o = 0
        self.addperm = 0
	if kw.get('native_ends',False):
	    self.addh2o = 2*self.mt['H'] + self.mt['O']
	if kw.get('permethylated_ends',False):
            self.addperm = self.mt['C'] + 3*self.mt['H'] + self.mt['O'] + self.mt['C'] + 3*self.mt['H']
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
	    mw = 0.0; pmw = 0.0
	    bad = False; pbad = False;
	    for m in gr.glycan.all_nodes():
		try:
		    ec = m.composition(self.cmp)
		    mw += ec.mass(self.mt)
		except KeyError:
		    bad = True
                try:
		    ec = m.composition(self.perm)
		    pmw += ec.mass(self.mt)
		except KeyError:
		    pbad = True
	    if not bad:
	        gr['molecular_weight'] = mw + self.addh2o
	    if not pbad:
	        gr['permethylated_molecular_weight'] = pmw + self.addperm
	    yield gr
