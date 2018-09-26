
class GlyRecord(dict):
    def __init__(self,source,accession,glycan,**kw):
        self.source    = source
	self.accession = accession
	self.glycan    = glycan
	self.update(kw)
    def __str__(self):
	s = ""
	s += "source = %s\n"%(self.source,)
	s += "accession = %s\n"%(self.accession,)
	for k,v in self.iteritems():
	    if k in ('image',):
		continue
	    s += "%s = %s\n"%(k,v)
	try:
	    s += "glycan = \n%s\n"%(self.glycan,)
	except KeyError:
	    pass
	return s
