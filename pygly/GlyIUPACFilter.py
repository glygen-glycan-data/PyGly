
from GlycanFormatter import IUPACLinearFormat

class GlyIUPACFilter:
    def __init__(self,glydb):
	self.glydb = glydb
        self.fmt = IUPACLinearFormat()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
            try:
                iupac = self.fmt.toStr(gr.glycan)
                gr['IUPAC'] = iupac
            except KeyError:
                pass
	    yield gr
	
    
