
from GlycanFormatter import LinearCodeFormat

class GlyLinCodeFilter:
    def __init__(self,glydb):
	self.glydb = glydb
        self.fmt = LinearCodeFormat()
    def __iter__(self):
	return self.next()
    def next(self):
	for gr in self.glydb:
            lincode = True
            try:
                lc = self.fmt.toStr(gr.glycan)
                gr['LinearCode'] = lc
            except KeyError:
                lincode = False
	    gr['lincode'] = lincode
	    yield gr
	
    
