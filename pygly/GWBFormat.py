
from JavaProgram import GWBFormatter
from GlycanFormatter import GlycoCTFormat

class GWBFormat(object):
    def __init__(self):
        self.fmt = GlycoCTFormat()
        
    def toStr(self,glycan):
	glystr = glycan
	if not isinstance(glystr,basestring):
	    glystr = self.fmt.toStr(glycan)
        writer = GWBFormatter(glystr)
	seq = writer().strip()
	if seq.startswith("Exception"):
	    return None
        return seq
        
        
                          
        
        

    
