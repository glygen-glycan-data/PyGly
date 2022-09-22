
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
	seq = writer().splitlines()[-1].strip()
	if "Exception" in seq or "org.glycoinfo" in seq or "GlycoCT2GWB.main" in seq:
	    return None
        return seq
        
        
                          
        
        

    
