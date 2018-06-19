
from JavaProgram import GlycoCT2Image
from GlycanFormatter import GlycoCTFormat

class GlycanImage(object):
    def __init__(self):
        self._scale = 1.0
        self._orientation = "RL"
        self._display = "normalinfo"
        self._notation = "cfg"
        self._redend = True
        self._format = "png"
        self.fmt = GlycoCTFormat()
        
    def scale(self,value=None):
        if value == None:
            return self._scale
        self._scale=float(value)

    def reducing_end(self,value=None):
        if value == None:
            return self._redend
        self._redend=bool(value)

    def notation(self,value=None):
        if value == None:
            return self._notation
        self._notation=value

    def format(self,value=None):
        if value == None:
            return self._format
        self._format=value

    def orientation(self,value=None):
        if value == None:
            return self._orientation
        self._orientation=value

    def display(self,value=None):
        if value == None:
            return self._display
        self._display=value

    def set(self,key,value):
	if not hasattr(self,key):
	    raise KeyError()
	getattr(self,key)(value)	

    def writeImage(self,glycan,filename):
	glystr = glycan
	if not isinstance(glystr,basestring):
	    glystr = self.fmt.toStr(glycan)
        imageWriter = GlycoCT2Image(glystr,
                                    filename,
                                    format=self._format,
                                    scale=self._scale,
                                    redend=self._redend,
                                    orient=self._orientation,
                                    display=self._display,
                                    notation=self._notation)
        return imageWriter()
        
        
                          
        
        

    
