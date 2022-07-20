
from . JavaProgram import GlycoCT2Image
from . GlycanFormatter import GlycoCTFormat

try:
  basestring
except NameError:
  basestring = str

class GlycanImage(object):
    def __init__(self):
        self._scale = 1.0
        self._orientation = "RL"
        self._display = "normalinfo"
        self._notation = "snfg"
        self._redend = True
        self._format = "png"
        self._opaque = True
        self._force = False
        self._verbose = False
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

    def opaque(self, value=None):
        if value == None:
            return self._opaque
        self._opaque=value

    def force(self, value=None):
        if value == None:
            return self._force
        self._force=value

    def verbose(self, value=None):
        if value == None:
           return self._verbose
        self._verbose = value

    def set(self,key,value):
        if not hasattr(self,key):
            raise KeyError(key)
        getattr(self,key)(value)

    def get(self,key):
        if not hasattr(self,key):
            raise KeyError(key)
        return getattr(self,key)()

    def writeImage(self,glycan,filename):
        glystr = glycan
        if not isinstance(glystr,basestring):
            glystr = self.fmt.toStr(glycan)
        imageWriter = GlycoCT2Image(glystr,
                                    filename,
                                    format=self._format,
                                    force=str(self._force).lower(),
                                    scale=self._scale,
                                    redend=str(self._redend).lower(),
                                    orient=self._orientation,
                                    display=self._display,
                                    notation=self._notation,
                                    opaque=str(self._opaque).lower(),
                                    verbose=self._verbose)
        return imageWriter()
