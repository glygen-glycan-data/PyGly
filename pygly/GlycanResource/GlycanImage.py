
try:
    from pygly.GlycanImage import GlycanImage as GI
except:
    # from GlycanImage import GlycanImage as GI
    pass

from GlyTouCan import GlyTouCan

class GlycanImage(object):
    def __init__(self):
	self.imageWriter = GI()
	self.gtc = GlyTouCan(prefetch=False)
    def write(self,*args,**kw):
	for k,v in kw.items():
	    self.imageWriter.set(k,v)
	self.imageWriter.set('format','png')
	for acc in args:
	    glycoct = self.gtc.glycoct(acc)
	    if not glycoct:
	        continue
            self.imageWriter.writeImage(glycoct,acc+".png")
