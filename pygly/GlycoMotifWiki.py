
from SMW import SMWSite, SMWClass
from GlyTouCan import GlyTouCan

class MotifCollection(SMWClass):
    template = 'MotifCollection'

    def normalize(self):
        for k in self.keys():
            if self.get(k) in (None,""):
                self.delete(k)

class Motif(SMWClass):
    template = 'Motif'
    directsubmit = {'wurcs': "WURCS", 'glycoct': "GlycoCT"}

    def normalize(self):
        for k in self.keys():
            if self.get(k) in (None,""):
                self.delete(k)
        if self.has('name') and isinstance(self.get('name'),basestring):
            self.set('name',self.get('name').split('\n'))
        if self.has('sameas') and isinstance(self.get('sameas'),basestring):
            self.set('sameas',self.get('sameas').split(','))
	if self.has('redend'):
	    self.set('redend',self.get('redend') == "True")
        if self.has('wurcs'):
            self.set('wurcs',self.get('wurcs').lstrip('<pre>').rstrip('</pre>'))
        if self.has('glycoct'):
            self.set('glycoct',self.get('glycoct').lstrip('<pre>').rstrip('</pre>'))

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('collection') and kwargs.get('accession')
        return "%s.%s"%(kwargs.get('collection'),kwargs.get('accession'))

    def astemplate(self):
        l = [ "{{%s"%(self.template,) ] 
        for k in sorted(self.keys()):
            v = self.get(k)
	    if k in ('pagename','template','class'):
		continue
            if v in (None,""):
                continue
            if k == 'name':
                v = "\n".join(v)
	    elif k == 'sameas':
		v = ",".join(v)
            elif k == 'wurcs':
                v = "<pre>%s</pre>"%(v,)
            elif k == 'glycoct':
                v = "<pre>%s</pre>"%(v,)
	    elif k == 'redend':
		v = str(v)
            l.append("|%s=%s"%(k,v))
        l.append("}}")
        return "\n".join(l)

class GlyTouCanMotif(Motif):
    gtc = None
    def __init__(self,**kwargs):
        assert 'accession' in kwargs
        kwargs['collection'] = 'GlyTouCan'
        if 'glytoucan' not in kwargs:
            kwargs['glytoucan'] = kwargs['accession']
        if 'wurcs' not in kwargs or 'glycoct' not in kwargs:
            if not self.gtc:
                self.gtc = GlyTouCan()
            if 'wurcs' not in kwargs:
                kwargs['wurcs'] = self.gtc.getseq(kwargs['glytoucan'],'wurcs')
            if 'glycoct' not in kwargs:
                kwargs['glycoct'] = self.gtc.getseq(kwargs['glytoucan'],'glycoct')
        super(GlyTouCanMotif,self).__init__(**kwargs)

class GlyGenMotif(Motif):
    gtc = None
    def __init__(self,**kwargs):
        assert 'accession' in kwargs
	# This should be the pagename
        kwargs['collection'] = 'GlyGen'
        if 'glytoucan' not in kwargs:
            kwargs['glytoucan'] = kwargs['accession']
        if 'wurcs' not in kwargs or 'glycoct' not in kwargs:
            if not self.gtc:
                self.gtc = GlyTouCan()
            if 'wurcs' not in kwargs:
                kwargs['wurcs'] = self.gtc.getseq(kwargs['glytoucan'],'wurcs')
            if 'glycoct' not in kwargs:
                kwargs['glycoct'] = self.gtc.getseq(kwargs['glytoucan'],'glycoct')
        super(GlyGenMotif,self).__init__(**kwargs)

class CCRCMotif(Motif):
    gtc = None
    def __init__(self,**kwargs):
        assert 'accession' in kwargs
        assert 'glytoucan' in kwargs
        kwargs['collection'] = 'CCRC'
        if 'wurcs' not in kwargs or 'glycoct' not in kwargs:
            if not self.gtc:
                self.gtc = GlyTouCan()
            if 'wurcs' not in kwargs:
                kwargs['wurcs'] = self.gtc.getseq(kwargs['glytoucan'],'wurcs')
            if 'glycoct' not in kwargs:
                kwargs['glycoct'] = self.gtc.getseq(kwargs['glytoucan'],'glycoct')
        super(CCRCMotif,self).__init__(**kwargs)

class GlycoMotifWiki(SMWSite):
    # Suitable for running from dalton
    prefix = 'glycomotif'
    host = 'codon'
    port = 8282

    template2class = {'Motif': Motif, 
		      'MotifCollection': MotifCollection}

    def get(self,pagename=None,collection=None,accession=None):
	if pagename:
	    return super(GlycoMotifWiki,self).get(pagename)
        return super(GlycoMotifWiki,self).get(Motif.pagename(collection=collection,accession=accession))

    def itermotif(self):
	for pagename in self.itercat('Motif'):
	    yield self.get(pagename)

    def itercollection(self):
	for pagename in self.itercat('MotifCollection'):
	    yield self.get(pagename)

if __name__ == "__main__":
    import sys, os

    smw = GlycoMotifWiki(username='edwardsnj',password='XXXXXXXXXXXXXX_FIX_XXXXXXXXXXXXXX')
    print smw.get('GlyTouCan','G00026MO')
    print smw.get('GlyTouCan','G00029MO')

    motif = GlyTouCanMotif(accession="G00028MO",
                           name="N-Glycan high mannose")
    smw.put(motif)
    print smw.get('GlyTouCan','G00028MO')