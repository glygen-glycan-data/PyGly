
from SMW import SMWSite, SMWClass
from GlyTouCan import GlyTouCan

class Collection(SMWClass):
    template = 'Collection'
    directsubmit = {}

    def normalize(self):
        for k in self.keys():
            if self.get(k) in (None,""):
                self.delete(k)

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('pagename')
        return kwargs.get('pagename')

    def astemplate(self):
        l = [ "{{%s"%(self.template,) ] 
        for k in sorted(self.keys()):
            v = self.get(k)
	    if k in ('pagename','template','class'):
		continue
            if v in (None,""):
                continue
            l.append("|%s=%s"%(k,v))
        l.append("}}")
        return "\n".join(l)

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
	if self.has('redend') and isinstance(self.get('redend'),basestring):
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
    collection = 'GlyTouCan'
    def __init__(self,**kwargs):
        assert 'accession' in kwargs
        kwargs['collection'] = self.collection
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

class GlyGenMotif(GlyTouCanMotif):
    collection = 'GlyGen'

class CCRCMotif(GlyTouCanMotif):
    collection = 'UGA-CCRC'
    def __init__(self,**kwargs):
        assert 'accession' in kwargs
        assert 'glytoucan' in kwargs
        super(CCRCMotif,self).__init__(**kwargs)

class GlycoEpitopeMotif(CCRCMotif):
    collection = 'GlycoEpitope'

class GlycoMotifWiki(SMWSite):
    # Suitable for running from dalton
    prefix = 'glycomotif'
    host = 'codon'
    port = 8282

    template2class = {'Motif': Motif, 
		      'Collection': Collection}
    dump_exclude_categories = set(['Motif'])

    def get(self,pagename=None,collection=None,accession=None):
	if pagename:
	    return super(GlycoMotifWiki,self).get(pagename)
        return super(GlycoMotifWiki,self).get(Motif.pagename(collection=collection,accession=accession))

    def itermotif(self):
	for pagename in self.itercat('Motif'):
	    yield self.get(pagename)

    def itercollection(self):
	for pagename in self.itercat('Collection'):
	    yield self.get(pagename)

if __name__ == "__main__":
    import sys, os

    smw = GlycoMotifWiki()
    print smw.get(collection='GlyTouCan',accession='G00026MO')
    print smw.get(collection='GlyTouCan',accession='G00029MO')

    motif = GlyTouCanMotif(accession="G00028MO",
                           name="N-Glycan high mannose")
    smw.put(motif)
    print smw.get(collection='GlyTouCan',accession='G00028MO')
