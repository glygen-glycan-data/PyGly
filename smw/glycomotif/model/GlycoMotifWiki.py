
__all__ = [ "GlycoMotifWiki", "Collection", "Motif",
            "GlyTouCanMotif", "AllMotif", "CCRCMotif", "GlycoEpitopeMotif" ]

from smw import SMW

class Collection(SMW.SMWClass):
    template = 'Collection'

class Motif(SMW.SMWClass):
    template = 'Motif'
    directsubmit = {'wurcs': "WURCS", 'glycoct': "GlycoCT"}
    aglyconvalues = ['Cer','Ser/Thr','Asn','R']

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('collection') and kwargs.get('accession')
        return "%s.%s"%(kwargs.get('collection'),kwargs.get('accession'))

    def asaglycon(self,value):
	if value not in Motif.aglyconvalues:
	    raise ValueError("Bad aglycon value: %r"%(value,))
	return value

    def toPython(self,data):
	data = super(Motif,self).toPython(data)
	
	# name is newline separated
        if isinstance(data.get('name'),basestring):
            data['name'] = map(lambda s: s.strip(),data.get('name').split('\n'))

	# aglycon is comma separated, with specific possible values
        if isinstance(data.get('aglycon'),basestring):
            data['aglycon'] = sorted(map(self.asaglycon,map(lambda s: s.strip(),data.get('aglycon').split(','))))
	elif data.get('aglycon') != None:
	    data['aglycon'] = sorted(map(self.asaglycon,data.get('aglycon')))

	# sameas is comma separated
        if isinstance(data.get('sameas'),basestring):
            data['sameas'] = map(lambda s: s.strip(),data.get('sameas').split(','))

	# redend is a boolean
	if isinstance(data.get('redend'),basestring):
	    data['redend'] = self.asboolean(data.get('redend'))

	# Strip <pre> and </pre> if it is there 
	if 'wurcs' in data:
            data['wurcs'] = data.get('wurcs').lstrip('<pre>').rstrip('</pre>')

	# Strip <pre> and </pre> if it is there 
	if 'glycoct' in data:
            data['glycoct'] = data.get('glycoct').lstrip('<pre>').rstrip('</pre>')

	return data

    def toTemplate(self,data):
	data = super(Motif,self).toTemplate(data)
	
	if 'name' in data:
	    data['name'] = "\n".join(data['name'])

	if 'sameas' in data:
	    data['sameas'] = ",".join(data['sameas'])

	if 'aglycon' in data:
	    data['aglycon'] = ",".join(sorted(data['aglycon']))

	if 'redend' in data:
	    data['redend'] = ("true" if data['redend'] else "false")

	if 'wurcs' in data:
	    data['wurcs'] = "<pre>" + data['wurcs'] + "</pre>"
	    
	if 'glycoct' in data:
	    data['glycoct'] = "<pre>" + data['glycoct'] + "</pre>"
	    
	return data


import findpygly
from pygly.GlyTouCan import GlyTouCan

class GlyTouCanMotif(Motif):
    gtc = None
    id = 'GTC'
    def __init__(self,**kwargs):
        assert 'accession' in kwargs
        kwargs['collection'] = self.id
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

class AllMotif(GlyTouCanMotif):
    id = 'GM'

class CCRCMotif(GlyTouCanMotif):
    id = 'CCRC'
    def __init__(self,**kwargs):
        assert 'accession' in kwargs
        assert 'glytoucan' in kwargs
        super(CCRCMotif,self).__init__(**kwargs)

class GlycoEpitopeMotif(CCRCMotif):
    id = 'GE'

class GlycoMotifWiki(SMW.SMWSite):
    _name = 'glycomotif'

    def get(self,pagename=None,collection=None,accession=None):
	if pagename:
	    return super(GlycoMotifWiki,self).get(pagename)
        return super(GlycoMotifWiki,self).get(Motif.pagename(collection=collection,accession=accession))

    def itermotif(self,collection=None):
        if collection != None and not isinstance(collection,basestring):
            collection = collection.id
	for pagename in self.itercat('Motif'):
	    m  = self.get(pagename)
            if not collection or m.get('collection') == collection:
                yield m

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
