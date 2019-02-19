
__all__ = [ "GlycoMotifWiki", "Collection", "Motif",
            "GlyTouCanMotif", "AllMotif", "CCRCMotif", "GlycoEpitopeMotif",
            "GlydinMotif",
            "GlydinCummingsMotif", "GlydinHayesMotif", "GlydinCermavMotif",
            "GlydinSugarbindMotif", "GlydinBioligoMotif",
            "UniCarbMotif"]

import sys

from smw import SMW

class Collection(SMW.SMWClass):
    template = 'Collection'

    def __init__(self,**kwargs):
	if kwargs.get('primary') == None:
	    kwargs['primary'] = True
        super(Collection,self).__init__(**kwargs)

    def toPython(self,data):
	data = super(Collection,self).toPython(data)

	if isinstance(data.get('primary'),basestring):
	    data['primary'] = self.asboolean(data.get('primary'))
	
	return data

    def toTemplate(self,data):
	data = super(Collection,self).toTemplate(data)

	if 'primary' in data:
	    data['primary'] = ("true" if data['primary'] else "false")
	
	return data
	

class Motif(SMW.SMWClass):
    template = 'Motif'
    aglyconvalues = ['Cer','Ser/Thr','Asn','R','Other']

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

        # aglycon is comma separated, with specific possible values, sorted, so behaves as set
        if isinstance(data.get('aglycon'),basestring):
            data['aglycon'] = sorted(map(self.asaglycon,map(lambda s: s.strip(),data.get('aglycon').split(','))))
        elif data.get('aglycon') != None:
            data['aglycon'] = sorted(map(self.asaglycon,data.get('aglycon')))

        # sameas is comma separated
        if isinstance(data.get('sameas'),basestring):
            data['sameas'] = sorted(map(lambda s: s.strip(),data.get('sameas').split(',')))

        # topology is comma separated
        if isinstance(data.get('topology'), basestring):
            data['topology'] = map(lambda s: s.strip(), data.get('topology').split(','))
        
        if isinstance(data.get("displayhgv"), basestring):
            data["displayhgv"] = self.asboolean(data.get("displayhgv"))

        # redend is a list of booleans, sorted, so behaves as set
        if isinstance(data.get('redend'),basestring):
            data['redend'] = sorted(map(self.asboolean,map(lambda s: s.strip(),data.get('redend').split(','))))
        elif data.get('redend') in (True,False):
            data['redend'] = [ data.get('redend') ]
        elif data.get('redend') != None:
            data['redend'] = sorted(map(self.asboolean,data.get('redend')))

        # Strip <pre> and </pre> if it is there
        # if 'wurcs' in data:
            #     data['wurcs'] = data.get('wurcs').lstrip('<pre>').rstrip('</pre>')

        # Strip <pre> and </pre> if it is there
        # if 'glycoct' in data:
            #     data['glycoct'] = data.get('glycoct').lstrip('<pre>').rstrip('</pre>')

        return data

    def toTemplate(self,data):
        data = super(Motif,self).toTemplate(data)

        if 'name' in data:
            data['name'] = "\n".join(data['name'])

        if 'sameas' in data:
            data['sameas'] = ",".join(sorted(data['sameas']))

        if 'aglycon' in data:
            data['aglycon'] = ",".join(sorted(data['aglycon']))

        if 'redend' in data:
            data['redend'] = ",".join(sorted(("true" if redend else "false") for redend in data['redend']))

        if 'topology' in data:
            data['topology'] = ",".join(data['topology'])

        if "displayhgv" in data:
            data["displayhgv"] = ("true" if data["displayhgv"] else "false")

        # if 'wurcs' in data:
        #     data['wurcs'] = "<pre>" + data['wurcs'] + "</pre>"

        # if 'glycoct' in data:
        #     data['glycoct'] = "<pre>" + data['glycoct'] + "</pre>"

        return data


import findpygly
from pygly.GlyTouCan import GlyTouCan

class GlyTouCanMotif(Motif):
    gtc = None
    id = 'GTC'
    def __init__(self,**kwargs):
        assert kwargs.get('accession') != None
        kwargs['collection'] = self.id
        if kwargs.get('glytoucan') == None:
            kwargs['glytoucan'] = kwargs['accession']
        if kwargs.get('wurcs') == None or kwargs.get('glycoct') == None:
            if not self.gtc:
                self.gtc = GlyTouCan()
            if kwargs.get('wurcs') == None:
                kwargs['wurcs'] = self.gtc.getseq(kwargs['glytoucan'],'wurcs')
            if kwargs.get('glycoct') == None:
                kwargs['glycoct'] = self.gtc.getseq(kwargs['glytoucan'],'glycoct')
        super(GlyTouCanMotif,self).__init__(**kwargs)

class AllMotif(GlyTouCanMotif):
    id = 'GM'

class CCRCMotif(GlyTouCanMotif):
    id = 'CCRC'
    def __init__(self,**kwargs):
        assert kwargs.get('accession') != None
        assert kwargs.get('glytoucan') != None
        super(CCRCMotif,self).__init__(**kwargs)

class GlycoEpitopeMotif(CCRCMotif):
    id = 'GE'

class GlydinMotif(CCRCMotif):
    id = 'GD'

class GlydinCummingsMotif(CCRCMotif):
    id = 'GDC'

class GlydinHayesMotif(CCRCMotif):
    id = 'GDH'

class GlydinCermavMotif(CCRCMotif):
    id = 'GDV'

class GlydinSugarbindMotif(CCRCMotif):
    id = 'GDSB'

class GlydinBioligoMotif(CCRCMotif):
    id = 'GDB'

class UniCarbMotif(CCRCMotif):
    id = 'UCM'

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
