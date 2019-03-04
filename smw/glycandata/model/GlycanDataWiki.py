
__all__ = [ "GlycanDataWiki", "Glycan", "Annotation" ]

import sys
from operator import itemgetter

from smw import SMW

class Glycan(SMW.SMWClass):
    template = 'Glycan'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('accession')
        return kwargs.get('accession')
    
    def toPython(self,data):
	data = super(Glycan,self).toPython(data)

        # mw is a float
        if isinstance(data.get('mw'),basestring):
            data['mw'] = float(data.get('mw'))

        # monocount is an integer
        if isinstance(data.get('monocount'),basestring):
            data['mw'] = int(data.get('monocount'))

	return data

    def toTemplate(self,data):
	data = super(Glycan,self).toTemplate(data)

        if 'mw' in data:
            data['mw'] = str(data.get('mw'))
        
        if 'monocount' in data:
            data['monocount'] = str(data.get('monocount'))
        
	return data

class Annotation(SMW.SMWSite):
    template = 'Annotation'

class GlycanDataWiki(SMW.SMWSite):
    _name = 'glycandata'

    def __init__(self):
        super(GlycanDataWiki,self).__init__()

    # Assumes accession == pagename, see above
    def get(self,accession):
        return super(GlycanDataWiki,self).get(accession)

    def iterglycan(self):
	for pagename in self.itercat('Glycan'):
	    m  = self.get(pagename)
            yield m
