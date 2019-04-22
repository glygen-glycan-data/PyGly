
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

        if '_subobjs' in data:
	    data['annotations'] = sorted(data['_subobjs'],key=lambda an: (an.get('type'),an.get('property'),an.get('source')))
            del data['_subobjs']

	return data

    def toTemplate(self,data):
	data = super(Glycan,self).toTemplate(data)

        if 'annotations' in data:
            data['_subobjs'] = sorted(data['annotations'],key=lambda an: (an.get('type'),an.get('property'),an.get('source')))
            del data['annotations']
        
	return data

    def add_annotation(self,**kwargs):
        assert kwargs.get('type')
        assert kwargs.get('property')
        assert kwargs.get('source')
        self.append('annotations',Annotation(**kwargs))

    def set_annotation(self,**kwargs):
        assert kwargs.get('type')
        assert kwargs.get('property')
        assert kwargs.get('source')        
        ans = list(self.annotations(type=kwargs.get('type'),
                                    property=kwargs.get('property'),
                                    source=kwargs.get('source')))
        goodvalue = (kwargs.get('value') not in (None,"",[]))
        if len(ans) == 0 and goodvalue:
            self.add_annotation(**kwargs)
        elif len(ans) == 1 and goodvalue:
            ans[0].update(**kwargs)
        else:
	    self.delete_annotations(type=kwargs.get('type'),property=kwargs.get('property'),source=kwargs.get('source'))
            if goodvalue:
                self.add_annotation(**kwargs)

    def delete_annotations(self,**kwargs):
        todel = list(self.annotations(**kwargs))
        ans = self.get('annotations',[])
        for i in range(len(ans)-1,-1,-1):
            if ans[i] in todel:
                del ans[i]
        self.set('annotations',ans)

    def annotations(self,type=None,property=None,source=None):
        for an in self.get('annotations',[]):
            if (type == None or an.get('type') == type) and \
               (property == None or an.get('property') == property) and \
               (source == None or an.get('source') == source):
                yield an

class Annotation(SMW.SMWClass):
    template = 'Annotation'

    @staticmethod
    def intstrvalue(v):
        try:
            return int(v),""
        except:
           pass
        try:
           return float(v),""
        except:
           pass
        return 1e+20,v

    def toPython(self,data):
        data = super(Annotation,self).toPython(data)

        # value may be a list
        if isinstance(data.get('value'),basestring):
            data['value'] = sorted(map(lambda s: s.strip(),data.get('value').split(';')),key=self.intstrvalue)
        elif isinstance(data.get('value'),float) or isinstance(data.get('value'),int):
            data['value'] = [ str(data['value']) ]
        
        return data

    def toTemplate(self,data):
        data = super(Annotation,self).toTemplate(data)

        if 'value' in data:
            data['value'] = ";".join(map(str,sorted(data['value'],key=self.intstrvalue)))

        return data

class GlycanDataWiki(SMW.SMWSite):
    _name = 'glycandata'

    def __init__(self):
        super(GlycanDataWiki,self).__init__()

    # Assumes accession == pagename, see above
    def get(self,accession):
        return super(GlycanDataWiki,self).get(accession)

    def iterglycan(self):
        for pagename in sorted(self.itercat('Glycan')):
            m  = self.get(pagename)
            yield m
