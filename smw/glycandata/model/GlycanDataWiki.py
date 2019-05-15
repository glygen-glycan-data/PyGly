
__all__ = [ "GlycanDataWiki", "Glycan", "Annotation" ]

import sys
from operator import itemgetter

from smw import SMW

import findpygly
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format, GlycanParseError

class Glycan(SMW.SMWClass):
    template = 'Glycan'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('accession')
        return kwargs.get('accession')
    
    def toPython(self,data):
	data = super(Glycan,self).toPython(data)

	if '_subobjs' in data:
	    del data['_subobjs']

	return data

    def toTemplate(self,data):
	data = super(Glycan,self).toTemplate(data)

	return data

    def set_annotation(self,**kwargs):
        if 'annotation' in kwargs:
            ann = kwargs.get('annotation')
        else:
            ann = Annotation(**kwargs)
        if not ann.goodvalue():
            return 
        if not hasattr(self,'_annotations'):
            self._annotations = dict()
        self._annotations[ann.key()] = ann

    def delete_annotations(self,**kwargs):
        if not hasattr(self,'_annotations'):
            return
        for an in list(self.annotations(**kwargs)):
            del self._annotations[an.key()]

    def count_annotations(self,**kwargs):
	return len(list(self.annotations(**kwargs)))

    def get_annotation_value(self,property=None,**kwargs):
	return str(self.get_annotation(property=property,**kwargs).get('value'))

    def get_annotation(self,property=None,**kwargs):
	if property != None:
	    kwargs['property'] = property
	anns = list(self.annotations(**kwargs))
	if len(anns) == 0:
	    raise LookupError("No matching annotations")
	elif len(anns) > 1:
	    raise LookupError("Too many annotations")
	return anns[0]

    def has_annotations(self,**kwargs):
	for an in self.annotations(**kwargs):
	    return True
        return False

    def annotations(self,type=None,property=None,source=None):
        if hasattr(self,'_annotations'):
            for key,an in sorted(self._annotations.items()):
                if (type == None or an.get('type') == type) and \
                       (property == None or an.get('property') == property) and \
                       (source == None or an.get('source') == source):
                    yield an

    def __str__(self):
        sl = [ super(Glycan,self).__str__() ]
        for an in self.annotations():
           sl.append(str(an))
        return "\n".join(sl)

    glycoct_format = GlycoCTFormat()
    wurcs_format = WURCS20Format()
    def getGlycan(self):
        try:
            sequence = self.get_annotation_value('GlycoCT')
            return self.glycoct_format.toGlycan(sequence)
        except (LookupError,GlycanParseError):
            pass
        try:
            sequence = self.get_annotation_value('WURCS')
            return self.wurcs_format.toGlycan(sequence)
        except (LookupError,GlycanParseError):
            pass
        return None

class Annotation(SMW.SMWClass):
    template = 'Annotation'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('glycan')
        assert kwargs.get('type')
        assert kwargs.get('property')
        assert kwargs.get('source') 
        pagename = ".".join(tuple(map(kwargs.get,('glycan','type','property','source'))))
	assert ':' not in pagename
	return pagename

    def key(self):
        assert self.get('type')
        assert self.get('property')
        assert self.get('source') 
        return tuple(map(self.get,('type','property','source')))

    def goodvalue(self):
        return (self.get('value') not in (None,"",[]))

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

	if data.get('type') in ['CrossReference','Motif','Taxonomy','Publication'] or \
           data.get('property') in ['Compositions','Topologies','Saccharides','SubsumedBy','Subsumes']:
            if isinstance(data.get('value'),basestring):
                data['value'] = sorted(map(lambda s: s.strip(),data.get('value').split(';')),key=self.intstrvalue)
        
        return data

    def toTemplate(self,data):

        data = super(Annotation,self).toTemplate(data)
        
	if data.get('value'):
	  if data.get('type') in ['CrossReference','Motif','Taxonomy','Publication'] or \
             data.get('property') in ['Compositions','Topologies','Saccharides','SubsumedBy','Subsumes']:
	    if isinstance(data['value'],list):
		if len(data['value']) > 1:
                    data['value'] = ";".join(map(str,sorted(data['value'],key=self.intstrvalue)))
	            data['multivaluesep'] = ";"
	        else:
	            data['value'] = str(data['value'][0])
	    else:
		data['value'] = str(data['value'])

        return data

class GlycanDataWiki(SMW.SMWSite):
    _name = 'glycandata'

    def __init__(self):
        super(GlycanDataWiki,self).__init__()

    # Assumes accession == pagename, see above
    def get(self,accession):
        g = super(GlycanDataWiki,self).get(accession)
        if g:
            for anpage in self.site.allpages(prefix='%s.'%(accession)):
                an = super(GlycanDataWiki,self).get(anpage.name)
                g.set_annotation(annotation=an)
        return g

    def put(self,g):
        accession = g.get('accession')
        keys2page = dict()
	for anpage in self.site.allpages(prefix='%s.'%(accession)):
            an = super(GlycanDataWiki,self).get(anpage.name)
	    keys2page[an.key()] = anpage
        changed = 0
        for an in g.annotations():
            an.set('glycan',accession)
            if super(GlycanDataWiki,self).put(an):
                changed += 1
	    if an.key() in keys2page:
	        del keys2page[an.key()]
	for anpage in keys2page.values():
	    if anpage.exists:
	        anpage.delete()
                changed += 1
        if super(GlycanDataWiki,self).put(g):
            changed += 1
        return (changed>0)

    def iterglycan(self):
        for pagename in self.itercat('Glycan'):
            m = self.get(pagename)
            yield m
