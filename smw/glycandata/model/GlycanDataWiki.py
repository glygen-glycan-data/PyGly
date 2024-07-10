from __future__ import print_function
__all__ = [ "GlycanDataWiki", "GlycanDataWikiNew", "GlycanDataDiskCache", "GlycanData", "Glycan", "Annotation" ]

import sys, re
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

        # if '_subobjs' in data:
        #     del data['_subobjs']

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

    def add_annotation(self,**kwargs):
        assert kwargs.get('value')
        value = kwargs.get('value')
        del kwargs['value']
        if self.has_annotations(**kwargs):
            values = self.get_annotation_values(**kwargs)
            values.append(value)
        else:
            values = [value]
        self.set_annotation(value=values,**kwargs)

    def delete_annotations(self,**kwargs):
        if not hasattr(self,'_annotations'):
            return
        for an in list(self.annotations(**kwargs)):
            del self._annotations[an.key()]

    def count_annotations(self,**kwargs):
        return len(list(self.annotations(**kwargs)))

    def get_annotation_values(self,property=None,**kwargs):
        return [ str(x) for x in self.get_annotation(property=property,**kwargs).get('value')]

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

    def annotations(self,type=None,property=None,source=None,sourceid=None):
        if hasattr(self,'_annotations'):
            for key,an in sorted(self._annotations.items()):
                if (type == None or an.get('type') == type) and \
                       (property == None or an.get('property') == property) and \
                       (source == None or an.get('source') == source) and \
                       (sourceid == None or an.get('sourceid') == sourceid):
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
            sequence = self.get_annotation_value('WURCS')
            return self.wurcs_format.toGlycan(sequence)
        except (LookupError,GlycanParseError,RuntimeError):
            pass
        try:
            sequence = self.get_annotation_value('GlycoCT')
            return self.glycoct_format.toGlycan(sequence)
        except (LookupError,GlycanParseError,RuntimeError):
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
        pagename = ".".join(tuple(filter(None,map(kwargs.get,('glycan','type','property','source','sourceid')))))
        assert ':' not in pagename
        return pagename

    def key(self):
        assert self.get('type')
        assert self.get('property')
        assert self.get('source') 
        return tuple(filter(None,map(self.get,('type','property','source','sourceid'))))

    def goodvalue(self):
        return (self.get('value') not in (None,"",[],set([])))

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

        if data.get('type') in ['CrossReference','Motif','Taxonomy','Publication','Enzyme','Name'] or \
           data.get('property') in ['Compositions','Topologies','Saccharides','SubsumedBy','Subsumes','Ancestor','Descendant','FullyDetermined','Leaf','SequenceHash','ReducingEnd','GlycanType','GlycanSubtype','HasMonosaccharide'] or \
           data.get('property').endswith(' Evidence'):
            if isinstance(data.get('value'),str):
                data['value'] = sorted(map(lambda s: s.strip(),data.get('value').split(';')),key=self.intstrvalue)
        
        return data

    def toTemplate(self,data):

        data = super(Annotation,self).toTemplate(data)
        
        if data.get('value'):
          if data.get('type') in ['CrossReference','Motif','Taxonomy','Publication','Enzyme','Name'] or \
             data.get('property') in ['Compositions','Topologies','Saccharides','SubsumedBy','Subsumes','Ancestor','Descendant','FullyDetermined','Leaf','SequenceHash','ReducingEnd','GlycanType','GlycanSubtype','HasMonosaccharide'] or \
             data.get('property').endswith(' Evidence'):
            if isinstance(data['value'],list) or isinstance(data['value'],set):
                if len(set(data['value'])) > 1:
                    data['value'] = ";".join(map(str,sorted(set(data['value']),key=self.intstrvalue)))
                    data['multivaluesep'] = ";"
                else:
                    data['value'] = str([ x for x in data['value']][0])
            else:
                data['value'] = str(data['value'])

        return data

class GlycanDataWikiNew(SMW.SMWSite):
    _name = 'glycandata'

    def get(self,accession):
        g = super(GlycanDataWikiNew,self).get(accession)
        if g:
            for so in g.get('_subobjs'):
                g.set_annotation(annotation=so)
        return g

    def put(self,g):
        accession = g.get('accession')
        # annotationpages = []
        # for anpage in self.site.allpages(prefix='%s.'%(accession)):
        #     annotationpages.append(anpage)
        g.set('_subobjs',[])
        for an in g.annotations():
            an.set('glycan',accession)
            g.append('_subobjs',an)
        g.sort('_subobjs',Annotation.key)
        # for anpage in annotationpages:
        #     if anpage.exists:
        #         anpage.delete()
        if super(GlycanDataWikiNew,self).put(g):
            return True
        return False

    def iterglycan(self):
        for pagename in self.itercat('Glycan'):
            m = self.get(pagename)
            yield m


class GlycanDataWiki(SMW.SMWSite):
    _name = 'glycandata'

    # Assumes accession == pagename, see above
    def get(self,accession):
        g = super(GlycanDataWiki,self).get(accession)
        if g:
            for anpagename in self.site.allpages(prefix='%s.'%(accession),generator=False):
                an = super(GlycanDataWiki,self).get(anpagename)
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

    def delete(self,acc):
        super(GlycanDataWiki,self).delete(acc)
        for anpagename in self.site.allpages(prefix='%s.'%(acc),generator=False):
             super(GlycanDataWiki,self).delete(anpagename)

    def iterglycan(self):
        for pagename in self.itercat('Glycan'):
            m = self.get(pagename)
            yield m

    def iterglycanid(self,**kwargs):
        if kwargs.get('property') == None and \
           kwargs.get('type') == None and \
           kwargs.get('source') == None:
            regex = kwargs.get('regex',r'^(G\d{5}[A-Z]{2})$')
        else:
            values = [kwargs.get('type',r'[^.]+'),
                      kwargs.get('property',r'[^.]+'),
                      kwargs.get('source',r'[^.]+')]
            regex = r'^(G\d{5}[A-Z]{2})\.' + r'\.'.join(values) + '(\.[^.]+)?$'
        regex = re.compile(regex)
        for pagename in self.site.allpages(generator=False):
            m = regex.search(pagename)
            if m:
                yield m.group(1)

    def iterannotation(self):
        for pagename in self.itercat('Annotation'):
            m = self.get(pagename)
            yield m

import os.path, shutil, os, glob

class GlycanDataDiskCache(object):

    def __init__(self,location=".cache"):
        self.path = '/' + os.path.abspath(location).strip('/')
        assert os.path.isdir(os.path.split(self.path)[0])
        if not os.path.isdir(self.path):
            os.makedirs(self.path)
        print("Connected to disk cache: %s"%(self.path,),file=sys.stderr)

    def acc2path(self,acc):
        assert re.search(r'^G\d{5}[A-Z]{2}$',acc), acc
        return "/".join([self.path,acc[1],acc[2],acc[3],acc])

    def acc2glypath(self,acc):
        return "/".join([self.acc2path(acc),acc+".txt"])

    def acc2annpath(self,acc):
        return "/".join([self.acc2path(acc),acc+".*.txt"])

    def get(self,accession):

        glypath = self.acc2glypath(accession)
        if not os.path.exists(glypath):
            return None
        g = Glycan(**SMW.SMWSite.parse_template_text(accession,open(glypath).read()))
        for anpath in glob.glob(self.acc2annpath(accession)):
            an = Annotation(**SMW.SMWSite.parse_template_text(accession,open(anpath).read()))
            g.set_annotation(annotation=an)
        return g

    def write(self,filename,text):
        data = None
        if os.path.exists(filename):
            data = open(filename).read()
        if text == data:
            return False
        wh = open(filename,'w')
        wh.write(text)
        wh.close()
        return True

    def put(self,g):

        changes = 0

        accession = g.get('accession')
        path = self.acc2path(accession)
        if not os.path.exists(path):
            os.makedirs(path)

        glypath = self.acc2glypath(accession)
        if self.write(glypath,g.astemplate()):
            changes += 1

        ankeys = set()
        for an in g.annotations():
            an.set('glycan',accession)
            ankey = ".".join(an.key())
            ankeys.add(ankey)
            anpath = "/".join([path,accession+"."+ankey+".txt"])
            if self.write(anpath,an.astemplate()):
                changes += 1

        for anfile in glob.glob(self.acc2annpath(accession)):
            ankey = os.path.split(anfile)[1].rsplit('.',1)[0].split('.',1)[1]
            if ankey not in ankeys:
                os.unlink(anfile)
                changes += 1

        return (changes>0)
        
    def delete(self,accession):
        
        path = self.acc2path(accession)
        shutil.rmtree(path)

    def iterglycan(self,fr=None,to=None):
        for acc in self.iterglycanid(fr,to):
            yield self.get(acc)
        
    def iterglycanid(self,fr=None,to=None):

        for root, dirs, files in os.walk(self.path):
            dirs.sort()
            any = False
            for d in dirs:
                if re.search(r'^G\d{5}[A-Z]{2}$',d):
                    any = False
                    if ((fr==None) or (fr <= d)) and ((to==None) or (to >= d)):
                        yield d
            if any:
                dirs = []
            
    def tocache(self,wiki):

        for g in wiki.iterglycan():
            if self.put(g):
                print(g.get('accession'))

    def towiki(self,wiki,fr=None,to=None):

        for g in self.iterglycan(fr,to):
            if wiki.put(g):
                # wiki.refresh(g)
                print(g.get('accession'))

def GlycanData():
    if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
        dir = sys.argv[1]
        sys.argv.pop(1)
        return GlycanDataDiskCache(dir)
    else:
        return GlycanDataWikiNew()
