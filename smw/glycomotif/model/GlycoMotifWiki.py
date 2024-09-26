
from __future__ import print_function
from past.builtins import basestring

__all__ = [ "GlycoMotifWiki", "Collection", "Motif", "Publication", "Keyword", "Enzyme",
            "GlyTouCanMotif", "AllMotif", "CCRCMotif", "GlycoEpitopeMotif",
            "GlydinMotif",
            "GlydinCummingsMotif", "GlydinHayesMotif", "GlydinCermavMotif",
            "GlydinSugarbindMotif", "GlydinBioligoMotif",
            "UniCarbMotif", "GlyGenMotif", "GlyGenGlycanDict"]

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
    alignmentvalues = [
        "Substructure",
        "Core",
        "Whole-Glycan",
        "Nonreducing-End"
    ]

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
            data['name'] = filter(None,map(lambda s: s.strip(),data.get('name').split('\n')))

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
        
        # redend_alignments is comma separated
        # if isinstance(data.get('redend_alignments'), basestring):
        #     data['redend_alignments'] = map(lambda s: s.strip(), data.get('redend_alignments').split(','))
        # data['redend_alignments'] = sorted(set(data['redend_alignments']))
        
        # other_alignments is comma separated (these two keys should partition the alignments)
        # if isinstance(data.get('other_alignments'), basestring):
        #     data['other_alignments'] = map(lambda s: s.strip(), data.get('other_alignments').split(','))
        # data['other_alignments'] = sorted(set(data['other_alignments'])-set(data['redend_alignments']))
        
        if isinstance(data.get("displayhgv"), basestring):
            data["displayhgv"] = self.asboolean(data.get("displayhgv"))

        # redend is a list of booleans, sorted, so behaves as set
        if isinstance(data.get('redend'),basestring):
            data['redend'] = sorted(map(self.asboolean,map(lambda s: s.strip(),data.get('redend').split(','))))
        elif data.get('redend') in (True,False):
            data['redend'] = [ data.get('redend') ]
        elif data.get('redend') != None:
            data['redend'] = sorted(map(self.asboolean,data.get('redend')))

        # alignment is a list of strings, sorted, so behaves as set
        if isinstance(data.get('alignment'), basestring):
            data['alignment'] = map(lambda s: s.strip(), data.get('alignment').split(','))
            data['alignment'] = sorted(filter(lambda x: x in self.alignmentvalues, data['alignment']))

        if isinstance(data.get('keyword'),basestring):
            data['keyword'] = set(map(lambda s: s.strip(),data.get('keyword').split(';')))

        for k in data:
            if ('enzyme' in k or 'sandbox' in k) and isinstance(data.get(k),basestring):
                data[k] = set(map(lambda s: s.strip(),data.get(k).split(';')))

        if isinstance(data.get('motifcanonres'),basestring):
            data['motifcanonres'] = set(map(lambda s: s.strip(),data.get('motifcanonres').split(';')))

        if isinstance(data.get('dbxref'),basestring):
            data['dbxref'] = set(map(lambda s: tuple(s.strip().split(':')),data.get('dbxref').split(';')))

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
            data['name'] = "\n".join(filter(None,data['name']))

        if 'keyword' in data:
            data['keyword'] = ";".join(filter(None,sorted(set(map(lambda s: s.strip(),data['keyword'])))))

        for k in data:
            if 'enzyme' in k or 'sandbox' in k:
                data[k] = ";".join(filter(None,sorted(set(map(lambda s: s.strip(),data[k])))))

        if 'dbxref' in data:
            data['dbxref'] = ";".join(map(lambda t: "%s:%s"%t,sorted(data['dbxref'])))

        if 'sameas' in data:
            data['sameas'] = ",".join(sorted(data['sameas']))

        if 'aglycon' in data:
            data['aglycon'] = ",".join(sorted(data['aglycon']))

        if 'redend' in data:
            data['redend'] = ",".join(sorted(("true" if redend else "false") for redend in data['redend']))

        if 'topology' in data:
            data['topology'] = ",".join(data['topology'])

        # rea = set()
        # if 'redend_alignments' in data:
        #     # rea = sorted(set(data['redend_alignments']))
        #     # data['redend_alignments'] = ",".join(rea)
        #     del data['redend_alignments']

        # if 'other_alignments' in data:
        #     # oa = sorted(set(data['other_alignments'])-rea)
        #      # data['other_alignments'] = ",".join(oa)
        #     del data['other_alignments']

        if "displayhgv" in data:
            data["displayhgv"] = ("true" if data["displayhgv"] else "false")

        if 'alignment' in data:
            data['alignment'] = ",".join(sorted(data['alignment']))

        # if 'wurcs' in data:
        #     data['wurcs'] = "<pre>" + data['wurcs'] + "</pre>"

        # if 'glycoct' in data:
        #     data['glycoct'] = "<pre>" + data['glycoct'] + "</pre>"

        return data


import findpygly
from pygly.GlycanResource import GlyTouCan

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
                self.gtc = GlyTouCan(usecache=False,prefetch=False)
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

class GlyGenMotif(GlyTouCanMotif):
    id = 'GGM'
    def __init__(self,**kwargs):
        assert kwargs.get('accession') != None
        assert kwargs.get('glytoucan') != None
        super(GlyGenMotif,self).__init__(**kwargs)

class GlyGenGlycanDict(GlyTouCanMotif):
    id = 'GGD'
    def __init__(self,**kwargs):
        assert kwargs.get('accession') != None
        assert kwargs.get('glytoucan') != None
        super(GlyGenMotif,self).__init__(**kwargs)

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

class Publication(SMW.SMWClass):
    template = 'Publication'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('authors') and kwargs.get('title') and kwargs.get('year')
        authlist = kwargs.get('authors').split(None,3)
        if len(authlist) > 2 and authlist[2] != "":
            authslug = authlist[0]+' et al.'
        else:
            authslug = authlist[0]
        year = kwargs.get('year')
        titlewords = kwargs.get('title').split(None,10)
        titleslug = ' '.join(titlewords[:-1])
        if titlewords[-1]:
            titleslug += ' ...'
        return ("%s (%s) %s"%(authslug,year,titleslug)).replace(" ",'_').replace("[","").replace(']','')

    def toPython(self,data):
        data = super(Publication,self).toPython(data)

        if isinstance(data.get('citedby'),basestring):
            data['citedby'] = data.get('citedby').split(';')

        return data

    def toTemplate(self,data):
        data = super(Publication,self).toTemplate(data)

        if 'citedby' in data:
            data['citedby'] = ';'.join(map(str,data.get('citedby')))

        return data

class Keyword(SMW.SMWClass):
    template = 'Keyword'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('keyword')
        return kwargs.get('keyword')

class Enzyme(SMW.SMWClass):
    template = 'Enzyme'

    def toPython(self,data):
        data = super(Enzyme,self).toPython(data)

        for key in ('phenotype','disease','tissue','celltype'):
            if isinstance(data.get(key),basestring):
                data[key] = data.get(key,'').split(';')

        for key in ('ortholog',):
            if isinstance(data.get(key),basestring):
                data[key] = data.get(key,'').split(',')

        return data

    def toTemplate(self,data):
        data = super(Enzyme,self).toTemplate(data)

        for key in ('phenotype','disease','tissue','celltype'):
            if key in data:
                data[key] = ';'.join(sorted(map(str,data.get(key,''))))

        if 'ortholog' in data:
            data['ortholog'] = ",".join(sorted(data['ortholog']))

        return data

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('genename')
        return kwargs.get('genename')

class GlycoMotifWiki(SMW.SMWSite):
    _name = 'glycomotif'

    def get(self,pagename=None,collection=None,accession=None):
        if pagename:
            return super(GlycoMotifWiki,self).get(pagename)
        return super(GlycoMotifWiki,self).get(Motif.pagename(collection=collection,accession=accession))

    def itermotif(self,collection=None):
        if collection != None and not isinstance(collection,basestring):
            collection = collection.id
        if collection:
            regex = r'^%s\.'%(collection,)
        else:
            regex = r'^[A-Z]{2,4}\.'
        for page in self.iterpages(regex=regex,include_categories=['Motif']):
            m = self.get(page.name)
            if not collection or m.get('collection') == collection:
                yield m

    def itercollection(self):
        for pagename in self.itercat('Collection'):
            yield self.get(pagename)

    def iterenzyme(self,species=None):
        for pagename in self.itercat('Enzyme'):
            e = self.get(pagename)
            if not species or e.get('species') == species:
                yield e

if __name__ == "__main__":
    import sys, os

    smw = GlycoMotifWiki()
    print(smw.get(collection='GlyTouCan',accession='G00026MO'))
    print(smw.get(collection='GlyTouCan',accession='G00029MO'))

    motif = GlyTouCanMotif(accession="G00028MO",
                           name="N-Glycan high mannose")
    smw.put(motif)
    print(smw.get(collection='GlyTouCan',accession='G00028MO'))
