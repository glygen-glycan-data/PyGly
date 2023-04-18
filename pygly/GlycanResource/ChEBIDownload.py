
import sys, os, gzip

from .WebServiceResource import WebServiceResource

class ChEBIDownload(WebServiceResource):
    apiurl = 'http://ftp.ebi.ac.uk/pub/databases/chebi/ontology'
    singlevalues = ['id','name']
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"chebidownload.ini")
        super(ChEBIDownload,self).__init__(**kw)

    def allterms(self):
        data = None
        for l in self.query_obogz():
            l = l.strip()
            if l == "[Term]":
                data = {}
                continue
            if l == "" :
                if data != None:
                    yield data
                    data = None
                continue
            if data != None:
                sl = l.split(': ',1)
                if sl[0] in data:
                    assert sl[0] not in self.singlevalues, sl[0]
                    data[sl[0]].append(sl[1])
                else:
                    if sl[0] in self.singlevalues:
                        data[sl[0]] = sl[1]
                    else:
                        data[sl[0]] = [ sl[1] ]
        if data != None:
            yield data

    def allgtcacc(self):
        for r in self.allterms():
            gtcacc = set()
            for xf in r.get('xref',[]):
                if xf.startswith('GlyTouCan:'):
                    gtcacc.add(xf.split(':',1)[1])
            if len(gtcacc) == 0:
                continue
            for acc in gtcacc:
                yield dict(chebiterm=r['id'],gtcacc=acc,name=r['name'])
