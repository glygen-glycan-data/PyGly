
import sys, os, gzip

from .WebServiceResource import WebServiceResource

class PubChemDownload(WebServiceResource):
    apiurl = 'https://ftp.ncbi.nlm.nih.gov/pubchem'
    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"pubchemdownload.ini")
        super(PubChemDownload,self).__init__(**kw)

    def allrows(self):
        for l in self.query_sidmapgz():
            sl = [ s.strip() for s in l.split('\t') ]
            if sl[1] == 'GlyTouCan Project':
                sid = sl[0]
                gtc = sl[2]
                cid = sl[3] if len(sl) > 3 else None
                yield dict(sid=sid,accession=gtc,cid=cid)

    def allgtc(self):
        for r in self.allrows():
            yield "SID"+r['sid'],r['accession']
            if r['cid']:
                yield "CID"+r['cid'],r['accession']
