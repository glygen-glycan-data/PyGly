
import sys, os, gzip
from collections import defaultdict

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
                yield dict(sid=sid,accession=gtc,cid=cid,type='GlyTouCan')
            elif sl[1] == 'ChEBI':
                sid = sl[0]
                gtc = sl[2]
                cid = sl[3] if len(sl) > 3 else None
                yield dict(sid=sid,accession=gtc,cid=cid,type='ChEBI')

    def allgtc(self):
        for r in self.allrows():
            if r['type'] != 'GlyTouCan':
                continue
            yield "SID"+r['sid'],r['accession']
            if r['cid']:
                yield "CID"+r['cid'],r['accession']

    def chebi_tocid(self):
        result = self.query_chebi_annotation()
        npages = result['Annotations']['TotalPages']
        while True:
            for r in result['Annotations']['Annotation']:
                if 'LinkedRecords' not in r:
                    continue
                for cid in r['LinkedRecords']['CID']:
                    yield r['URL'].split('=')[1],cid
            if result['Annotations']['Page'] == npages:
                break
            result = self.query_chebi_annotation(result['Annotations']['Page']+1)

    def allchebigtc(self):
        chebi2pmcid = defaultdict(set)
        pmcid2gtc = defaultdict(set)
        for r in self.allrows():
            if r['type'] not in ("GlyTouCan","ChEBI"):
                continue
            if not r['cid']:
                continue
            if r['type'] == "GlyTouCan":
                pmcid2gtc[r['cid']].add(r['accession'])
            else:
                chebi2pmcid[r['accession']].add(r['cid'])
        for r in self.chebi_tocid():
            chebi2pmcid[r[0]].add(r[1])
        for chebi in chebi2pmcid:
            for cid in chebi2pmcid[chebi]:
                for gtc in pmcid2gtc[cid]:
                    yield chebi,gtc
