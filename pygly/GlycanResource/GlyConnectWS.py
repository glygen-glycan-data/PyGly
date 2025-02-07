
from .WebServiceResource import WebServiceResource

import os, os.path, re, time, urllib
import hashlib 

class GlyConnectWS(WebServiceResource):
    apiurl = 'https://glyconnect.expasy.org/api'

    def __init__(self,**kw):
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glyconnect.ini")
        super(GlyConnectWS,self).__init__(**kw)

    def allrows(self):
        for species in self.taxonomy():
            for row in self.rowsbyspecies(species=species['name']):
                yield row

    def rowsbyspecies(self,species):
        seen = set()
        for data in self.query_glycosylations(species)['results']:
            # print data
            row = dict(taxid=data['taxonomy']['taxonomy_id'],
                       accession=data['structure']['id'],
                           gtcacc=data['structure'].get('glytoucan_id'),
                       compaccession=data['composition']['id'],
                       compgtcacc=data['composition'].get('glytoucan_id'),
                       compbyonic=data['composition'].get('format_byonic'))
            if 'Other' in row['compbyonic']:
                continue
            pmids = set()
            for ref in data['references']:
                if 'pmid' in ref:
                    pmids.add(ref['pmid'])
            row['pmids'] = ";".join(map(str,sorted(pmids)))
            row['hash'] = hashlib.sha1(";".join(["%s:%s"%t for t in sorted(row.items())]).encode()).hexdigest().lower()
            if row['hash'] in seen:
                continue
            seen.add(row['hash'])
            # print data
            yield row

    def taxonomy(self):
        delay = 1
        for id in range(1,300):
            try:
                data = self.query_taxonomy_single(id=id)
            except urllib.error.HTTPError:
                time.sleep(delay)
                continue
            for l in data:
                l = l.decode()
                m = re.search(r'<h1>(.*)</h1>',l)
                if m and 'Oops' not in m.group(1):
                     yield dict(id=id,name=m.group(1))
            time.sleep(delay)

    def alltaxa(self):
        for row in self.allrows():
            if row.get('accession') and row.get('gtcacc') and row.get('taxid'):
                yield "S"+str(row.get('accession')),row.get('gtcacc'),row.get('taxid')
            elif row.get('compaccession') and row.get('compgtcacc') and row.get('taxid'):
                yield "C"+str(row.get('compaccession')),row.get('compgtcacc'),row.get('taxid')

    def gtcbyspecies(self,species):
        seen = set()
        for row in self.rowsbyspecies(species):
            gtcacc = row.get('gtcacc')
            if not gtcacc:
                gtcacc = row.get('compgtcacc')
            if not gtcacc:
                continue
            taxid = row.get('taxid')
            if not taxid:
                continue
            if gtcacc in seen:
                continue
            seen.add(gtcacc)
            yield gtcacc
