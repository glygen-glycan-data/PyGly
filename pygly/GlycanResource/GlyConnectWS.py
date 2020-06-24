
from WebServiceResource import WebServiceResource

import os, os.path, re, time
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
	    row['hash'] = hashlib.sha1(";".join(map(lambda t: "%s:%s"%t,sorted(row.items())))).hexdigest().lower()
	    if row['hash'] in seen:
		continue
	    seen.add(row['hash'])
	    # print data
	    yield row

    def taxonomy(self):
	for id in range(300):
	    data = self.query_taxonomy_single(id=id)
	    for l in data.splitlines():
		m = re.search(r'<h1>(.*)</h1>',l)
		if m and 'Oops' not in m.group(1):
		     yield dict(id=id,name=m.group(1))
	    time.sleep(1)

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
