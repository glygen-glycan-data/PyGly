
import csv
import urllib2

class UniCarbKBDump(object):
    dumpfileurl = "https://gitlab.com/matthew.campbell1980/Unicarb-Glygen/-/raw/master/data_files/unicarbkb/DATA_RELEASE/STABLE/mammalian/%s.csv"
    species2taxa = {'human': '9606', 'mouse': '10090', 'rat': '10116'}
    species2filename = {'human': 'human06022020', 'mouse': 'mouse06022020', 'rat': 'rat06022020'}

    def records(self):
	hdr = {'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.11 (KHTML, like Gecko) Chrome/23.0.1271.64 Safari/537.11',
               'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
               'Accept-Charset': 'ISO-8859-1,utf-8;q=0.7,*;q=0.3',
               'Accept-Encoding': 'none',
               'Accept-Language': 'en-US,en;q=0.8',
               'Connection': 'keep-alive'}
	for species in ('human','mouse','rat'):
            url = self.dumpfileurl%(self.species2filename[species],)
	    req = urllib2.Request(url, headers=hdr)
            for row in csv.DictReader(urllib2.urlopen(req)):
		row['taxid'] = self.species2taxa[species]
		yield row

    def alltaxa(self):
	seen = set()
	for row in self.records():
	    data = (row['Id'],row['taxid'])
	    if data in seen:
		continue
	    seen.add(data)
	    yield data

    def allgtc(self):
	seen = set()
	for row in self.records():
	    if not row['Toucan']:
		continue
	    data = (row['Id'],row['Toucan'])
	    if data in seen:
		continue
	    seen.add(data)
	    yield data

    def allpub(self):
	seen = set()
	for row in self.records():
	    if not row['Pmid'] or row['Pmid'] == "0":
		continue
	    data = (row['Id'],row['Pmid'])
	    if data in seen:
		continue
	    seen.add(data)
	    yield data
