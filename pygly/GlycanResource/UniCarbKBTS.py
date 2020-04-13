
import os.path

from TripleStoreResource import TripleStoreResource

class UniCarbKBTS(TripleStoreResource):

    # endpt = "http://130.56.249.35:40935/unicarbkb/query"
    endpt = "http://203.101.226.128:40935/unicarbkb/query"
    defns = "http://rdf.unicarbkb.org/structure/"

    def __init__(self,**kw):
	kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"unicarbkb.ini")
	super(UniCarbKBTS,self).__init__(**kw)

    def alltaxa(self):
	for row in self.query_taxonomy():
	    yield row['accession'],row['taxon']

    def allgtc(self):
	for row in self.query_gtcacc():
	    yield row['accession'],row['glytoucan']

    def allpub(self):
	for row in self.query_publication():
	    if row['pmid'] == '0':
		continue
	    yield row['accession'],row['pmid']
