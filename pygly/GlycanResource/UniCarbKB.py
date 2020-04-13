
from UniCarbKBDump import UniCarbKBDump
from UniCarbKBTS import UniCarbKBTS

def union(methodname):
    def wrapper(self):
	seen = set()
	for cls in (UniCarbKBTS,UniCarbKBDump):
	    try:
	        for row in getattr(cls,methodname)(self):
	            if row in seen:
		        continue
		    seen.add(row)
		    yield row
	    except IOError:
		pass
    return wrapper

from collections import defaultdict

class UniCarbKB(UniCarbKBDump,UniCarbKBTS):

    alltaxa = union("alltaxa")
    allgtc = union("allgtc")
    allpub = union("allpub")

    def gtcbytaxa(self,taxon):
	accmap = defaultdict(set)
	for acc,gtc in self.allgtc():
	    accmap[acc].add(gtc)
	seen = set()
	for acc,taxid in self.alltaxa():
	    if int(taxid) == int(taxon):
		for gtc in accmap[acc]:
		    if gtc not in seen:
			yield gtc
			seen.add(gtc)
