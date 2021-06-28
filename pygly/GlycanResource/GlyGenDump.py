
from .WebServiceResource import WebServiceResource
import os, os.path, re

def uniqueify(iterable):
    def wrapper(*args,**kwargs):
        seen = set()
        for it in iterable(*args,**kwargs):
            if it not in seen:
                yield it
                seen.add(it)
    return wrapper

# https://docs.google.com/document/d/1GbfX3JLWPP56cAcBFSGvwAgn4PpG3pc4qnvkXsMbwes

class GlyGenSourceFile(WebServiceResource):
    apiurl = "https://data.glygen.org/ln2downloads"
    # verbose = True
    def __init__(self,**kw):
        self._taxidlookup = None
        kw['delaytime'] = 60
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygendump.ini")
        super(GlyGenSourceFile,self).__init__(**kw)

    def glyconnect(self,name):
        for row in self.query_jsonsourcefile(source="glyconnect",filename="glyconnect_"+name)['results']:
            yield row

    def glyconnect_human(self):
        for row in self.glyconnect("human"):
            yield row

    def glyconnect_mouse(self):
        for row in self.glyconnect("mouse"):
            yield row

    def glyconnect_rat(self):
        for row in self.glyconnect("rat"):
            yield row

    def glyconnect_dros(self):
        for row in self.glyconnect("drosophila"):
            yield row

    def glyconnect_covid(self):
        for row in self.glyconnect("sarscov2"):
            yield row

    def glyconnect_all(self):
        for sp in ("human","mouse","rat","drosophila","sarscov2"):
            for row in self.glyconnect(sp):
                yield row

    @uniqueify
    def glyconnect_allgtc(self):
        for row in self.glyconnect_all():
            if row['structure'].get('glytoucan_id'):
                yield ("S"+str(row['structure']['id']),row['structure'].get('glytoucan_id'))
            if row['composition'].get('glytoucan_id'):
                yield ("C"+str(row['composition']['id']),row['composition'].get('glytoucan_id'))

    @uniqueify
    def glyconnect_alltaxa(self):
        for row in self.glyconnect_all():
            yield ("S"+str(row['structure']['id']),row['taxonomy'].get('taxonomy_id'))
            yield ("C"+str(row['composition']['id']),row['taxonomy'].get('taxonomy_id'))

    @uniqueify
    def glyconnect_allpubs(self):
        for row in self.glyconnect_all():
            for pub in row['references']:
                if pub.get('pmid'):
                    yield ("S"+str(row['structure']['id']),pub['pmid'])
                    yield ("C"+str(row['composition']['id']),pub['pmid'])

    @uniqueify
    def glygen_protein_taxa(self):
        for sp,taxid in (("human","9606"),("mouse","10090"),("rat","10116")):
            for row in self.query_protein_masterlist(sp):
                for key in ("uniprotkb_canonical_ac","unreviewed_isoforms","reviewed_isoforms"):
                    if row.get(key):
                        acc = row.get(key)
                        yield acc,taxid
                        acc = acc.split('-',1)[0]
                        yield acc,taxid

    def unicarbkb(self,name,loadtaxid=False):
        taxid = None
        if name == "sarscov2":
            taxid = '2697049'
        elif loadtaxid and self._taxidlookup == None:
            self._taxidlookup = dict(self.glygen_protein_taxa())
        for row in self.query_csvsourcefile(source="unicarbkb",filename="unicarbkb_"+name):
            if taxid:
                row['taxid'] = taxid
            elif self._taxidlookup and self._taxidlookup.get(row['Protein']):
                row['taxid'] = self._taxidlookup[row['Protein']]
            # Permit integers or comp_ strings in the Id column
	    try:
		uckbid = int(row['Id'])
		if uckbid < 1 or uckbid > 100000:
		    continue
	    except ValueError:
		# comp_HexNAc1Hex0dHex0NeuAc0NeuGc0Pent0S0P0KDN0HexA0
		if not re.search(r'^comp_HexNAc\d+Hex\d+dHex\d+NeuAc\d+NeuGc\d+Pent\d+S\d+P\d+KDN\d+HexA\d+$',row['Id'].strip()):
		    continue
            yield row

    def unicarbkb_human_mouse_rat(self):
        for row in self.unicarbkb("human_mouse_rat"):
            yield row

    def unicarbkb_covid(self):
        for row in self.unicarbkb("sarscov2"):
            yield row

    def unicarbkb_all(self,**kw):
        for sp in ("human_mouse_rat","sarscov2"):
            for row in self.unicarbkb(sp,**kw):
                yield row

    @uniqueify
    def unicarbkb_allgtc(self):
        for row in self.unicarbkb_all():
            if row.get('Toucan') and row.get('Id'):
		if re.search('^G\d{5}[A-Z]{2}$',row.get('Toucan').strip()):
                    try:
			dummy = int(row['Id'])
		    except:
			continue
                    yield (row['Id'],row.get('Toucan').strip())

    @uniqueify
    def unicarbkb_alltaxa(self):
        for row in self.unicarbkb_all(loadtaxid=True):
            if row.get('taxid') and row.get('Id'):
                yield (row['Id'],row['taxid'])

    @uniqueify
    def unicarbkb_allpubs(self):
        for row in self.unicarbkb_all():
            if row.get('Id') and row.get('Pmid') and int(row['Pmid']) > 0:
                yield (row['Id'],row['Pmid'])

    def ncfg(self,name):
        for row in self.query_csvsourcefile(source="ncfg",filename=name+"_ncfg"):
            yield row

    def ncfg_synthetic(self):
        for row in self.ncfg("glycan_evidence"):
            yield row

    def ncfg_all(self):
        for sp in ("glycan_evidence",):
            for row in self.ncfg(sp):
                yield row

    def ncfg_allgtc(self):
        for row in self.ncfg_all():
            yield None,row['glytoucan_ac']

    def gptwiki(self,name):
        for row in self.query_csvsourcefile(source="gptwiki",filename=name):
            yield row

    def gptwiki_sites(self):
        for row in self.gptwiki("glycosites"):
            yield row

    def gptwiki_all(self):
        for sp in ("glycosites",):
            for row in self.gptwiki(sp):
                yield row

    @uniqueify
    def gptwiki_allgtc(self):
        for row in self.gptwiki_all():
            yield (row['GlyTouCan'],row['GlyTouCan'])

    @uniqueify
    def gptwiki_alltaxa(self):
        for row in self.gptwiki_all():
            yield (row['GlyTouCan'],"9606")

    def sandbox(self,name):
        for row in self.query_csvsourcefile(source="sandbox",filename="sandbox_"+name):
            yield row

    def sandbox_glycans(self):
        for row in self.sandbox("accessions"):
            yield row

    def sandbox_all(self):
        for sp in ("accessions",):
            for row in self.sandbox(sp):
                yield row

    @uniqueify
    def sandbox_allgtc(self):
        for row in self.sandbox_all():
            yield None,row['glytoucan_ac']

    @uniqueify
    def allgtc(self):
        for source in ("unicarbkb","glyconnect","ncfg","gptwiki","sandbox"):
            for row in getattr(self,source+"_allgtc")():
                yield row[1]

