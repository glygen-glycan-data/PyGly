
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
    apiurl = "https://data.glygen.org/ln2data/downloads"
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

    glygen_species = {
	"human": ('Homo sapiens','9606'),
        "mouse": ('Mus musculus','10090'),
        "rat": ('Rattus norvegicus','10116'),
        "hcv1a": ('Hepatitis C virus (genotype 1a, isolate H)','11108'),
        "hcv1b": ('Hepatitis C virus (genotype 1b, isolate Japanese)','11116'),
        "sarscov1": ('SARS coronavirus (SARS-CoV-1)','694009'),
        "sarscov2": ('SARS coronavirus (SARS-CoV-2 or 2019-nCoV)','2697049')
    }

    @uniqueify
    def glygen_protein_taxa(self):
        for spname in self.glygen_species:
            taxid = self.glygen_species[spname][1]
            for row in self.query_protein_masterlist(spname):
                for key in ("uniprotkb_canonical_ac","unreviewed_isoforms","reviewed_isoforms"):
                    if row.get(key):
                        acc = row.get(key)
                        yield acc,taxid
                        acc = acc.split('-',1)[0]
                        yield acc,taxid

    unicarbkb_filenames = {
        'human_mouse_rat': 'known_legacy_human_mouse_rat_glygen',
        'sarscov2': 'known_legacy_sarscov2_glygen',
        'sarscov1': 'known_33064451_glygen',
        'glycomics_study': 'known_unicarbkb_glycomics_study'
    }

    def unicarbkb(self,name,loadtaxid=False):
        taxid = None
        if name in self.glygen_species:
            taxid = self.glygen_species[name][1]
        elif loadtaxid and self._taxidlookup == None:
            self._taxidlookup = dict(self.glygen_protein_taxa())
            for sp,taxid in self.glygen_species.values():
                self._taxidlookup[sp] = taxid
        filename = self.unicarbkb_filenames[name]
        for row in self.query_csvsourcefile(source="unicarbkb",filename=filename):
            if taxid:
                row['taxid'] = taxid
            elif self._taxidlookup:
                for key in ('protein','species'):
                    if row.get(key) and self._taxidlookup.get(row[key].strip()):
                        row['taxid'] = self._taxidlookup[row[key].strip()]
                        break
            uckbid = None
            try:
                uckbid = int(row['Id'])
            except ValueError:
                pass
            if not uckbid:
                # comp_HexNAc1Hex0dHex0NeuAc0NeuGc0Pent0S0P0KDN0HexA0
                if re.search(r'^comp_HexNAc\d+Hex\d+dHex\d+NeuAc\d+NeuGc\d+Pent\d+S\d+P\d+KDN\d+HexA\d+$',row['Id'].strip()):
                    uckbid = row['Id'].strip()
                elif re.search(r'^G\d{5}[A-Z]{2}$',row['Id'].strip()):
                    uckbid = row['Id'].strip()
                else:
                    continue
            row['Id'] = uckbid
            yield row

    def unicarbkb_human_mouse_rat(self):
        for row in self.unicarbkb("human_mouse_rat"):
            yield row

    def unicarbkb_sarscov2(self):
        for row in self.unicarbkb("sarscov2"):
            yield row

    def unicarbkb_sarscov1(self):
        for row in self.unicarbkb("sarscov1"):
            yield row

    def unicarbkb_glycomics_study(self):
        for row in self.unicarbkb("glycomics_study"):
            yield row

    def unicarbkb_all(self,**kw):
        for sp in ("human_mouse_rat","sarscov2","sarscov1","glycomics_study"):
            for row in self.unicarbkb(sp,**kw):
                yield row

    @uniqueify
    def unicarbkb_allgtc(self):
        for row in self.unicarbkb_all():
            if row.get('toucan') and row.get('Id'):
                if re.search('^G\d{5}[A-Z]{2}$',row.get('toucan').strip()):
                    yield (row['Id'],row.get('toucan').strip())

    @uniqueify
    def unicarbkb_alltaxa(self):
        for row in self.unicarbkb_all(loadtaxid=True):
            if row.get('taxid') and row.get('Id'):
                if int(row['taxid']) > 0:
                    yield (row['Id'],row['taxid'].strip())

    @uniqueify
    def unicarbkb_allpubs(self):
        for row in self.unicarbkb_all():
            if row.get('Id') and row.get('pmid'):
                if int(row['pmid']) > 0:
                    yield (row['Id'],row['pmid'].strip())

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

