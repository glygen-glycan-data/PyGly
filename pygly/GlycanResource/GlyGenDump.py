
from .WebServiceResource import WebServiceResource
import os, os.path, re, traceback, sys

def uniqueify(iterable):
    def wrapper(*args,**kwargs):
        seen = set()
        for it in iterable(*args,**kwargs):
            if it not in seen:
                yield it
                seen.add(it)
    return wrapper

class GlyGenSourceFile(WebServiceResource):
    apiurl = "https://data.glygen.org/ln2data/downloads"
    # verbose = True
    def __init__(self,**kw):
        self._taxidlookup = None
        # kw['delaytime'] = 60
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
        for row in self.glyconnect("fruitfly"):
            yield row

    def glyconnect_covid(self):
        for row in self.glyconnect("sarscov2"):
            yield row

    def glyconnect_yeast(self):
        for row in self.glyconnect("yeast"):
            yield row

    def glyconnect_all(self):
        for sp in ("human","mouse","rat","fruitfly","sarscov2","yeast"):
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
            if row['taxonomy'].get('taxonomy_id'):
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
        "sarscov2": ('SARS coronavirus (SARS-CoV-2 or 2019-nCoV)','2697049'),
        "fruitfly": ('Drosophila melanogaster','7227'),
        "yeast": ('Saccharomyces cerevisiae', '4932'),
        "slimemold": ('Dictyostelium discoideum', '44689'),
        "pig": ('Sus scrofa', '9823')
    }

    @uniqueify
    def glygen_protein_taxa(self):
        for row in self.query_glygen_protein_masterlist():
            taxid = row['tax_id']
            for key in ("uniprotkb_canonical_ac","unreviewed_isoforms","reviewed_isoforms"):
                if row.get(key):
                    acc = row.get(key)
                    yield acc,taxid
                    acc = acc.split('-',1)[0]
                    yield acc,taxid

    unicarbkb_filenames = """ 
        known_legacy_dbogap_dbogap2_glygen
        known_34106099_glygen known_unicarbkb_glycomics_study
        fuzzy_legacy_cph_start_end_positions
        known_legacy_cph_single_positions
        unknown_legacy_human_proteoform_glycosylation_sites_unicarbkb_ogalnac_data
        known_legacy_human_mouse_rat_glygen
        known_34004174_glygen known_33806155_glygen
        known_33167210_glygen known_legacy_reis_gastric_ags
        known_legacy_sarscov2_glygen unknown_20131323_glygen
        fuzzy_simplecell_all_glygen_qc mixed_34089345_glygen
        known_33064451_glygen unknown_33993882_glygen
        fuzzy_legacy_human_proteoform_glycosylation_sites_unicarbkb_qin_isoforms
        unknown_34083726_glygen known_legacy_larsson_xyl
        unknown_34327564_glygen fuzzy_legacy_qin_oglanac
        unknown_34191522_glygen
    """.split()

    badupacc = dict(map(str.split,filter(lambda s: s.strip() != "","""
       P02470	9913
       P02488	9544
       P02505	8797
       P03070	1891767
       P08318	10360
       P17767	12214
       P54939	9031
       Q8BHY1	10090
       Q9QYX7 	10090
       Q9QYX6	10090
       A4D997	9606
       Q9UBGc0Pent0S0P0KDN0HexA0	-
       Q7Z7Gc0Pent0S0P0KDN0HexA0	-
       Q4KMGc0Pent0S0P0KDN0HexA0	-
       #N/A	-
       Q7TQ25	10116
       B4F7C9	10116
       B5DEQ0	10116
       F1LTJ9	10116
       F1LUR0	10116
       F1M8K1	10116
       D3ZWD0	10116
       D3ZGJ7	10116
       G3V3K2	9606
       B7Z1H8	9606
       B4E1H2	9606
       P90489	868565
       K9N5Q8	1263720
       Protein	-
    """.splitlines())))

    def unicarbkb(self,filename,loadtaxid=False):
        if self._taxidlookup == None:
            self._taxidlookup = {}
            for key,(sp,taxid) in self.glygen_species.items():
                self._taxidlookup[sp] = taxid
                self._taxidlookup[key] = taxid
        if loadtaxid and len(set(self._taxidlookup.values())) == len(self.glygen_species):
            self._taxidlookup.update(dict(self.glygen_protein_taxa()))
            for k,v in self.badupacc.items():
                if v == "-":
                    self._taxidlookup[k] = None
                else:
                    self._taxidlookup[k] = v
        for i,row in enumerate(self.query_csvsourcefile(source="unicarbkb",filename=filename)):
            row['filename'] = filename
            row['lineno'] = str(i+2);
            # deal with case inconsistency
            for k in list(row.keys()):
                row[k] = row[k].strip()
                if k.lower() != k:
                    row[k.lower()] = row[k]
                    del row[k]
            thetaxid = None
            if row.get('taxid'):
                thetaxid = row.get('taxid')
            prottaxid = None
            if row.get('protein'):
                val = row.get('protein')
                if val in self._taxidlookup:
                    prottaxid = self._taxidlookup[val]
                elif val in ("","0","NULL"):
                    pass
                elif loadtaxid:
                    raise RuntimeError("Bad protein accession, %s: %s"%(val,row))
            spectaxid = None
            if row.get('species'):
                val = row.get('species')
                if val in self._taxidlookup:
                    spectaxid = self._taxidlookup[val]
                elif val.lower() in self._taxidlookup:
                    spectaxid = self._taxidlookup[val.lower()]
                elif val in ("","0"):
                    pass
                else:
                    raise RuntimeError("Bad species, %s: %s"%(val,row))

            if loadtaxid and len(set(filter(None,[thetaxid,prottaxid,spectaxid]))) > 1:
                # print("WARNING: Species conflict %s:%s, taxid=%s,protein=%s,spectaxid=%s,pmid=%s"%(row['filename'],row['lineno'],thetaxid,prottaxid,spectaxid,row['pmid']),file=sys.stderr)
                # continue
                # sys.exit(1)
                pass

            if thetaxid:
                row['taxid'] = thetaxid
            elif prottaxid:
                row['taxid'] = prottaxid
            elif spectaxid:
                row['taxid'] = spectaxid

            uckbid = None
            try:
                uckbid = int(row['id'])
            except KeyError:
                if 'glytoucan' in row:
                    row['id'] = row['glytoucan']
                else:
                    continue
            except ValueError:
                pass
            if not uckbid:
                # comp_HexNAc1Hex0dHex0NeuAc0NeuGc0Pent0S0P0KDN0HexA0
                if re.search(r'^comp_HexNAc\d+Hex\d+dHex\d+NeuAc\d+NeuGc\d+Pent\d+S\d+P\d+KDN\d+HexA\d+$',row['id']):
                    uckbid = row['id']
                elif re.search(r'^G\d{5}[A-Z]{2}$',row['id']):
                    uckbid = row['id']
                else:
                    continue
            row['id'] = uckbid
            # print(row)
            yield row

    def unicarbkb_all(self,**kw):
        for filename in self.unicarbkb_filenames:
            for row in self.unicarbkb(filename=filename,**kw):
                yield row

    @uniqueify
    def unicarbkb_allgtc(self):
        for row in self.unicarbkb_all():
            if row.get('toucan') and row.get('id'):
                if re.search('^G\d{5}[A-Z]{2}$',row.get('toucan').strip()):
                    yield (row['id'],row.get('toucan').strip())

    @uniqueify
    def unicarbkb_alltaxa(self):
        for row in self.unicarbkb_all(loadtaxid=True):
            if row.get('taxid') and row.get('id'):
                if int(row['taxid']) > 0:
                    yield (row['id'],row['taxid'].strip())

    @uniqueify
    def unicarbkb_allpubs(self):
        for row in self.unicarbkb_all():
            if row.get('id') and row.get('pmid'):
                if int(row['pmid']) > 0:
                    yield (row['id'],row['pmid'].strip())

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
            if row.get('glytoucan_ac'):
                yield row['glytoucan_ac'],row['glytoucan_ac']

    def ncfg_allpubs(self):
        for row in self.ncfg_all():
            if row.get('glytoucan_ac') and row.get('evidence'):
                yield row['glytoucan_ac'],row['evidence']

    def synthetic(self):
        for row in self.query_csvfullurl(url='https://data.glygen.org/ln2releases/data/current/compiled/glycan_synthesized.csv'):
            yield row

    def synthetic_all(self):
        for row in self.synthetic():
            yield row

    def synthetic_allgtc(self):
        for row in self.synthetic_all():
            if row.get('glytoucan_ac'):
                yield row['glytoucan_ac'],row['glytoucan_ac']

    def synthetic_allpubs(self):
        for row in self.synthetic_all():
            if row.get('glytoucan_ac') and row.get('evidence'):
                yield row['glytoucan_ac'],row['evidence']

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
            yield row['glytoucan_ac'],row['glytoucan_ac']

    def mcw_oglcnac(self,name):
        for row in self.query_csvsourcefile(source="mcw_oglcnac",filename=name+"_o_glcnacome_mcw"):
            yield row

    def mcw_oglcnac_human(self):
        for row in self.mcw_oglcnac("human"):
            yield row

    def mcw_oglcnac_mouse(self):
        for row in self.mcw_oglcnac("mouse"):
            yield row

    def mcw_oglcnac_rat(self):
        for row in self.mcw_oglcnac("rat"):
            yield row

    def mcw_oglcnac_fruitfly(self):
        for row in self.mcw_oglcnac("fruitfly"):
            yield row

    def mcw_oglcnac_yeast(self):
        for row in self.mcw_oglcnac("yeast"):
            yield row

    def mcw_oglcnac_all(self):
        for sp in ("human","mouse","rat","fruitfly","yeast"):
            for row in self.mcw_oglcnac(sp):
                yield row

    @uniqueify
    def mcw_oglcnac_allgtc(self):
        for row in self.mcw_oglcnac_all():
            yield "G49108TO","G49108TO"

    @uniqueify
    def mcw_oglcnac_alltaxa(self):
        taxa_lookup = dict()
        for k,v in self.glygen_species.values():
            taxa_lookup[k] = v
        for row in self.mcw_oglcnac_all():
            if 'organism' in row:
                org = row['organism'].strip()
                if org in taxa_lookup:
                    yield "G49108TO",taxa_lookup[org]
                else:
                    assert False, row

    @uniqueify
    def mcw_oglcnac_allpubs(self):
        for row in self.mcw_oglcnac_all():
            if row.get('PMIDS'):
                for pmid in row.get('PMIDS').split(';'):
                    try:
                        pmid = int(pmid)
                        if pmid > 0:
                            yield "G49108TO",pmid
                    except ValueError:
                        pass

    @uniqueify
    def allgtc(self):
        for source in ("unicarbkb","glyconnect","ncfg","gptwiki","sandbox","synthetic","mcw_oglcnac"):
            for row in getattr(self,source+"_allgtc")():
                yield row[0],row[1]

    @uniqueify
    def alltaxa(self):
        for source in ("unicarbkb","glyconnect","gptwiki","mcw_oglcnac"):
            for row in getattr(self,source+"_alltaxa")():
                yield row[0],row[1],source
    @uniqueify
    def allpubs(self):
        for source in ("unicarbkb","glyconnect","ncfg","synthetic","mcw_oglcnac"):
            for row in getattr(self,source+"_allpubs")():
                yield row[0],row[1],source

class GlyGenDataset(WebServiceResource):
    apiurl = "https://data.glygen.org/ln2releases/data/current/reviewed"
    def __init__(self,**kw):
        kw['delaytime'] = 60
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygends.ini")
        self.datasets = {}
        for l in self.datasets_data.splitlines():
            sl = l.split()
            if len(sl) < 3:
                continue
            ds = dict(zip(["dsid","filename","taxid","pmid"],sl))
            self.datasets[ds["dsid"]] = ds
        super(GlyGenDataset,self).__init__(**kw)

    def dataset(self,dsid):
        assert dsid in self.datasets, "Bad GlyGen dataset accession: "+dsid
        for row in self.query_csvdataset(filename=self.datasets[dsid]['filename']):
            if 'taxid' not in row and 'taxid' in self.datasets[dsid]:
                row['taxid'] = self.datasets[dsid]['taxid']
            if row.get('xref_key') == 'protein_xref_pubmed':
                row['pmid'] = row['xref_id']
            if 'pmid' not in row and 'pmids' in row:
                row['pmid'] = row['pmids']
            if 'pmid' not in row and 'pmid' in self.datasets[dsid]:
                row['pmid'] = self.datasets[dsid]['pmid']
            if 'accession' not in row and 'saccharide' in row:
                row['accession'] = row['saccharide']
            yield row

    datasets_data = """
        GLY_000142      human_proteoform_glycosylation_sites_harvard                9606
        GLY_000335      hcv1a_proteoform_glycosylation_sites_literature             11108
        GLY_000632      rat_proteoform_glycosylation_sites_oglcnac_mcw              10116
        GLY_000633      fruitfly_proteoform_glycosylation_sites_oglcnac_mcw         7227
        GLY_000631      mouse_proteoform_glycosylation_sites_oglcnac_mcw            10090
        GLY_000517      human_proteoform_glycosylation_sites_oglcnac_mcw            9606
        GLY_000708      human_proteoform_glycosylation_sites_oglcnac_atlas          9606
        GLY_000709      mouse_proteoform_glycosylation_sites_oglcnac_atlas          10090
        GLY_000710      rat_proteoform_glycosylation_sites_oglcnac_atlas            10116
        GLY_000711      fruitfly_proteoform_glycosylation_sites_oglcnac_atlas       7227
        GLY_000716      human_proteoform_glycosylation_sites_o_gluc                 9606
    """
    # yeast_proteoform_glycosylation_sites_oglcnac_mcw

    @uniqueify
    def dataset_allgtc(self,dsid):
        for row in self.dataset(dsid):
            yield row['accession'],row['accession']
            
    @uniqueify
    def dataset_alltaxa(self,dsid):
        for row in self.dataset(dsid):
            yield row['accession'],row['taxid']
            
    @uniqueify
    def dataset_allpubs(self,dsid):
        for row in self.dataset(dsid):
            if 'pmid' in row:
                for pmid in row['pmid'].split(';'):
                    yield row['accession'],pmid
            
    @uniqueify
    def allgtc(self):
        for dsid in self.datasets:
            for row in self.dataset_allgtc(dsid):
                yield row,dsid

    @uniqueify
    def alltaxa(self):
        for dsid in self.datasets:
            for row in self.dataset_alltaxa(dsid):
                yield row[0],row[1],dsid

    @uniqueify
    def allpubs(self):
        for dsid in self.datasets:
            for row in self.dataset_allpubs(dsid):
                yield row[0],row[1],dsid

class GlyGenFile:

    def __init__(self,**kw):
        self._kw = kw
    
    @uniqueify
    def allgtc(self):
        ggsf = GlyGenSourceFile(**self._kw)
        for row in ggsf.allgtc():
            yield row
        ggds = GlyGenDataset(**self._kw)
        for row in ggds.allgtc():
            yield row

    @uniqueify
    def alltaxa(self):
        ggsf = GlyGenSourceFile(**self._kw)
        for row in ggsf.alltaxa():
            yield row
        ggds = GlyGenDataset(**self._kw)
        for row in ggds.alltaxa():
            yield row

    @uniqueify
    def allpubs(self):
        ggsf = GlyGenSourceFile(**self._kw)
        for row in ggsf.allpubs():
            yield row
        ggds = GlyGenDataset(**self._kw)
        for row in ggds.allpubs():
            yield row
