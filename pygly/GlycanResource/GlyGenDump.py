
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
    apiurl = "https://data.glygen.org/ln2downloads"
    # verbose = True
    def __init__(self,**kw):
        self._taxidlookup = None
        # kw['delaytime'] = 60
        kw['iniFile'] = os.path.join(os.path.dirname(os.path.realpath(__file__)),"glygendump.ini")
        super(GlyGenSourceFile,self).__init__(**kw)

    def section2filename(self,section):
        return section

    def json2rows(self,json):
        return json

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

    def jsonrows(self,section,filename=None):
        if not filename:
            filename = section
        for row in self.json2rows(self.query_jsonsourcefile(source=self.source,filename=filename)):
            for k,v in list(row.items()):
                if isinstance(v,str):
                    row[k] = v.strip()
                if not row[k]:
                    del row[k]
            yield row

    def csvrows(self,section,filename=None,path=None):
        if filename:
            path = self.source + "/current/" + filename
        for row in self.query_csvfullurl(url=self.apiurl+'/'+path+".csv"):
            for k,v in list(row.items()):
                if isinstance(v,str):
                    row[k] = v.strip()
                if not row[k]:
                    del row[k]
            yield row

    def allrows(self,**kw):
        sections = self.sections.split()
        if not hasattr(self,'glygen_sourceid'):
            ggsids = [None]*len(sections)
        elif self.glygen_sourceid == None:
            ggsids = [None]*len(sections)
        else:
            ggsids = self.glygen_sourceid.split()
        assert(len(sections) == len(ggsids))
        if not hasattr(self,'glygen_pubmedid'):
            ggpmids = [None]*len(sections)
        else:
            ggpmids = self.glygen_pubmedid.split()
        assert(len(sections) == len(ggpmids))
        for sec,ggsid,ggpmid in zip(sections,ggsids,ggpmids):
            for row in self.rows(sec,**kw):
                if ggsid:
                    row['_glygen_sourceid'] = ggsid
                if ggpmid:
                    row['_glygen_pubmedid'] = ggpmid
                yield row

    gtcfield = 'glytoucan_ac'
    gtcfixed = None
    taxidfield = 'taxid'
    taxidfixed = None
    pmidfield = 'pmid'
    pmidfixed = None

    def getgtc(self,row):
        if row.get(self.gtcfield,self.gtcfixed):
            yield row.get(self.gtcfield,self.gtcfixed),row.get(self.gtcfield,self.gtcfixed)

    def gettaxa(self,row):
        if row.get(self.gtcfield,self.gtcfixed) and row.get(self.taxidfield,self.taxidfixed):
            yield row.get(self.gtcfield,self.gtcfixed),row.get(self.gtcfield,self.gtcfixed),row.get(self.taxidfield,self.taxidfixed),row.get(self.pmidfield,self.pmidfixed)

    def getpubs(self,row):
        if row.get(self.gtcfield,self.gtcfixed) and row.get(self.pmidfield,self.pmidfixed):
            yield row.get(self.gtcfield,self.gtcfixed),row.get(self.gtcfield,self.gtcfixed),row.get(self.pmidfield,self.pmidfixed)

    @uniqueify
    def allgtc(self):
        for row in self.allrows():
            for acc,gtc in self.getgtc(row):
                yield acc,gtc,self.glygen_source,row.get('_glygen_sourceid')

    @uniqueify
    def alltaxa(self):
        for row in self.allrows():
            for acc,gtc,taxid,pmid in self.gettaxa(row):
                yield acc,gtc,taxid,self.glygen_source,row.get('_glygen_sourceid')
                if pmid:
                    yield acc,gtc,taxid,"PubMed",pmid
                if row.get('_glygen_pubmedid'):
                    yield acc,gtc,taxid,"PubMed",row['_glygen_pubmedid']

    @uniqueify
    def allpubs(self):
        for row in self.allrows():
            for acc,gtc,pmid in self.getpubs(row):
                yield acc,gtc,pmid,self.glygen_source,row.get('_glygen_sourceid')

    def allsourcegtc(self):
        for src in self.sources:
            ggsf = eval(src+"SourceFile")(verbose=self._verbose)
            for row in ggsf.allgtc():
                yield row

    def allsourcetaxa(self):
        for src in self.sources:
            ggsf = eval(src+"SourceFile")(verbose=self._verbose)
            for row in ggsf.alltaxa():
                yield row

    def allsourcepubs(self):
        for src in self.sources:
            ggsf = eval(src+"SourceFile")(verbose=self._verbose)
            for row in ggsf.allpubs():
                yield row

    sources = """
        GlyConnect
        UniCarbKB
        NCFG
        Synthetic
        GPTWiki
        MCWOGlcNAc
        OGlcNAcAtlas
        Harvard
        OGluC
        Literature
        BiomarkerDB
        EMBL
    """.split()

class GlyConnectSourceFile(GlyGenSourceFile):
    source = 'glyconnect'
    glygen_source = "GlyConnect"
    glygen_sourceid = None
    sections = 'human mouse rat fruitfly sarscov2 yeast dicty pig'

    def json2rows(self,json):
        return json['results']

    def rows(self,section):
        return self.jsonrows(section=section,filename=self.source + "_" + section)

    def getgtc(self,row):
        if row['structure'].get('glytoucan_id'):
            yield "S"+str(row['structure']['id']),row['structure'].get('glytoucan_id')
        if row['composition'].get('glytoucan_id'):
            yield "C"+str(row['composition']['id']),row['composition'].get('glytoucan_id')

    def gettaxa(self,row):
        if row['taxonomy'].get('taxonomy_id'):
          for pub in row['references']:
            pmid = pub.get('pmid')
            yield "S"+str(row['structure']['id']),row['structure'].get('glytoucan_id',"-"),row['taxonomy'].get('taxonomy_id'),pmid
            yield "C"+str(row['composition']['id']),row['composition'].get('glytoucan_id',"-"),row['taxonomy'].get('taxonomy_id'),pmid

    def getpubs(self,row):
        for pub in row['references']:
            if pub.get('pmid'):
                yield "S"+str(row['structure']['id']),row['structure'].get('glytoucan_id',"-"),pub['pmid']
                yield "C"+str(row['composition']['id']),row['structure'].get('glytoucan_id',"-"),pub['pmid']

class UniCarbKBSourceFile(GlyGenSourceFile):
    source = 'unicarbkb'
    glygen_source = "UniCarbKB"
    glygen_sourceid = None
    sections = """ 
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
    """
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
       B2RYF6	10116
    """.splitlines())))
    uckb2gtcacc = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/master/data/uckbcomp2glytoucan.txt"

    def rows(self,filename,loadtaxid=False):
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
        if not hasattr(self,'_comp_lookup'):
            self._comp_lookup = dict()
            for l in self.query_txtfullurl(url=self.uckb2gtcacc).splitlines():
                gtcacc,comp = l.decode('utf8').split()
                self._comp_lookup["comp_"+comp] = gtcacc
        for i,row in enumerate(self.csvrows(filename,filename=filename)):
            row['filename'] = filename
            row['lineno'] = str(i+2);
            # deal with case inconsistency
            for k in list(row.keys()):
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
            if str(row['id']).startswith("comp_") and not row.get('toucan') and row['id'] in self._comp_lookup:
                row['toucan'] = self._comp_lookup[row['id']]
            yield row

    def getgtc(self,row):
        if row.get('toucan') and row.get('id'):
            if re.search('^G\d{5}[A-Z]{2}$',row.get('toucan')):
                yield row['id'],row.get('toucan')

    @uniqueify
    def alltaxa(self):
        for row in self.allrows(loadtaxid=True):
            for acc,gtc,taxid,pmid in self.gettaxa(row):
                yield acc,gtc,taxid,self.glygen_source,acc
                if pmid:
                    yield acc,gtc,taxid,"PubMed",pmid

    def gettaxa(self,row):
        if row.get('taxid') and row.get('id'):
            gtc = row.get('toucan',"")
            if not re.search('^G\d{5}[A-Z]{2}$',gtc):
                gtc = "-"
            if int(row['taxid']) > 0:
                pmid = row.get('pmid')
                yield row['id'],gtc,row['taxid'],pmid

    @uniqueify
    def getpubs(self,row):
        if row.get('id') and row.get('pmid'):
            gtc = row.get('toucan',"")
            if not re.search('^G\d{5}[A-Z]{2}$',gtc):
                gtc = "-"
            if int(row['pmid']) > 0:
                yield row['id'],gtc,row['pmid']

class NCFGSourceFile(GlyGenSourceFile):
    source="ncfg"
    glygen_source="NCFG"
    glygen_sourceid="GLY_000600"
    glygen_pubmedid="30745240"
    sections="glycan_evidence"
    pmidfield="evidence"

    def rows(self,section):
        return self.csvrows(section=self.source,filename=section+"_ncfg")

class SyntheticSourceFile(GlyGenSourceFile):
    source="synthetic"
    glygen_source="GlyGen"
    glygen_sourceid="GLY_000309"
    sections="synthetic"
    pmidfield="evidence"

    def rows(self,section):
        return self.csvrows(section=self.source,path="export_files/current/glycan_synthesized")

class GPTWikiSourceFile(GlyGenSourceFile):
    source="gptwiki"
    glygen_source="GPTwiki"
    glygen_sourceid=None
    sections="glycosites"
    gtcfield="GlyTouCan"
    taxidfixed="9606"

    def rows(self,section):
        return self.csvrows(self.source,filename=section)

class MCWOGlcNAcSourceFile(GlyGenSourceFile):
    source = "mcw_oglcnac"
    glygen_source = "OGlcNAcDB"
    glygen_sourceid = """
       GLY_000517
       GLY_000631
       GLY_000632
       GLY_000633
       GLY_000788
    """
    sections = "human mouse rat fruitfly yeast"
    gtcfixed = "G49108TO"

    def rows(self,section):
        return self.csvrows(section,filename=section+"_o_glcnacome_mcw")

    def gettaxa(self,row):
        if not hasattr(self,'_taxa_lookup'):
            self._taxa_lookup = dict()
            for k,v in self.glygen_species.values():
                self._taxa_lookup[k] = v
        if 'organism' in row:
            org = row['organism']
            assert(org in self._taxa_lookup)
            if not row.get('PMIDS'):
                yield self.gtcfixed,self.gtcfixed,self._taxa_lookup[org],None
            else:
                for pmid in row.get('PMIDS').split(';'):
                    try:
                        pmid = int(pmid)
                    except ValueError:
                        continue
                    if pmid > 0:
                        yield self.gtcfixed,self.gtcfixed,self._taxa_lookup[org],pmid
                    else:
                        yield self.gtcfixed,self.gtcfixed,self._taxa_lookup[org],None

    def getpubs(self,row):
        if row.get('PMIDS'):
            for pmid in row.get('PMIDS').split(';'):
                try:
                    pmid = int(pmid)
                    if pmid > 0:
                        yield self.gtcfixed,self.gtcfixed,pmid,row['section']
                except ValueError:
                    pass

class OGlcNAcAtlasSourceFile(GlyGenSourceFile):
    source = 'atlas_oglcnac'
    glygen_source = "OGlcNAcAtlas"
    sections = "ambiguous_sites unambiguous_sites"
    gtcfixed = "G49108TO"
    taxid2dsid = { 
                  "9606": "GLY_000708",
                  "10090": "GLY_000709",
                  "10116": "GLY_000710",
                  "7227": "GLY_000711",
                  "4932": "GLY_000800",
                 }
    speciesmap = {
                  "Drosophila": "fruitfly",
                 }

    def rows(self,section):
        for row in self.csvrows(self.source,filename=section+'_version_2.0'):
            species = self.speciesmap.get(row['species'],row['species'])
            if species in self.glygen_species:
                row['taxid'] = self.glygen_species[species][1]
            if row.get('taxid') and row['taxid'] in self.taxid2dsid:
                row['_glygen_sourceid'] = self.taxid2dsid[row['taxid']]
            if not row.get('_glygen_sourceid'):
                continue
            if row['pmid']:
                try:
                    int(row['pmid'])
                except ValueError:
                    del row['pmid']
            yield row

class HarvardSourceFile(GlyGenSourceFile):
    source = "harvard"
    glygen_source = "Harvard"
    glygen_sourceid = "GLY_000142"
    glygen_pubmedid = "29351928"
    sections = "T_Cell_O-GlcNAc_sites_CWoo_8-23-18"
    gtcfixed = "G70994MS"
    taxidfixed = "9606"
    pmidfixed = "29351928"
  
    def rows(self,section):
        return self.csvrows(self.source,filename=section)

class OGluCSourceFile(GlyGenSourceFile):
    source = "ogluc"
    glygen_source = "GlyGen"
    glygen_sourceid = "GLY_000716"
    glygen_pubmedid = "34411563"
    sections = "human"
    gtcfield = "saccharide"
    pmidfield = "evidence"
    
    def rows(self,section):
        return self.csvrows(section,path="export_files/current/%s_proteoform_glycosylation_sites_o-gluc"%(section,))

class LiteratureSourceFile(GlyGenSourceFile):
    source = "literature"
    glygen_source = "GlyGen"
    glygen_sourceid = "GLY_000335"
    glygen_pubmedid = "18187336"
    sections = "hcv1a"
    # glygen_sourceid = "GLY_000612 GLY_000335"
    # glygen_pubmedid = "16442106 18187336"
    # sections = "sarscov1 hcv1a"
    gtcfield = "saccharide"
    taxidfield = "tax_id_uniprotkb_ac"
    pmidfield = "evidence"
    
    def rows(self,section):
        return self.csvrows(section,path="export_files/current/%s_proteoform_glycosylation_sites_literature"%(section,))

class BiomarkerDBSourceFile(GlyGenSourceFile):
    source = "biomarkerdb"
    glygen_source = "BiomarkerKB"
    glygen_sourceid = "GLY_000737"
    sections = "allbiomarkers-all"
    gtcfield = "GlyTouCan"
    taxidfixed = "9606"

    def rows(self,section):
        for row in self.csvrows(section=self.source,filename=section):
            for xref in map(str.strip,row.get('Main x-ref',"").split("|")):
                if xref.startswith("GTC:"):
                    row['GlyTouCan'] = xref[4:]
                    yield row

class EMBLSourceFile(GlyGenSourceFile):
    source = "embl"
    glygen_source = "EMBL"
    sections = "glygen_upload"
    taxid2dsid = { 
                  "9606": "GLY_000888",
                  "10090": "GLY_000889",
                 }
    byonic2gtcacc = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/master/data/byonic2glytoucan.txt"

    def rows(self,section):
        if not hasattr(self,'_comp_lookup'):
            self._comp_lookup = dict()
            for l in self.query_txtfullurl(url=self.byonic2gtcacc).splitlines():
                gtcacc,byonic = l.decode('utf8').split()
                self._comp_lookup[byonic] = gtcacc
        for row in self.csvrows(self.source,filename=section):
            comp = row.get('composition').split(' % ')[0].strip()
            taxid = row.get('taxonomy_id')
            if not comp or not taxid:
                continue
            if comp not in self._comp_lookup:
                continue
            row['glytoucan_ac'] = self._comp_lookup[comp]
            row['taxid'] = taxid
            row['_glygen_sourceid'] = self.taxid2dsid[taxid]
            yield row

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
    """
    #     yeast_proteoform_glycosylation_sites_oglcnac_mcw
    #     GLY_000142      human_proteoform_glycosylation_sites_harvard                9606
    #     GLY_000335      hcv1a_proteoform_glycosylation_sites_literature             11108
    #     GLY_000632      rat_proteoform_glycosylation_sites_oglcnac_mcw              10116
    #     GLY_000633      fruitfly_proteoform_glycosylation_sites_oglcnac_mcw         7227
    #     GLY_000631      mouse_proteoform_glycosylation_sites_oglcnac_mcw            10090
    #     GLY_000517      human_proteoform_glycosylation_sites_oglcnac_mcw            9606
    #     GLY_000708      human_proteoform_glycosylation_sites_oglcnac_atlas          9606
    #     GLY_000709      mouse_proteoform_glycosylation_sites_oglcnac_atlas          10090
    #     GLY_000710      rat_proteoform_glycosylation_sites_oglcnac_atlas            10116
    #     GLY_000711      fruitfly_proteoform_glycosylation_sites_oglcnac_atlas       7227
    #     GLY_000716      human_proteoform_glycosylation_sites_o_gluc                 9606

    @uniqueify
    def dataset_allgtc(self,dsid):
        for row in self.dataset(dsid):
            yield row['accession'],row['accession']
            
    @uniqueify
    def dataset_alltaxa(self,dsid):
        for row in self.dataset(dsid):
            yield row['accession'],row['accession'],row['taxid'],None
            
    @uniqueify
    def dataset_allpubs(self,dsid):
        for row in self.dataset(dsid):
            if 'pmid' in row:
                for pmid in row['pmid'].split(';'):
                    yield row['accession'],row['accession'],pmid
            
    @uniqueify
    def allgtc(self):
        for dsid in self.datasets:
            for row in self.dataset_allgtc(dsid):
                yield row[0],row[1],dsid

    @uniqueify
    def alltaxa(self):
        for dsid in self.datasets:
            for row in self.dataset_alltaxa(dsid):
                yield row[0],row[1],row[2],dsid

    @uniqueify
    def allpubs(self):
        for dsid in self.datasets:
            for row in self.dataset_allpubs(dsid):
                yield row[0],row[1],row[2],dsid

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
