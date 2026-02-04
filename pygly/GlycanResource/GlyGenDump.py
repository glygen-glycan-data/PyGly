
from .WebServiceResource import WebServiceResource
import os, os.path, re, traceback, sys
from collections import defaultdict

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
        "pig": ('Sus scrofa', '9823'),
        "chicken": ('Gallus gallus', '9031'),
        "arabidopsis": ('Arabidopsis thaliana', '3702'),
        "bovine": ('Bos taurus', '9913'),
        "zebrafish": ('Danio rerio', '7955'),
        "hamster": ('Cricetulus griseus', '10029'),
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
        anyrows = False
        for row in self.json2rows(self.query_jsonsourcefile(source=self.source,filename=filename)):
            for k,v in list(row.items()):
                if isinstance(v,str):
                    row[k] = v.strip()
                if not row[k]:
                    del row[k]
            anyrows = True
            yield row
        if not anyrows:
            raise RuntimeError("No rows from JSON source %s:%s"%(self.source,filename))

    def csvrows(self,section,filename=None,path=None,extn="csv"):
        if filename:
            path = self.source + "/current/" + filename
        if extn is None:
            extn = ""
        else:
            extn = "."+extn
        anyrows = False
        for row in self.query_csvfullurl(url=self.apiurl+'/'+path+extn):
            for k,v in list(row.items()):
                if isinstance(v,str):
                    row[k] = v.strip()
                if not row[k]:
                    del row[k]
            anyrows = True
            yield row
        if not anyrows:
            raise RuntimeError("No rows from CSV source %s:%s"%(self.source,filename))

    def tsvrows(self,section,filename=None,path=None,extn="tsv"):
        if filename:
            path = self.source + "/current/" + filename
        anyrows = False
        for row in self.query_tsvfullurl(url=self.apiurl+'/'+path+"."+extn):
            for k,v in list(row.items()):
                if isinstance(v,str):
                    row[k] = v.strip()
                if not row[k]:
                    del row[k]
            anyrows = True
            yield row
        if not anyrows:
            raise RuntimeError("No rows from TSV source %s:%s"%(self.source,filename))

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
        if not hasattr(self,'glygen_doi'):
            ggdois = [None]*len(sections)
        else:
            ggdois = self.glygen_doi.split()
        assert(len(sections) == len(ggpmids))
        for sec,ggsid,ggpmid,ggdoi in zip(sections,ggsids,ggpmids,ggdois):
            for row in self.rows(sec,**kw):
                if ggsid:
                    row['_glygen_sourceid'] = ggsid
                if ggpmid:
                    row['_glygen_pubmedid'] = ggpmid
                if ggdoi:
                    row['_glygen_doi'] = ggdoi
                yield row

    idfield = 'glytoucan_ac'
    idfixed = None
    gtcfield = 'glytoucan_ac'
    gtcfixed = None
    taxidfield = 'taxid'
    taxidfixed = None
    tissuefield = None
    tissuefixed = None
    pmidfield = 'pmid'
    pmidfixed = None
    doifield = None
    doifixed = None

    def getgtc(self,row):
        if row.get(self.gtcfield,self.gtcfixed):
            yield row.get(self.idfield,self.idfixed),row.get(self.gtcfield,self.gtcfixed)

    def gettaxa(self,row):
        if row.get(self.idfield,self.idfixed) and row.get(self.taxidfield,self.taxidfixed):
            yield row.get(self.idfield,self.idfixed),row.get(self.gtcfield,self.gtcfixed),row.get(self.taxidfield,self.taxidfixed),row.get(self.pmidfield,self.pmidfixed),row.get(self.doifield,self.doifixed)

    def gettissue(self,row):
        if row.get(self.idfield,self.idfixed) and row.get(self.taxidfield,self.taxidfixed) and row.get(self.tissuefield,self.tissuefixed):
            yield row.get(self.idfield,self.idfixed),row.get(self.gtcfield,self.gtcfixed),row.get(self.taxidfield,self.taxidfixed),row.get(self.tissuefield,self.tissuefixed),row.get(self.pmidfield,self.pmidfixed),row.get(self.doifield,self.doifixed)

    def getpubs(self,row):
        if row.get(self.idfield,self.idfixed) and (row.get(self.pmidfield,self.pmidfixed) or row.get(self.doifield,self.doifixed)):
            yield row.get(self.idfield,self.idfixed),row.get(self.gtcfield,self.gtcfixed),row.get(self.pmidfield,self.pmidfixed),row.get(self.doifield,self.doifixed)

    @uniqueify
    def allgtc(self):
        for row in self.allrows():
            for acc,gtc in self.getgtc(row):
                yield acc,gtc,self.glygen_source,row.get('_glygen_sourceid')

    @uniqueify
    def alltaxa(self):
        for row in self.allrows():
            for acc,gtc,taxid,pmid,doi in self.gettaxa(row):
                yield acc,gtc,taxid,self.glygen_source,row.get('_glygen_sourceid')
                if pmid:
                    yield acc,gtc,taxid,"PubMed",pmid
                if doi:
                    yield acc,gtc,taxid,"DOI",doi
                if row.get('_glygen_pubmedid'):
                    yield acc,gtc,taxid,"PubMed",row['_glygen_pubmedid']
                if row.get('_glygen_doi'):
                    yield acc,gtc,taxid,"DOI",row['_glygen_doi']

    @uniqueify
    def alltissue(self):
        for row in self.allrows():
            for acc,gtc,taxid,tissue,pmid,doi in self.gettissue(row):
                yield acc,gtc,taxid,tissue,self.glygen_source,row.get('_glygen_sourceid')
                if pmid:
                    yield acc,gtc,taxid,tissue,"PubMed",pmid
                if doi:
                    yield acc,gtc,taxid,tissue,"DOI",doi
                if row.get('_glygen_pubmedid'):
                    yield acc,gtc,taxid,tissue,"PubMed",row['_glygen_pubmedid']
                if row.get('_glygen_doi'):
                    yield acc,gtc,taxid,tissue,"DOI",row['_glygen_doi']

    @uniqueify
    def allpubs(self):
        for row in self.allrows():
            for acc,gtc,pmid,doi in self.getpubs(row):
                if pmid:
                    yield acc,gtc,"PubMed",pmid,self.glygen_source,row.get('_glygen_sourceid')
                if doi:
                    yield acc,gtc,"DOI",doi,self.glygen_source,row.get('_glygen_sourceid')

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

    def allsourcetissue(self):
        for src in self.sources:
            ggsf = eval(src+"SourceFile")(verbose=self._verbose)
            for row in ggsf.alltissue():
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
        MWMetRef
        Synthetic
        GPTWiki
        MCWOGlcNAc
        PlateletOLinked
        TableMaker
        OGlcNAcAtlas
        Harvard
        OGluC
        Literature
        Diabetes
        TwinsUK
        BiomarkerDB
        EMBL
        GlyProtDB
        GlycomeAtlas
        MatrixDB
        PDCCCRCC
        PDBGlycan
    """.split()

class GlyConnectSourceFile(GlyGenSourceFile):
    source = 'glyconnect'
    glygen_source = "GlyConnect"
    glygen_sourceid = None
    sections = 'human mouse rat fruitfly sarscov2 yeast dicty pig chicken arabidopsis bovine hamster'

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
            yield "S"+str(row['structure']['id']),row['structure'].get('glytoucan_id',"-"),row['taxonomy'].get('taxonomy_id'),pmid,None
            yield "C"+str(row['composition']['id']),row['composition'].get('glytoucan_id',"-"),row['taxonomy'].get('taxonomy_id'),pmid,None

    def getpubs(self,row):
        for pub in row['references']:
            if pub.get('pmid'):
                yield "S"+str(row['structure']['id']),row['structure'].get('glytoucan_id',"-"),pub['pmid'],None
                yield "C"+str(row['composition']['id']),row['structure'].get('glytoucan_id',"-"),pub['pmid'],None

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
    doifield = "doi"

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
            for l in self.query_txtfullurl(url=self.uckb2gtcacc):
                gtcacc,comp = l.split()
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
            for acc,gtc,taxid,pmid,doi in self.gettaxa(row):
                yield acc,gtc,taxid,self.glygen_source,acc
                if pmid:
                    yield acc,gtc,taxid,"PubMed",pmid
                if doi:
                    yield acc,gtc,taxid,"DOI",doi

    def gettaxa(self,row):
        if row.get('taxid') and row.get('id'):
            gtc = row.get('toucan',"")
            if not re.search('^G\d{5}[A-Z]{2}$',gtc):
                gtc = "-"
            if int(row['taxid']) > 0:
                pmid = row.get('pmid')
                yield row['id'],gtc,row['taxid'],pmid,None

    @uniqueify
    def getpubs(self,row):
        if row.get('id') and row.get('pmid'):
            gtc = row.get('toucan',"")
            if not re.search('^G\d{5}[A-Z]{2}$',gtc):
                gtc = "-"
            if int(row['pmid']) > 0:
                yield row['id'],gtc,row['pmid'],None

class NCFGSourceFile(GlyGenSourceFile):
    source="ncfg"
    glygen_source="NCFG"
    glygen_sourceid="GLY_000600"
    glygen_pubmedid="30745240"
    sections="glycan_evidence"
    pmidfield="evidence"

    def rows(self,section):
        return self.csvrows(section=self.source,filename=section+"_ncfg")

class MWMetRefSourceFile(GlyGenSourceFile):
    source="mw"
    glygen_source="MetabolomicsWorkbench"
    glygen_sourceid="GLY_001169"
    sections="mw_refmet_mapping_result"
    idfield="regno"
    gtcfield="glytoucan_ac"

    def rows(self,section):
        return self.tsvrows(section=self.source,filename=section)

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
    # glygen_sourceid="GLY_000480"
    sections="glycosites"
    idfield="GlyTouCan"
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
       GLY_000959
       GLY_001049
       GLY_001171
       GLY_001250
    """
    sections = "human mouse rat fruitfly yeast chicken arabidopsis bovine zebrafish"
    gtcfixed = "G49108TO"
    idfixed = "G49108TO"

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
                yield self.gtcfixed,self.gtcfixed,self._taxa_lookup[org],None,None
            else:
                for pmid in row.get('PMIDS').split(';'):
                    try:
                        pmid = int(pmid)
                    except ValueError:
                        continue
                    if pmid > 0:
                        yield self.gtcfixed,self.gtcfixed,self._taxa_lookup[org],pmid,None
                    else:
                        yield self.gtcfixed,self.gtcfixed,self._taxa_lookup[org],None,None

    def getpubs(self,row):
        if row.get('PMIDS'):
            for pmid in row.get('PMIDS').split(';'):
                try:
                    pmid = int(pmid)
                    if pmid > 0:
                        yield self.gtcfixed,self.gtcfixed,pmid,None
                except ValueError:
                    pass

class PlateletOLinkedSourceFile(GlyGenSourceFile):
    source = "user_submission/platelet_o_linked"
    glygen_source = "UserSubmission"
    glygen_sourceid = """
        GLY_001051
        GLY_001051
    """
    sections = """
        mmc5
        mmc7
    """
    taxidfixed = '9606'
    pmidfixed = '38237698'

    byonic2gtcacc = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/master/data/byonic2glytoucan.txt"

    def rows(self,section):
        if not hasattr(self,'_comp_lookup'):
            self._comp_lookup = dict()
            for l in self.query_txtfullurl(url=self.byonic2gtcacc):
                gtcacc,byonic = l.split()
                self._comp_lookup[byonic] = gtcacc
        for row in self.tsvrows(section=section,filename='1-s2.0-S1535947624000070-'+section):
            for glycomp in map(str.strip,row['Glycans NHFAGNa'].split(',')):
                if glycomp in self._comp_lookup:
                    row['glytoucan_ac'] = self._comp_lookup[glycomp]
                    yield row

class TableMakerSourceFile(GlyGenSourceFile):
    source = "tablemaker"
    glygen_source = "TableMaker"
    taxid2taxid = {
        "9605": "9606"
    }
    taxid2dsid = { 
        "9606": "GLY_001408",  #Human
        "10090": "GLY_001485", #Mouse
        "10116": "GLY_001486", #Rat
        "7227": "GLY_001253",  #FruitFly
        "4932": "GLY_001494",  #Yeast
        "9823": "GLY_001490",  #Pig
        "9031": "GLY_001489",  #Chicken
        "3702": "GLY_001487",  #Arabidopsis
        "9913": "GLY_001497",  #Bovine
        "7955": "GLY_001252",  #Zebrafish
        "10029": "GLY_001484", #Hamster
        "44689": "GLY_001488", #Dicty
        }
    sections = """
        TG1468274
	    TG1127346
	    TG1040297
	    TG2415036
	    TG2986138
	    TG4951079
	    TG7181079
	    TG7732158
	    TG8189974
	    TG8741863
        TG4253352
        TG5664582
        TG2727671
        TP6640575
        TG4964642
        TG7680256
        TG1131507
        TP5868669
        TG10381046
        TG5577100
    """
    taxidfield = 'Species'
    # pmidfield = 'Evidence' sometimes contains a doi!
    doifield = "doi"
    tissuefield = 'Tissue'
    idfield = 'GlyTouCan ID'
    gtcfield = 'GlyTouCan ID'    

    def rows(self,section):
        for row in self.csvrows(section=section,filename=section):
            # print(row)
            if row.get('Evidence',"").strip():
                if re.search(r'^\d+$',row['Evidence'].strip()):
                    row['pmid'] = row['Evidence'].strip()
                elif row['Evidence'].strip().startswith('10.'):
                    row['doi'] = row['Evidence'].strip()
                else:
                    raise RuntimeError(str(row))
            if row.get(self.taxidfield,"").strip():
                taxid = row[self.taxidfield]
                taxid = self.taxid2taxid.get(taxid,taxid)
                if taxid in self.taxid2dsid:
                    row['_glygen_sourceid'] = self.taxid2dsid[taxid]
                    yield row

class OGlcNAcAtlasSourceFile(GlyGenSourceFile):
    source = 'atlas_oglcnac'
    glygen_source = "OGlcNAcAtlas"
    sections = "ambiguous_sites unambiguous_sites"
    gtcfixed = "G49108TO"
    idfixed = "G49108TO"
    taxid2dsid = { 
                  "9606": "GLY_000708",  #Human
                  "10090": "GLY_000709", #Mouse
                  "10116": "GLY_000710", #Rat
                  "7227": "GLY_000711",  #FruitFly
                  "4932": "GLY_000800",  #Yeast
                  "9823": "GLY_000955",  #Pig
                  "9031": "GLY_000962",  #Chicken
                  "3702": "GLY_001050",  #Arabidopsis
                  "9913": "GLY_001172",  #Bovine
                 }
    speciesmap = {
                  "Drosophila": "fruitfly",
                  "Arabidopsis": "arabidopsis",
                 }

    def rows(self,section):
        for row in self.csvrows(self.source,filename=section):
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
    idfield = "saccharide"
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
    idfield = "saccharide"
    gtcfield = "saccharide"
    taxidfield = "tax_id_uniprotkb_ac"
    pmidfield = "evidence"
    
    def rows(self,section):
        return self.csvrows(section,path="export_files/current/%s_proteoform_glycosylation_sites_literature"%(section,))

class DiabetesSourceFile(GlyGenSourceFile):
    glygen_source = "GornikLab"
    source = "export_files"
    sections = "none"
    idfield = "ID"
    gtcfield = "GlyTouCan"
    taxidfixed = "9606"
    glygen_sourceid = "GLY_000960"
    glygen_pubmedid = "28905229"
    pmidfixed = "28905229"

    def rows(self,section):
        return self.csvrows("export_files",filename="diabetes_glycomic_mapping")

class TwinsUKSourceFile(GlyGenSourceFile):
    glygen_source = "TwinsUK"
    source = "twinsuk"
    sections = "twinsuk_plasma_nglycan_mapping"
    idfield = "GlyToucan ID"
    gtcfield = "GlyToucan ID"
    taxidfixed = "9606"
    glygen_sourceid = "GLY_001461"
    pmidfield = "PMID"
    tissuefield = "Tissue ID"

    def rows(self,section):
        for row in self.tsvrows(section=self.source,filename=section):
            yield row

class BiomarkerDBSourceFile(GlyGenSourceFile):
    source = "biomarkerdb"
    glygen_source = "BiomarkerKB"
    glygen_sourceid = "GLY_000737"
    sections = "oncomx"
    idfield = "GlyTouCan"
    gtcfield = "GlyTouCan"
    taxidfixed = "9606"

    def rows(self,section):
        for row in self.tsvrows(section=self.source,filename=section):
            # print(row)
            pmid = row.get("evidence_source","")
            if pmid.lower().startswith("pubmed:"):
                row['pmid'] = pmid.split(':',1)[1]
            xref = row.get("assessed_biomarker_entity_id","")
            if xref.startswith("GTC:"):
                row['GlyTouCan'] = xref.split(':',1)[1].split()[0]
                yield row

class EMBLSourceFile(GlyGenSourceFile):
    source = "embl"
    glygen_source = "EMBL"
    glygen_doi = "10.1101/2023.09.13.557529"
    sections = "glygen_upload"
    taxid2dsid = { 
                  "9606": "GLY_000888",
                  "10090": "GLY_000889",
                 }
    byonic2gtcacc = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/master/data/byonic2glytoucan.txt"

    def rows(self,section):
        if not hasattr(self,'_comp_lookup'):
            self._comp_lookup = dict()
            for l in self.query_txtfullurl(url=self.byonic2gtcacc):
                gtcacc,byonic = l.split()
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

class GlyProtDBSourceFile(GlyGenSourceFile):
    source = "glycosmos"
    glygen_source = "GlyCosmos"
    glygen_sourceid = "GLY_001251"
    sections = "filtered_glycoproteins_glycan_list"

    def rows(self,section):
        if self._taxidlookup == None:
            self._taxidlookup = {}
            for key,(sp,taxid) in self.glygen_species.items():
                self._taxidlookup[sp] = taxid
        for row in self.csvrows(self.source,filename=section):
            for gtcacc in filter(str.strip,row['GlyTouCan IDs'].split(",")):
                yield dict(glytoucan_ac=gtcacc,taxid=self._taxidlookup[row['Organism']])

class GlycomeAtlasSourceFile(GlyGenSourceFile):
    source = "glycosmos"
    glygen_source = "GlyCosmos"
    taxid2dsid = { 
                  "9606":  "GLY_001409",  #Human
                  "10090": "GLY_001411", #Mouse
                  "7955":  "GLY_001410",  #Zebrafish
    }
    sections = "glycosmos_glycomeatlas_list"
    idfield = "GlyTouCan ID"
    gtcfield = "GlyTouCan ID"
    taxidfield = "Taxonomy ID"
    tissuefield = "Tissue"
    pmidfield = "PMID"
    doifield = "DOI"
    
    def rows(self,section):
        for row in self.csvrows(self.source,filename=section):
            if row.get(self.taxidfield) and row[self.taxidfield] in self.taxid2dsid:
                row['_glygen_sourceid'] = self.taxid2dsid[row[self.taxidfield]]
            if not row.get('_glygen_sourceid'):
                continue
            if row.get(self.tissuefield):
                for ts in list(row[self.tissuefield].split(',')):
                    ts = ts.strip().split('(',1)[1].rstrip(')').replace('_',':')
                    row1 = dict(row.items())
                    row1[self.tissuefield] = ts
                    yield row1
            else:
                yield row
            

class MatrixDBSourceFile(GlyGenSourceFile):
    source = "matrixdb"
    glygen_source = "MatrixDB"
    sections = "matrixdb_CORE"
    gagurl = "https://matrixdb.univ-lyon1.fr/download/Custom_MatrixDB_biomolecules.tsv"
    gagurl = "https://matrixdb.univ-lyon1.fr/downloads/biomolecules.csv"
    idfield = "id"
    gtcfield = "gtcacc"

    def rows(self,section):
        if not hasattr(self,'_gag_lookup'):
            self._gag_lookup = dict()
            header = None
            for l in self.query_txtfullurl(url=self.gagurl):
                sl = l.split('\t')
                if not header:
                    header = sl
                    continue
                row = dict(zip(header,sl))
                if not row['MatrixDB identifier'].startswith('GAG_'):
                    continue
                gtcid = None
                for exid in row['External identifiers'].split(','):
                    if exid.startswith('GlyTouCan:'):
                        gtcid = exid.split(':',1)[1]
                        break
                if not gtcid:
                    continue
                self._gag_lookup[row['MatrixDB identifier'].strip()] = gtcid
        gagseen = set()
        for row in self.tsvrows(self.source,filename=section,extn='tab'):
            for key in ("#ID(s) interactor A","ID(s) interactor A","ID(s) interactor B","Alt. ID(s) interactor A","Alt. ID(s) interactor B"):
                if key in row and row[key].startswith('matrixdb:GAG_'):
                    gagid = row[key].split(':',1)[1]
                    if gagid not in gagseen and gagid in self._gag_lookup:
                        gagseen.add(gagid)
                        yield dict(id=gagid,gtcacc=self._gag_lookup[gagid])

class PDCCCRCCSourceFile(GlyGenSourceFile):
    source = "pdc"
    glygen_source = "PDC-CCRCC"
    sections = "none"
    taxidfixed = "9606"
    glygen_sourceid = "GLY_000961"
    glygen_pubmedid = "37074911"
    pmidfixed = "37074911"
    shortcomp2gtcacc = "https://raw.githubusercontent.com/glygen-glycan-data/GNOme/master/data/shortcomp2glytoucan.txt"

    def rows(self,section):
        if not hasattr(self,'_comp_lookup'):
            self._comp_lookup = dict()
            for l in self.query_txtfullurl(url=self.shortcomp2gtcacc):
                gtcacc,shcomp = l.split()
                self._comp_lookup[shcomp] = gtcacc
        for row in self.tsvrows(self.source,filename="ccRCC_TMT_intact_glycopeptide_abundance_MD-MAD"):
            # print(row['Nglycan'])
            splshcomp = re.split(r'([A-Z])',row['Nglycan'])
            shcompdict = defaultdict(int)
            for i in range(1,len(splshcomp),2):
                if int(splshcomp[i+1]) > 0:
                    shcompdict[splshcomp[i]] = int(splshcomp[i+1])
            normshcomp = ""
            for m in "HNFS":
                if shcompdict[m] > 1:
                    normshcomp += (m + str(shcompdict[m]))
                elif shcompdict[m] == 1:
                    normshcomp += m
            # print(normshcomp)
            if normshcomp not in self._comp_lookup:
                continue
            gtcacc = self._comp_lookup[normshcomp]
            row['glytoucan_ac'] = gtcacc
            yield row

class PDBGlycanSourceFile(GlyGenSourceFile):
    source = "pdb"
    glygen_source = "PDB"
    sections = "glycosites_rcsb_pdb"
    taxidfixed = "9606"
    glygen_sourceid = "GLY_001412"
    idfield = "src_xref_id"
    doifield = "doi"
    gtcfield = "saccharide"
    
    def rows(self,section):
        for row in self.csvrows(self.source,filename=section):
            if row.get('xref_id'):
                if re.search(r'^\d+$',row['xref_id'].strip()):
                    row['pmid'] = row['xref_id'].strip()
                elif row['xref_id'].startswith('10.'):
                    row['doi'] = row['xref_id'].strip()
            # print(row)
            yield row

class GlycoShapeSourceFile(GlyGenSourceFile):
    source = "glycoshape"
    glygen_source = "GlycoShape"
    sections = "glycan_list"
    glygen_sourceid = "GLY_001491"
   
    def json2rows(self,json):
        for acc in json:
            yield dict(glytoucan_ac=acc)

    def rows(self,section):
        return self.jsonrows(self.source,filename=section)
        
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
