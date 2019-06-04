#!/bin/env python27

import sys
import traceback
import urllib2
import time
import csv
import re
from pygly.memoize import memoize
from rdflib import ConjunctiveGraph, Namespace


class APIError(RuntimeError):
    pass

class UnicarbURLNotSupported(APIError):
    pass

class abstract_api(object):

    def __init__(self):
        self._lastrequesttime = 0
        self._lastrequestcount = 0
        self.delaytime = .2
        self.delaybatch = 1

    def _wait(self, delaytime=None):
        if delaytime != None:
            time.sleep(delaytime)
            return
        if (self._lastrequestcount % self.delaybatch) == 0 and self._lastrequestcount > 0:
            time.sleep(self.delaytime)
        self._lastrequesttime = time.time()
        self._lastrequestcount += 1


class Triple_store_api(abstract_api):
    endpt = ''

    def __init__(self):
        super(Triple_store_api, self).__init__()
        self.g = None
        self.ssg = None
        self.opener = None
        self.maxattempts = 3
        self.alphamap = None

    def setup_sparql(self):
        self.g = ConjunctiveGraph(store='SPARQLStore')
        self.g.open(self.endpt)

    @memoize()
    def query(self, sparql):
        self._wait()

        if self.g == None:
            self.setup_sparql()

        attempt = 0
        response = None
        while response == None and attempt < self.maxattempts:
            try:
                attempt += 1
                response = self.g.query(sparql)
            except:
                traceback.print_exc()
                self._wait(self.delaytime ** attempt)

        if response == None:
            raise IOError("Cannot query SPARQL endpoint")

        return response


class web_api(abstract_api):
    api = ''
    user = ''
    apikey = ''

    def __init__(self, user=None, apikey=None):
        super(web_api, self).__init__()
        self.user = user
        self.apikey = apikey

    def setup_api(self, user=None, apikey=None):
        if user == None:
            user = self.user
            apikey = self.apikey
        if user == None:
            user, apikey = self.getcredentials()
        # print user,apikey
        self.password_mgr = urllib2.HTTPPasswordMgrWithDefaultRealm()
        self.password_mgr.add_password(None, self.api, user, apikey)
        self.handler = urllib2.HTTPBasicAuthHandler(self.password_mgr)
        self.opener = urllib2.build_opener(self.handler)

    def getcredentials(self):
        raise NotImplemented


class UniCarb(Triple_store_api):
    # endpt = 'http://137.92.56.159:40935/unicarbkb'
    # endpt = 'http://137.92.56.159:40935//unicarbkb'
    # endpt = 'http://137.92.56.159:40935/unicarbkbv2_test_only/query'
    # endpt = 'http://137.92.56.159:40935/unicarbkb2.0/query'
    # endpt = 'http://ts.glytoucan.org/sparql'
    # endpt = 'http://203.101.226.16:40935/unicarbkbv2.0.1'
    # endpt = 'http://203.101.226.16:40935/unicarbkbv2.0.1/query'
    endpt = 'http://203.101.226.16:40935/unicarbkbv2.0.4/query'

    unicarb2glytoucan_query = """
    PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>

    SELECT DISTINCT ?structure ?glytoucan
    WHERE {
      ?structure a glycan:Saccharide .
      ?structure glycan:has_glytoucan_id ?glytoucan
    }
    """

    taxonomy_query = """
    PREFIX glycan: <http://purl.jp/bio/12/glyco/glycan#>
    PREFIX glyconj: <http://purl.jp/bio/12/glyco/conjugate#>
    
    SELECT DISTINCT ?structure ?taxon
    WHERE {
        ?structure a glycan:Saccharide .
    
        ?refsaccharide a glyconj:ReferencedSaccharide .
    
        ?refsaccharide glycan:has_glycan ?structure .
        ?refcompound a glyconj:ReferencedGlycoconjugate .
    
        ?refcompound glyconj:ReferencedSaccharide ?refsaccharide .
        ?refcompound glycan:is_from_source ?source .
        ?source a glycan:source_natural .
    
        ?source glycan:has_taxon ?taxon
    }
    """

    def __init__(self, useCache=True):
        pass
    
    def stripURL(self, url):
        return url.split("/")[-1]
    
    def intConversion(self, s):
        temp = re.compile(r"^\d+\^\^xsd\:int$").findall(s)
        if len(temp) == 1:
            t = temp[0]
            return int(re.compile(r"^\d+").findall(s)[0])
        else:
            raise RuntimeError
    
    def unicarbGlycanURL2id(self, url):
        glycanType = url.split("/")[-2]
        if glycanType == "structure":
            key = url.split("/")[-1]
        elif glycanType == "composition":
            key = url.split("/")[-1]
        else:
            print "A glycan is neither structure nor composition"
            raise ValueError
        return glycanType, key
    
    def getProtein(self, filepath):
        # The columns will be: UniProt ID,Position,PubMed ID,UniCarb ID,GlyTouCan ID,AminoAcid,TypeAminoAcid,Notes
        res = []
        with open(filepath) as f:
            for i, row in enumerate(csv.reader(f)):
                if i == 0:
                    continue
                processedRow = []
                row[0] = self.stripURL(row[0])
                row[1] = self.intConversion(row[1])
                res.append(row)
        return res
    
    def protein_human(self):
        return self.getProtein("./unicarb_temp/data_files_unicarbkb_QC_build_triplestore_2.0.4_human_unicarbkb_2_0_4_human.csv")
    
    def protein_mouse(self):
        return self.getProtein("./unicarb_temp/data_files_unicarbkb_QC_build_triplestore_2.0.4_mouse_unicarbkb_2_0_4_mouse.csv")
    
    def getAllProtein(self):
        return self.protein_human() + self.protein_mouse()
    
    def unicarb2glytoucan(self):
        g = ConjunctiveGraph(store='SPARQLStore')
        g.open(self.endpt)
        results = g.query(self.unicarb2glytoucan_query)
        res = {}
        for r in results.bindings:
            row = map(str, map(r.get, results.vars))
            if self.unicarbGlycanURL2id(row[0])[0] == 'structure':
                uid = self.unicarbGlycanURL2id(row[0])[1]
            else:
                continue
            gid = row[1]

            res[uid] = gid
        return res

    def taxonomy(self):
        g = ConjunctiveGraph(store='SPARQLStore')
        g.open(self.endpt)
        results = g.query(self.taxonomy_query)
        res = {}
        for r in results.bindings:
            row = map(str, map(r.get, results.vars))
            if self.unicarbGlycanURL2id(row[0])[0] == 'structure':
                uid = self.unicarbGlycanURL2id(row[0])[1]
            else:
                continue
            taxonomy = self.stripURL(row[1])
            
            if uid not in res:
                res[uid] = []
            
            res[uid].append(taxonomy)
        return res
    
    idmapping = False
    def IDmapping(self, UniCarbID):
        if not self.idmapping:
            self.idmapping = self.unicarb2glytoucan()
        UniCarbID = str(UniCarbID)
        return self.idmapping[UniCarbID]



if __name__ == "__main__":
    uc = UniCarb()
    #uc.q1()
    uc.IDmapping("6156")

    uc.protein_mouse()
    
    uc.protein_human()
    
    uc.taxonomy()
