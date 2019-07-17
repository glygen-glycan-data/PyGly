#!/bin/env python27

import sys
import traceback
import urllib2
import time
import csv
import re
import os
import gzip
import atexit
import cPickle as pickle
from memoize import memoize
from rdflib import ConjunctiveGraph, Namespace
from collections import defaultdict


class APIError(RuntimeError):
    pass


class UnicarbURLNotSupported(APIError):
    pass


class abstract_api(object):
    cachefile = ".abstractcache"

    def __init__(self, usecache=False):
        self._lastrequesttime = 0
        self._lastrequestcount = 0
        self.delaytime = .2
        self.delaybatch = 1

        self.usecache = usecache
        self.cacheupdated = False
        if usecache:
            self.cache_init()

    def cache_init(self):
        self.cachedata = None
        self.cacheupdated = False

        if os.path.exists(self.cachefile):
            try:
                self.cachedata = pickle.load(gzip.open(self.cachefile, 'rb'))
                # print >>sys.stderr, "Loaded cached data"
            except:
                self.cachedata = {}
        else:
            self.cachedata = {}
        atexit.register(self.__del__)

    def __del__(self):
        if self.cacheupdated:
            try:
                wh = gzip.open(self.cachefile, 'wb')
                pickle.dump(self.cachedata, wh, -1)
                wh.close()
                # print >>sys.stderr, "Saved cached data"
                self.cacheupdated = False
            except:
                pass

    def cachegetmany(self, valuekey, acc, iterable):
        if valuekey not in self.cachedata:
            self.cachedata[valuekey] = reduce(lambda d, x: d.setdefault(x[0], []).append(x[1]) or d,
                                              iterable, {})
            self.cacheupdated = True
        return self.cachedata[valuekey].get(acc, [])

    def cacheget(self, valuekey, acc, iterable):
        if valuekey not in self.cachedata:
            self.cachedata[valuekey] = dict(iterable)
            self.cacheupdated = True
        return self.cachedata[valuekey].get(acc)

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

    def __init__(self,usecache=False):
        super(Triple_store_api, self).__init__(usecache=usecache)
        self.g = None
        self.ssg = None
        self.opener = None
        self.maxattempts = 3
        self.alphamap = None

        self.setup_sparql()

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


class UniCarbKB(Triple_store_api):
    # endpt = 'http://137.92.56.159:40935/unicarbkb'
    # endpt = 'http://137.92.56.159:40935//unicarbkb'
    # endpt = 'http://137.92.56.159:40935/unicarbkbv2_test_only/query'
    # endpt = 'http://137.92.56.159:40935/unicarbkb2.0/query'
    # endpt = 'http://ts.glytoucan.org/sparql'
    # endpt = 'http://203.101.226.16:40935/unicarbkbv2.0.1'
    # endpt = 'http://203.101.226.16:40935/unicarbkbv2.0.1/query'
    endpt = 'http://203.101.226.16:40935/unicarbkbv2.0.4/query'

    cachefile = ".unicarbcache"

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

    protein_query = """
    prefix sio: <http://semanticscience.org/resource/>
    prefix owl:   <http://www.w3.org/2002/07/owl#>
    prefix gco: <http://purl.jp/bio/12/glyco/conjugate#>
    prefix faldo: <http://www.biohackathon.org/resource/faldo/>
    prefix dcterms: <http://purl.org/dc/terms/>
    
    SELECT distinct ?ProteinAcc ?Position ?Saccharide ?Id
    WHERE {
      ?Glyco a gco:Referenced_glycoconjugate ; gco:has_protein_part ?Protein .
      ?Protein gco:has_protein ?ProteinAcc .
      ?Protein gco:has_saccharide_set ?Set .
      ?Set sio:is-component-part-of ?SetItem .
      ?SetItem owl:sameAs ?Saccharide .
      ?Saccharide dcterms:identifier ?Id .
      ?Protein gco:glycosylated_at ?Region .
      ?Region faldo:ExactPosition ?Faldo .
      ?Faldo faldo:position ?Position .
    }
    """


    def __init__(self, usecache=False):
        super(UniCarbKB, self).__init__(usecache=usecache)

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
        glycanType, key = url.split("/")[-2:]
        if glycanType not in ('structure', 'composition'):
            raise ValueError("Glycan URI is neither structure nor composition")
        return glycanType, key

    def getProtein(self, filepath):
        # The columns will be: UniProt ID, Position, PubMed ID, UniCarb ID, GlyTouCan ID, AminoAcid, TypeAminoAcid, Notes
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
        return self.getProtein(
            "./unicarb_temp/data_files_unicarbkb_QC_build_triplestore_2.0.4_human_unicarbkb_2_0_4_human.csv")

    def protein_mouse(self):
        return self.getProtein(
            "./unicarb_temp/data_files_unicarbkb_QC_build_triplestore_2.0.4_mouse_unicarbkb_2_0_4_mouse.csv")

    def getAllProtein(self):
        return self.protein_human() + self.protein_mouse()

    def unicarb2glytoucan(self):
        querykey = "unicarb2glytoucan"
        if self.usecache and querykey in self.cachedata:
            return self.cachedata[querykey]

        results = self.g.query(self.unicarb2glytoucan_query)
        res = defaultdict(set)
        for r in results.bindings:
            row = map(str, map(r.get, results.vars))
            gtype, uid = self.unicarbGlycanURL2id(row[0])
            if gtype != 'structure' or not uid:
                continue
            gid = row[1]
            res[uid].add(gid)
        if self.usecache:
            self.cachedata[querykey] = res
            self.cacheupdated = True
        return res

    def taxonomy(self):
        querykey = "taxonomy"
        if self.usecache and querykey in self.cachedata:
            return self.cachedata[querykey]
        results = self.g.query(self.taxonomy_query)
        res = defaultdict(set)
        for r in results.bindings:
            row = map(str, map(r.get, results.vars))
            gtype, uid = self.unicarbGlycanURL2id(row[0])
            if gtype != 'structure' or not uid:
                continue
            taxonomy = self.stripURL(row[1])
            res[uid].add(taxonomy)
        if self.usecache:
            self.cachedata[querykey] = res
            self.cacheupdated = True
        return res

    def protein_query_search(self):
        querykey = "fullpriteindetail"
        if self.usecache and querykey in self.cachedata:
            return self.cachedata[querykey]
        results = self.g.query(self.protein_query)
        res = defaultdict(set)
        for r in results.bindings:
            row = map(str, map(r.get, results.vars))
            print row
        if self.usecache:
            self.cachedata[querykey] = res
            self.cacheupdated = True
        return res

    idmapping = False

    def IDmapping(self, UniCarbID):
        if not self.idmapping:
            self.idmapping = self.unicarb2glytoucan()
        UniCarbID = str(UniCarbID)
        return self.idmapping[UniCarbID]


if __name__ == "__main__":

    import sys

    cmd = sys.argv.pop(1)

    if cmd.lower() == "taxa":

        uc = UniCarbKB()
        taxa = uc.taxonomy()
        for k, v in taxa.items():
            for vi in v:
                print k, vi

    elif cmd.lower() == "glytoucan":

        uc = UniCarbKB()
        gtc = uc.unicarb2glytoucan()
        for k, v in gtc.items():
            for vi in v:
                print k, vi

    else:
        print >> sys.stderr, "Bad command: %s" % (cmd,)
        sys.exit(1)
