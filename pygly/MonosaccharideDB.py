import sys
import xml
import urllib
import urllib2
import copy
import xml.etree.ElementTree as ET
from GlyTouCan import *
from GlycanFormatter import GlycoCTFormat

GlycoctParser = GlycoCTFormat()
gtc = GlyTouCan()

class MonosaccharideDB:
    
    page_cache_by_id = {}
    page_cache_by_glycoct = {}

    def __init__(self):
        pass
    
    def encode(self, d):
        return urllib.urlencode(d)
    
    def get_document_by_glycoct(self, basetype, substlist):
        apiurl = "http://www.monosaccharidedb.org/display_monosaccharide.action?"
        
        param = {"scheme": "glycoct",
                 "name": basetype,
                 "output": "xml"}
        apiurl += urllib.urlencode(param)
        
        for subst in substlist:
            params = subst
            apiurl += "&%s" % urllib.urlencode(params)
        
        # print apiurl
        
        req = urllib2.Request(apiurl.replace("output=xml&", ""))
        httpdoc = urllib2.urlopen(req).read()
        req = urllib2.Request(apiurl)
        xmldoc = urllib2.urlopen(req).read()
        # Well, sometimes you have to request the URL twice to load IUPAC name
        # I don't know why what happens either.
    
        root = ET.fromstring(xmldoc)
        
        return root

    def get_document_by_mono_object(self, m):
        s = GlycoctParser.mtoStr(m)
        if s in self.page_cache_by_glycoct:
            return copy.deepcopy(self.page_cache_by_glycoct[s])
        else:
            basetype, substlist = GlycoctParser.mtodict(m)
            docroot = self.get_document_by_glycoct(basetype, substlist)
            self.page_cache_by_glycoct[s] = copy.deepcopy(docroot)
            
            return docroot
        
    def get_document_by_id(self, id):
        if id in self.page_cache_by_id:
            return copy.deepcopy(self.page_cache_by_id[id])
        else:
            apiurl = "http://www.monosaccharidedb.org/display_monosaccharide.action?id=%s&output=xml" % id
            req = urllib2.Request(apiurl.replace("&output=xml", ""))
            httpdoc = urllib2.urlopen(req).read()
            req = urllib2.Request(apiurl)
            xmldoc = urllib2.urlopen(req).read()
            root = ET.fromstring(xmldoc)
            self.page_cache_by_id[id] = copy.deepcopy(root)
            
            return root
    
    def parse_doc(self, root):
        res = {}
        for child in root:
            if child.tag in ["msdb_id", "msdb_name", "composition", "is_orientation_changed"]:
                res[child.tag] = child.text
            elif child.tag == "mass":
                res[child.tag] = child.attrib
            elif child.tag == "basetype":
                monobasetype = {}
                for c in child:
                    if c.tag in ["name", "size", "anomeric", "stereocode", "configuration", "composition"]:
                        monobasetype[c.tag] = c.text
                    elif c.tag in ["core_modifications", "mass"]:
                        monobasetype[c.tag] = c.attrib
                    elif c.tag == "ring":
                        ring = {}
                        for d in c:
                            ring[d.tag] = d.text
                        monobasetype["ring"] = ring
                    else:
                        raise IndexError 
                res[child.tag] = monobasetype
            elif child.tag == "synonyms":
                synonyms = []
                for c in child:
                    sym = c.attrib
                    sym["external_substituent"] = []
                    for gc in c:
                        sym["external_substituent"].append(gc.attrib)
                    synonyms.append(sym)
                res[child.tag] = synonyms
            elif child.tag == "substitutions":
                subst = []
                for sub in child:
                    d = {}
                    for sinfo in sub:
                        if sinfo.tag == "name":
                            d["name"] = sinfo.text
                        if sinfo.tag == "linkage":
                            for k, v in sinfo.attrib.items():
                                d[k] = v
                    subst.append(d)
                res[child.tag] = subst
            elif child.tag in ["atoms", "possible_linkage_positions"]:
                # Not really important at this moment
                pass
            else:
                print child.tag, child.attrib, child.text
                raise IndexError
        
        return res
    
    def get_IUPAC_synonyms(self, d):
        res = []
        try:
            for s in d['synonyms']:
                if s["scheme"] == "IUPAC":
                    res.append(s["name"])
            if res:
                return res
        except KeyError:
            return None

    def summary(self, acc):
        g = gtc.getGlycan(acc)
        for m in g.all_nodes():
            d = GlycoctParser.mtodict(m)
            root = self.get_document_by_glycoct(d[0], d[1])
            d = self.parse_doc(root)
            print "GlycoCT format: %s, MonosaccharideDB ID: %s, IUPAC synonyms: %s" % (GlycoctParser.mtoStr(m), d["msdb_id"], self.get_IUPAC_synonyms(d))


if __name__ == "__main__":
    m = MonosaccharideDB()
    gtcacc = sys.argv[1]
    # gtcacc = "G24536EZ"
    m.summary(gtcacc)


