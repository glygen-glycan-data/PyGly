#!/bin/env python27
import sys,os,gzip,string,csv
from xml.etree import ElementTree as ET
from optparse import OptionParser
from zipfile import ZipFile, ZIP_DEFLATED
from collections import defaultdict
from StringIO import StringIO

parser = OptionParser()
parser.add_option('--glycomedb',type='string',default=None,dest='glydb',
                  help="GlycomeDB XML. Required.")
parser.add_option('--resource',type='string',default="",dest='resource',
                  help="Comma-separated list of resource identifiers to extract. Default: Extract all glycans.")
parser.add_option('--taxid',type='string',default="",dest='taxid',
                  help="Comma-separated list of taxonomy ids (and descendants) to extract. Default: Extract all glycans.")
parser.add_option('--out',type='string',default=None,dest='out',
                  help="Name of zip-file for GlycoCT-condenced format glycan structures. Required.")
opt,args = parser.parse_args()

if opt.glydb == None:
    parser.error("GlycomeDB XML file required")
if opt.out == None:
    parser.error("No output zip-file specified")

taxids = set(map(int,filter(None,map(string.strip,opt.taxid.split(',')))))
tt = None
if len(taxids) > 0:
    from taxonomy import NCBITaxonomyTree
    tt = NCBITaxonomyTree()

md = dict()

def process_structure(zf,ele):
    assert ele.attrib['database'] == 'GlycomeDB'
    id = ele.attrib['id']
    tids = set(map(int,map(lambda t: t.attrib.get('ncbi'),ele.findall('taxon'))))
    res = map(lambda t: (t.attrib.get('db'),t.attrib.get('id')),ele.findall('resource'))
    if opt.resource and len(filter(lambda t: t[0] == opt.resource, res)) == 0:
	return
    if len(taxids) > 0:
        keep = False
        for tid in tids:
            if tid in taxids:
                keep = True
                break
            anc,dummy = tt.get_ancestor_taxid(tid)
            # print set(anc)
            # print taxids
            if set(anc) & taxids:
                keep = True
                break
        if not keep:
            return
    zf.writestr("%s.txt"%id,ele.findtext('sequence'))
    if len(tids) > 0 or len(res) > 0:
        s = StringIO()
        w = csv.writer(s)
        # w.writerow(['id',id])
        if len(tids) > 0:
            w.writerow(['taxon']+["ncbi:%s"%tid for tid in tids])
        if len(res) > 0:
            w.writerow(['resource']+["%s:%s"%t for t in res])
        zf.writestr("%s.att"%id,s.getvalue())

zf = ZipFile(opt.out,'w',ZIP_DEFLATED)
h = gzip.open(opt.glydb)
for event, ele in ET.iterparse(h):
    if ele.tag == 'structure':
        process_structure(zf,ele)
        ele.clear()
h.close()
zf.close()
