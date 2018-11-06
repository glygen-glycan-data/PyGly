#!/bin/env python27

import sys, traceback, os, csv

from getwiki import GlycoMotifWiki
from getwiki import GlydinMotif, GlydinCummingsMotif, GlydinHayesMotif, GlydinCermavMotif, GlydinSugarbindMotif, GlydinBioligoMotif

w = GlycoMotifWiki()

from dataset import XLSXFileTable
rows = XLSXFileTable(sys.argv[1])

from pygly.GlyTouCan import GlyTouCan
gtc = GlyTouCan()

from gtccache import GlyTouCanCache
gtccache = GlyTouCanCache()

"""
[x]collection structure
key: id (not row number)
value: {"aglycon":"String"/None, "name":""/None, "gtcid":""}
"""

ccollection = {}
hcollection = {}
vcollection = {}
sbcollection = {}
bcollection = {}
collections = {"c": ccollection, "h": hcollection, "v": vcollection, "sb": sbcollection, "b": bcollection}

def slash_handler(id, aglycon, name, gtcid):

    idlist = id.split("//")
    aglyconlist = aglycon.split("//")
    # Map the empty string to None
    aglyconlist = map(lambda x: None if not x else x, aglyconlist)

    if name is None:
        namelist = len(idlist) * [None]
    else:
        namelist = name.split("//")
        namelist = map(lambda x: None if not x else x, namelist)

    res = {}
    for i in range(len(idlist)):
        id0 = int(idlist[i])
        aglycon0 = aglyconlist[i]
        name0 = namelist[i]
        obj = {"aglycon": aglycon0, "name": name0, "gtcid": gtcid}
        res[id0] = obj

    return res

gtcid2glycoct = dict()
current = set()
for r in rows:

    index = int(r['ID'])  # which is also row number (it is row number in first column, not excel row number)
    glycoct = str(r['GlycoCT']).replace("\n\n", "\n").strip()

    cid = r['cummings ID']
    caglycon = r['Cummings Aglycon']
    cname = r['Cummings Name']

    hid = r['Hayes ID']
    haglycon = r['Hayes Aglycon']
    hname = r['Hayes Name']

    vid = r['cermav ID']
    vaglycon = r['cermav Aglycon']
    vname = r['Cermav Name']

    sbid = r['sugarbind ID']
    sbaglycon = r['sugarbind Aglycon']
    sbname = r['Sugarbind Name']

    bid = r['bioligo ID']
    baglycon = r['bioligo Aglycon']
    bname = None

    c = [cid, caglycon, cname]
    h = [hid, haglycon, hname]
    v = [vid, vaglycon, vname]
    sb = [sbid, sbaglycon, sbname]
    b = [bid, baglycon, bname]
    content = {"c": c, "h": h, "v": v, "sb": sb, "b": b}

    accession = "R%06d" % int(index)
    gtcid=gtccache.id2gtc(GlydinMotif.id + "." + accession)

    if gtcid == None:
        try:
            gtcid, isnew = gtc.register(glycoct)
        except:
            # traceback.print_exc()
            continue

    gtcid2glycoct[gtcid] = glycoct

    for key in content.keys():

        id0 = content[key][0]
        aglycon0 = content[key][1]
        name0 = content[key][2]
        collection = collections[key]

        if type(id0) == type(0):
            obj = {"aglycon": aglycon0, "name": name0, "gtcid": gtcid}
            collection[id0] = obj

        elif type(id0) == type(u"a"):
            objs = slash_handler(id0, aglycon0, name0, gtcid)
            for eachid0 in objs.keys():
                collection[eachid0] = objs[eachid0]

    # Load Glydin motif in this loop
    motif = GlydinMotif(accession=accession, name=None, glytoucan=gtcid, redend=None, aglycon=None, 
			glycoct=gtccache.gtc2glycoct(gtcid), wurcs=gtccache.gtc2wurcs(gtcid))
    # if GlyTouCan did not provide the GlycoCT, use ours...
    if not motif.get('glycoct'):
	motif.set('glycoct',gtcid2glycoct[gtcid])
    if w.put(motif):
        print >>sys.stderr, "Glydin:", accession
    current.add(accession)

for m in w.itermotif(collection=GlydinMotif):
    if m.get('accession') not in current:
        print >>sys.stderr, "Deleting:", m.get('pagename')
        w.delete(m.get('pagename'))


classnames2 = {"c": GlydinCummingsMotif, "h": GlydinHayesMotif, "v": GlydinCermavMotif, "sb": GlydinSugarbindMotif, "b": GlydinBioligoMotif}
aglycon2stdaglycon = {
    "Ser/Thr O-Mannose": "Ser/Thr",
    "Ceramide": "Cer",
    "Cer 9-O-Acetyl": "Cer",
}
stdaglycon = ["Ser/Thr", "Cer", "R", "Other"]
reaglycon = ["Ser/Thr", "Cer", "Other"]

for singleletter in classnames2.keys():
    motifClass = classnames2[singleletter]
    collection = collections[singleletter]
    current = set()
    print >>sys.stderr, "Loading %s Motif\n" % singleletter

    for acc in sorted(collection.keys()):
        data = collection[acc]

        accession = "%06d" % int(acc)
        glytoucan = data["gtcid"]
        aglycon = data["aglycon"]
        name = data["name"]
        redend = False

        if aglycon:
            aglycon = aglycon.strip()
            if aglycon in aglycon2stdaglycon:
                aglycon = aglycon2stdaglycon[aglycon]
            elif aglycon in stdaglycon:
                pass
            else:
                aglycon = "Other"

            if aglycon in reaglycon:
                redend = True
        else:
            aglycon = None
	    redend = None

        motif = motifClass(accession=accession, name=name, glytoucan=glytoucan, redend=redend, aglycon=aglycon, 
	                   glycoct=gtccache.gtc2glycoct(glytoucan), wurcs=gtccache.gtc2wurcs(glytoucan))
	# if GlyTouCan did not provide the GlycoCT, use ours...
        if not motif.get('glycoct'):
	    motif.set('glycoct',gtcid2glycoct[glytoucan])
        if w.put(motif):
            print >>sys.stderr, accession
        current.add(accession)

    for m in w.itermotif(collection=motifClass):
        if m.get('accession') not in current:
            print >>sys.stderr, "Deleting:", m.get('pagename')
            w.delete(m.get('pagename'))
