#!/bin/env python27

import os
import json
import copy
import findpygly
from pygly.GlycanFormatter import GlycoCTFormat, WURCS20Format



wurcs_parser = WURCS20Format()
glycoct_parser = GlycoCTFormat()

wpath = "dumps/wurcs"
gpath = "dumps/glycoct"

wlist = os.listdir(wpath)
glist = os.listdir(gpath)
alllist = list(set(wlist+glist))
print "Total glycan number %s" % len(alllist)

glycanobj = {}
for filename in alllist:
    acc = filename.rstrip(".txt")

    try:
        gseq = open(os.path.join(gpath, filename)).read().strip()
        obj = glycoct_parser.toGlycan(gseq)

    except:
        try:
            wseq = open(os.path.join(wpath, filename)).read().strip()
            obj = wurcs_parser.toGlycan(wseq)
        except:
            continue

    glycanobj[acc] = obj

print "Total parsable glycan number %s" % len(glycanobj)

from glycancomparison import GlycanCompatibleOneway, GlycanTopologySameAs
gco = GlycanCompatibleOneway()
gtsa = GlycanTopologySameAs()

# First compute the relationships
topomatchpool = []
compatibleres = {}

for acc1, g1 in glycanobj.items():

    ss = []
    comp = []

    for acc2, g2 in glycanobj.items():
        if acc1 == acc2:
            continue

        compflag = gco.get(g1, g2)
        topoflag = gtsa.get(g1, g2)

        if compflag:
            comp.append(acc2)

        if topoflag:
            topomatchpool.append(set([acc1, acc2]))

    if comp:
        compatibleres[acc1] = comp
"""
print len(topomatchpool)
print len(substructure)
print len(compatibleres)
"""

# Cluster the comparison results
anotherround = True
while (anotherround):
    anotherround = False
    for i, s1 in enumerate(topomatchpool):
        for j, s2 in enumerate(topomatchpool):

            if j <= i:
                continue

            if s1.intersection(s2):
                new = s1.union(s2)
                topomatchpool.append(new)
                anotherround = True
                topomatchpool.pop(j)
                topomatchpool.pop(i)
                break

        if anotherround:
            break

# print len(topomatchpool)


def viewerdatagenerater(component):
    data = {}
    edges = component["edges"]
    nodes = component["nodes"]
    root = component["root"]

    nodesd = {}
    edgesd = {}

    data["root"] = root
    data["name"] = root
    data["mw"] = "Unknown"

    for node in nodes:
        d = {}
        d["type"] = "Whatever"
        d["level"] = 0
        d["name"] = node
        if node == root:
            for node2 in nodes:
                if node2 != node:
                    break
            d["alternativeImageURL"] = "https://glytoucan.org/glycans/%s/image?style=compact&format=png&notation=cfg" %node2
            # print "https://glytoucan.org/glycans/%s/image?style=compact&format=png&notation=cfg"%node2
        nodesd[node] = d

    for node in edges.keys():
        l = []
        for edge in edges[node]:
            d = {}
            d["to"] = edge
            d["from"] = node
            d["type"] = "contains"
            l.append(d)

        edgesd[node] = l

    data["nodes"] = nodesd
    data["edges"] = edgesd

    json_data = json.dumps(data)
    return data

def removeShortcut(relationship):
    shortcut = {}
    res = copy.deepcopy(relationship)

    parents = relationship.keys()

    for p in parents:
        children = relationship[p]

        for c in children:
            try:
                grandchildren = relationship[c]
            except KeyError:
                continue

            for cc in relationship[c]:
                if cc in children:
                    if p in shortcut:
                        shortcut[p].append(cc)
                    else:
                        shortcut[p] = [cc]

    for k in shortcut.keys():
        shortcut[k] = sorted(list(set(shortcut[k])))

    for k, vchild in shortcut.items():
        for v in vchild:
            res[k].remove(v)

    return res



# For each cluster, generate the json
jsondata = {}
fakeroot = "Topology"
for i in topomatchpool:

    cluster = list(i)

    # Get the edge data
    link = {}
    for node in cluster:
        try:
            l = compatibleres[node]
            link[node] = l
        except KeyError:
            pass

    # Find all node in level 1 which is one level down from the fake composition
    allchildren = []
    for node in cluster:
        children = compatibleres.get(node, [])
        allchildren += children
    allchildren = list(set(allchildren))

    alllevel1node = []
    for node in cluster:
        if node not in allchildren:
            alllevel1node.append(node)

    alllevel1node.sort()
    representative = alllevel1node[0]

    # fakeroot = g.gettopo(cluster[0])
    link[fakeroot] = alllevel1node

    nodes = cluster + [fakeroot]

    link = removeShortcut(link)

    res = {
        "nodes": nodes,
        "edges": link,
        "root": fakeroot
    }

    data = viewerdatagenerater(res)
    jsondata[representative] = data

thingsinpool = set()
for s in topomatchpool:
    thingsinpool = thingsinpool.union(s)

for acc in glycanobj.keys():
    if acc in thingsinpool:
        continue

    nodes = [fakeroot, acc]
    link = {fakeroot: [acc]}

    res = {
        "nodes": nodes,
        "edges": link,
        "root": fakeroot
    }

    data = viewerdatagenerater(res)
    # print data
    jsondata[acc] = data


f1 = open("topology.json", "w")
# Make the JSON dump more predictable...
jsonstr = json.dumps(jsondata, sort_keys=True, indent=2, separators=(',', ': ')) + '\n'
f1.write(jsonstr)
f1.close()

print len(jsondata)
