#!/bin/env python27

import os, sys
import json
from GlycanFormatter import GlycoCTFormat, WURCS20Format
from GlyTouCan import GlyTouCan

g = GlyTouCan()

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

rpath = "dumps/redend"
reducing_end = {}
for filename in alllist:
    acc = filename.rstrip(".txt")
    r = open(os.path.join(rpath, filename)).read().strip()
    red = eval(r)
    if red is None:
        red = [True, False]

    reducing_end[acc] = red

from glycancomparison import GlycanCompatibleOneway, GlycanTopologySameAs, MotifSearchTopologicalSameAs
gco = GlycanCompatibleOneway()
gtsa = GlycanTopologySameAs()
mstsa = MotifSearchTopologicalSameAs()

topomatchpool = []

for acc1, g1 in glycanobj.items():

    for acc2, g2 in glycanobj.items():
        if acc1 == acc2:
            continue

        topoflag = gtsa.get(g1, g2)

        if topoflag:
            topomatchpool.append(set([acc1, acc2]))

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

# Add singleton to the pool
intopo = set()
for s in topomatchpool:
    intopo = intopo.union(s)

for glytoucanacc in glycanobj.keys():
    if glytoucanacc not in intopo:
        topomatchpool.append(set([glytoucanacc]))

# 3 informative dictionary
topologies = {}
cluster_re = {}
mono_num_dict = {}

acc = 1
for s in topomatchpool:
    topologies[acc] = list(s)
    acc += 1

for acc, l in topologies.items():
    r = []
    for node in l:
        r+=reducing_end[node]
    cluster_re[acc] = list(set(r))

for acc, l in topologies.items():
    glytoucan = l[0]
    num = len(list(glycanobj[glytoucan].all_nodes()))

    if num in mono_num_dict:
        mono_num_dict[num].append(acc)
    else:
        mono_num_dict[num] = [acc]

def emptyGraph(acc):
    data = {}
    root = "fakeroot"
    data["root"] = root
    data["name"] = root
    data["mw"] = "Unknown"

    nodes = {
        root: {"type": "whatever", "level": 0, "name": root, "hidden": True},
        acc: {"type": "whatever", "level": 0, "name": acc}
    }
    edges = {
        root: [{"to": acc, "from": root, "type": "contains", "hidden": True}]
    }

    data["nodes"] = nodes
    data["edges"] = edges
    return data


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
            d["hidden"] = True
        if node != gacc:
            a = [node]
            for l in topologies.values():
                if node in l:
                    a = l
                    break
            if len(a) < 2:
                d["label"] = "%s" % node
            else:
                d["label"] = "%s (Cluster)" % node
        nodesd[node] = d

    for node in edges.keys():
        l = []
        for edge in edges[node]:
            d = {}
            d["to"] = edge
            d["from"] = node
            d["type"] = "contains"
            if node == root:
                d["hidden"] = True
            l.append(d)

        edgesd[node] = l

    data["nodes"] = nodesd
    data["edges"] = edgesd

    #json_data = json.dumps(data)
    return data

# Compute the graph when extension and removal happens on the reducing end
json_data = {}
for gacc, gobj in glycanobj.items():

    num_mono = len(list(gobj.all_nodes()))

    parent_nodes = []
    if num_mono-1 in mono_num_dict:
        cluster_accs = mono_num_dict[num_mono-1]
        for cluster_acc in cluster_accs:
            cluster_representative = topologies[cluster_acc][0]
            motifobj = glycanobj[cluster_representative]
            #reducing_end_res = mstsa.reducingEndOnly(motifobj, gobj)
            non_reducing_end_res = mstsa.nonReducingEndOnly(motifobj, gobj)
            if non_reducing_end_res:
                parent_nodes.append(cluster_representative)

    child_nodes = []
    this_glycan_redend = reducing_end[gacc]
    if num_mono+1 in mono_num_dict:
        cluster_accs = mono_num_dict[num_mono + 1]

        for cluster_acc in cluster_accs:
            cluster_representative = topologies[cluster_acc][0]
            searchedobj = glycanobj[cluster_representative]

            #reducing_end_res = mstsa.reducingEndOnly(gobj, searchedobj)
            non_reducing_end_res = mstsa.nonReducingEndOnly(gobj, searchedobj)
            if non_reducing_end_res:
                child_nodes.append(cluster_representative)

    fakeroot = "fakeroot"
    data = {"root": fakeroot}
    data["nodes"] = parent_nodes + child_nodes + [gacc, fakeroot]

    data_edges = {}
    # print parent_nodes, gacc, child_nodes

    if len(parent_nodes) + len(child_nodes) == 0:
        json_data[gacc] = emptyGraph(gacc)
        continue
    else:
        if parent_nodes:
            data_edges[fakeroot] = parent_nodes
            for node in parent_nodes:
                data_edges[node] = [gacc]
        else:
            data_edges[fakeroot] = [gacc]
        if child_nodes:
            data_edges[gacc] = child_nodes
    data["edges"] = data_edges


    json_data[gacc] = viewerdatagenerater(data)


f1 = open("redonly.json", "w")
jsonstr = json.dumps(json_data) + '\n\n\n\n\n\n'
f1.write(jsonstr)
f1.close()

print len(json_data)


# Compute the graph when extension and removal happens on the non-reducing end
json_data = {}
for gacc, gobj in glycanobj.items():

    num_mono = len(list(gobj.all_nodes()))

    parent_nodes = []
    if num_mono-1 in mono_num_dict:
        cluster_accs = mono_num_dict[num_mono-1]
        for cluster_acc in cluster_accs:
            cluster_representative = topologies[cluster_acc][0]
            motifobj = glycanobj[cluster_representative]
            reducing_end_res = mstsa.reducingEndOnly(motifobj, gobj)
            #non_reducing_end_res = mstsa.nonReducingEndOnly(motifobj, gobj)
            if reducing_end_res:
                parent_nodes.append(cluster_representative)

    child_nodes = []
    this_glycan_redend = reducing_end[gacc]
    if num_mono+1 in mono_num_dict:
        cluster_accs = mono_num_dict[num_mono + 1]

        for cluster_acc in cluster_accs:
            cluster_representative = topologies[cluster_acc][0]
            searchedobj = glycanobj[cluster_representative]

            reducing_end_res = mstsa.reducingEndOnly(gobj, searchedobj)
            #non_reducing_end_res = mstsa.nonReducingEndOnly(gobj, searchedobj)
            if reducing_end_res:
                child_nodes.append(cluster_representative)

    fakeroot = "fakeroot"
    data = {"root": fakeroot}
    data["nodes"] = parent_nodes + child_nodes + [gacc, fakeroot]

    data_edges = {}
    # print parent_nodes, gacc, child_nodes

    if len(parent_nodes) + len(child_nodes) == 0:
        json_data[gacc] = emptyGraph(gacc)
        continue
    else:
        if parent_nodes:
            data_edges[fakeroot] = parent_nodes
            for node in parent_nodes:
                data_edges[node] = [gacc]
        else:
            data_edges[fakeroot] = [gacc]
        if child_nodes:
            data_edges[gacc] = child_nodes
    data["edges"] = data_edges


    json_data[gacc] = viewerdatagenerater(data)

f2 = open("nonredonly.json", "w")
jsonstr = json.dumps(json_data) + '\n\n\n\n\n\n'
f2.write(jsonstr)
f2.close()

print len(json_data)

