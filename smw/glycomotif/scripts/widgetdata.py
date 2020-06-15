#!/bin/env python27

import os
import sys
import os
import json
import copy
import pygly.alignment as alignment
from collections import defaultdict
from pygly.GlycanFormatter import WURCS20Format, GlycoCTFormat



topology_file_path = sys.argv[1]
non_file_path = sys.argv[2]
red_file_path = sys.argv[3]



wp = WURCS20Format()
gp = GlycoCTFormat()



class RootMonosaccharideTopoLeq(alignment.MonosaccharideImageEqual):

    def leq(self, a, b):
        return self.eq(a, b)

class MonosaccharideTopoLeq(alignment.MonosaccharideTopoEqual):

    def leq(self, a, b):
        return self.eq(a, b)

class LinkageTopoLeq(alignment.LinkageTopoEqual):

    def leq(self, a, b):
        return self.eq(a, b)

class GlycanTopoEqualWithRootTweak(alignment.GlycanEquivalence):

    def __init__(self, **kw):
        kw['substcmp']=alignment.SubstituentEqual(**kw)
        kw['linkcmp']=alignment.LinkageTopoEqual(**kw)
        kw['sublinkcmp']=kw['linkcmp']

        # monotest needs a subst test and a sublink test
        kw['monocmp']=alignment.MonosaccharideTopoEqual(**kw)
        kw['rootmonocmp'] = alignment.MonosaccharideImageEqual(**kw)
        super(GlycanTopoEqualWithRootTweak,self).__init__(**kw)



class GlycanTopoSubstructureSearch(alignment.SubstructureSearch):

    def __init__(self, **kw):
        kw["rootmonocmp"] = RootMonosaccharideTopoLeq(
            substcmp=alignment.SubstituentEqual(),
            sublinkcmp=alignment.LinkageEqual()
        )
        kw["monocmp"] = MonosaccharideTopoLeq(
            substcmp=alignment.SubstituentEqual(),
            sublinkcmp=alignment.LinkageEqual()
        )
        kw["linkcmp"] = LinkageTopoLeq()
        super(GlycanTopoSubstructureSearch, self).__init__(**kw)


topology_comparison = GlycanTopoEqualWithRootTweak()
gtcm = alignment.GlyTouCanMotif()
gtsss = GlycanTopoSubstructureSearch()





from getwiki import GlycoMotifWiki, AllMotif


glycans = {}
supported_acc = []

glycan_length = {}

w = GlycoMotifWiki()
AllMotifpageid = AllMotif.id

for m in w.itermotif():

    acc = m.get("glytoucan")
    if acc in glycans:
        continue

    try:
        glycans[acc] = wp.toGlycan(str(m.get("wurcs")))
    except:

        try:
            glycans[acc] = gp.toGlycan(m.get("glycoct"))
        except:
            continue

    g = glycans[acc]
    l = len(list(g.all_nodes()))

    glycan_length[acc] = l

print "%s motifs are supported" % len(glycans)

supported_acc = list(glycans.keys())
topology_pool = []

for i in range(len(supported_acc)):
    for j in range(i, len(supported_acc)):
        acc1 = supported_acc[i]
        acc2 = supported_acc[j]

        g1 = glycans[acc1]
        g2 = glycans[acc2]

        if topology_comparison.eq(g1, g2):
            topology_pool.append({acc1, acc2})


run = True
while run:
    run = False

    bk = False
    for i in range(len(topology_pool)):
        for j in range(i+1, len(topology_pool)):

            a = topology_pool[i]
            b = topology_pool[j]
            intersection = a.intersection(b)

            if len(intersection) > 0:
                bk = True
                run = True

                topology_pool.pop(j)
                topology_pool.pop(i)
                topology_pool.append(a.union(b))
                break

        if bk:
            break



def viewerdatageneraterTopo(component):
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

    return data


relationship_within_topology_cluster = {}
for i, cluster in enumerate(topology_pool):

    relation = {}
    for acc in cluster:
        relation[acc] = []

    for p in cluster:
        for c in cluster:
            if p == c:
                continue

            tmp = gtcm.leq(glycans[p], glycans[c])
            if tmp:
                relation[p].append(c)

    relation = dict(relation)
    for p in relation.keys():
        children = relation[p]
        for c in relation[p]:
            grandchildren = relation[c]
            for rm in set(children).intersection(grandchildren):
                relation[p].remove(rm)

    #relationship_within_topology_cluster[i] = relation
    allChildren = set(reduce(lambda x,y:x+y, relation.values()))
    rootDirectChildren = list(set(relation.keys()) - allChildren)

    TopoRoot = "Topology"
    data = {"root": TopoRoot}
    data["nodes"] = relation.keys() + [TopoRoot]
    data["edges"] = {TopoRoot: rootDirectChildren}
    for p, c in relation.items():
        if len(c) > 0:
            data["edges"][p] = c

    comp = viewerdatageneraterTopo(data)

    relationship_within_topology_cluster[relation.keys()[0]] = comp



topology_cluster_by_length = defaultdict(list)
topology_cluster_length = {}
for i, cluster in enumerate(topology_pool):

    accx = list(cluster)[0]
    l = glycan_length[accx]

    topology_cluster_by_length[l].append(i)
    topology_cluster_length[i] = l


relationship_red_end_topology_navigator = {}
relationship_non_red_end_topology_navigator = {}

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

        if node == PseudoRoot:
            d["hidden"] = True
        elif node != acc:
            a = None
            for l in topology_pool:
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


for acc in supported_acc:
    # red and non here refer to how the matches are done

    red_parents = []
    red_children = []

    non_parents = []
    non_children = []

    l = glycan_length[acc]

    if l-1 in topology_cluster_by_length:
        for cluster_id in topology_cluster_by_length[l-1]:
            represent_acc = list(sorted(topology_pool[cluster_id]))[0]

            flagr = gtsss.leq(glycans[represent_acc], glycans[acc], rootOnly=True, anywhereExceptRoot=False)
            flagn = gtsss.leq(glycans[represent_acc], glycans[acc], rootOnly=False, anywhereExceptRoot=True)

            if flagr:
                red_parents.append(represent_acc)

            if flagn:
                non_parents.append(represent_acc)

    if l+1 in topology_cluster_by_length:
        for cluster_id in topology_cluster_by_length[l+1]:
            represent_acc = list(sorted(topology_pool[cluster_id]))[0]

            flagr = gtsss.leq(glycans[acc], glycans[represent_acc], rootOnly=True, anywhereExceptRoot=False)
            flagn = gtsss.leq(glycans[acc], glycans[represent_acc], rootOnly=False, anywhereExceptRoot=True)

            if flagr:
                red_children.append(represent_acc)
            if flagn:
                non_children.append(represent_acc)

    #print
    #print red_parents, acc, red_children
    #print non_parents, acc, non_children

    parents = red_parents
    children = red_children

    PseudoRoot = "PseudoRoot"
    data = {"root": PseudoRoot}
    data["nodes"] = parents + children + [acc, PseudoRoot]
    data["edges"] = {}
    comp = emptyGraph(acc)

    if len(parents) + len(children) != 0:

        if len(parents) > 0:
            data["edges"][PseudoRoot] = parents
            for node in parents:
                data["edges"][node] = [acc]
        else:
            data["edges"][PseudoRoot] = [acc]
        if children:
            data["edges"][acc] = children

        comp = viewerdatagenerater(data)

    relationship_non_red_end_topology_navigator[acc] = comp




    parents = non_parents
    children = non_children

    PseudoRoot = "PseudoRoot"
    data = {"root": PseudoRoot}
    data["nodes"] = parents + children + [acc, PseudoRoot]
    data["edges"] = {}
    comp = emptyGraph(acc)

    if len(parents) + len(children) != 0:

        if len(parents) > 0:
            data["edges"][PseudoRoot] = parents
            for node in parents:
                data["edges"][node] = [acc]
        else:
            data["edges"][PseudoRoot] = [acc]
        if children:
            data["edges"][acc] = children

        comp = viewerdatagenerater(data)

    relationship_red_end_topology_navigator[acc] = comp


json.dump(relationship_within_topology_cluster, open(topology_file_path, "w"), sort_keys=True, indent=2)
json.dump(relationship_non_red_end_topology_navigator, open(red_file_path, "w"), sort_keys=True, indent=2)
json.dump(relationship_red_end_topology_navigator, open(non_file_path, "w"), sort_keys=True, indent=2)







# Load topology cluster attribute
w = GlycoMotifWiki()
AllMotifpageid = AllMotif.id
for m in w.itermotif():

    if m.get("collection") != AllMotifpageid:
        continue

    gtcid = m.get("glytoucan")
    equivalents = []
    for topos in topology_pool:
        if gtcid in topos:
            equivalents = topos[:]
            break

    equivalents = map(lambda x: str(AllMotifpageid + "." + x), equivalents)

    print "Set topology cluster attribute for ", m.get("id")
    m.set("topology", equivalents)
    w.update(m)


