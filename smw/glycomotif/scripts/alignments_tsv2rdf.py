#!/bin/env python2
import os
import sys
import csv
import copy
import hashlib
import rdflib
import rdflib.plugins.serializers.rdfxml


motif_tsv = sys.argv[1]
# prefix = sys.argv[2]
# rdf_file_path = sys.argv[3]
"""
motif_tsv = "../data/motif_alignment.tsv"
prefix = "sth"
rdf_file_path = "../data/alignment.rdf"
"""

alignments = {}
alignments_template = {
    "Core": [],
    "Substructure": [],
    "Whole-Glycan": [],
    "Nonreducing-End": []
}
for i, line in enumerate(csv.reader(open(motif_tsv), dialect="excel-tab")):
    if i == 0:
        continue

    # Motif Structure Core Substructure Whole Nonreducing_end
    macc = line[0]
    gacc = line[1]
    flags = map(lambda x: True if x == "Y" else False, line[2:])

    if macc not in alignments:
        alignments[macc] = copy.deepcopy(alignments_template)

    if flags[0]:
        alignments[macc]["Core"].append([gacc, flags[4]])

    if flags[1]:
        alignments[macc]["Substructure"].append([gacc, flags[5]])

    if flags[2]:
        alignments[macc]["Whole-Glycan"].append([gacc, flags[6]])

    if flags[3]:
        alignments[macc]["Nonreducing-End"].append([gacc, flags[7]])


rdfgraph = rdflib.Graph()

rdfns = rdflib.RDF
rdfsns = rdflib.RDFS
skosns = rdflib.Namespace('http://www.w3.org/2004/02/skos/core#')
swivtns = rdflib.Namespace('http://semantic-mediawiki.org/swivt/1.0#')
glycomotifns = rdflib.Namespace('http://glyomics.org/glycomotif#')


rdfgraph.bind("rdf", rdfns)
rdfgraph.bind("rdfs", rdfsns)
rdfgraph.bind("glycomotif", glycomotifns)
rdfgraph.bind("swivt", swivtns)
rdfgraph.bind("skos", skosns)


for motifacc, alignments_per_motif in alignments.items():
    # print motifacc
    # motif_rdf_node = rdflib.URIRef("http://glyomics.org/%s/Special:URIResolver/GM.%s" % (prefix, motifacc))
    # rdfgraph.add((motif_rdf_node, rdfns.type, swivtns["Subject"]))

    for alignment_type, structure_accs in alignments_per_motif.items():

        for pair in structure_accs:

            acc, strict = pair
            alignment_id = "-".join((motifacc, alignment_type, acc))

            matched_rdf_node = glycomotifns[alignment_id]

            rdfgraph.add((matched_rdf_node, glycomotifns["motif_accession"], rdflib.Literal(motifacc)))
            rdfgraph.add((matched_rdf_node, glycomotifns["alignment_type"], rdflib.Literal(alignment_type)))
            rdfgraph.add((matched_rdf_node, glycomotifns["structure_accession"], rdflib.Literal(acc)))
            if strict:
                rdfgraph.add((matched_rdf_node, glycomotifns["strict"], rdflib.Literal("true" , datatype=rdflib.XSD.boolean)))
            else:
                rdfgraph.add((matched_rdf_node, glycomotifns["strict"], rdflib.Literal("false", datatype=rdflib.XSD.boolean)))

writer = rdflib.plugins.serializers.rdfxml.PrettyXMLSerializer(rdfgraph, max_depth=2)
# writer.serialize(open(rdf_file_path, "w"))
writer.serialize(sys.stdout)








