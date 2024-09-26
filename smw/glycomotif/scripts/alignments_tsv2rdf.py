#!/bin/env python2
import os
import sys
import csv
import copy
import gzip
import hashlib
import rdflib
import rdflib.plugins.serializers.rdfxml


def get_indices(col):
    indices = False
    if col != 'N':
        indices = set(map(int,col[2:].split(',')))
    return(indices) 
        


motif_tsv = sys.argv[1]

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

ColumnMap = {
    "Core": ('Core_Inclusive','Core_Strict'),
    "Substructure": ('Substructure_Inclusive','Substructure_Strict'),
    "Whole-Glycan": ('Whole_Inclusive','Whole_Strict'),
    "Nonreducing-End": ('Non_Red_Inclusive','Non_Red_Strict')
}

for line in csv.DictReader(gzip.open(motif_tsv,'rt'), dialect="excel-tab"):
    motifacc = line['Motif']
    acc = line['Structure']
    for alignment_type in ColumnMap:
        loosevals = line[ColumnMap[alignment_type][0]].split(":")
        strictvals = line[ColumnMap[alignment_type][1]].split(":")
        if loosevals[0] == "Y":
            if strictvals[0] == "Y":
                strict = True
                indices = strictvals[1]
                linkindices = strictvals[2]
            else:
                strict = False
                indices = loosevals[1]
                linkindices = loosevals[2]

            alignment_id = "-".join((motifacc, alignment_type, acc))
            matched_rdf_node = glycomotifns[alignment_id]

            rdfgraph.add((matched_rdf_node, glycomotifns["motif_accession"], rdflib.Literal(motifacc)))
            rdfgraph.add((matched_rdf_node, glycomotifns["alignment_type"], rdflib.Literal(alignment_type)))
            rdfgraph.add((matched_rdf_node, glycomotifns["structure_accession"], rdflib.Literal(acc)))
            rdfgraph.add((matched_rdf_node, glycomotifns["structure_residue_ids"], rdflib.Literal(indices)))
            rdfgraph.add((matched_rdf_node, glycomotifns["structure_link_ids"], rdflib.Literal(linkindices)))
            rdfgraph.add((matched_rdf_node, glycomotifns["strict"], rdflib.Literal(str(strict).lower(), datatype=rdflib.XSD.boolean)))

writer = rdflib.plugins.serializers.rdfxml.PrettyXMLSerializer(rdfgraph, max_depth=2)
writer.serialize(sys.stdout.buffer)
