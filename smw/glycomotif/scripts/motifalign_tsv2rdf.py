#!/bin/env python27

import sys
import csv
import copy

motif_alignment_file_path_tsv = sys.argv[1]#"../data/motif_alignment.tsv"
motif_alignment_file_path_rdf = sys.argv[2]#"../data/motif_alignment.rdf"

motif_alignment = {}
motif_alignment_single = {
    "red_only": [],
    "other": []
}

with open(motif_alignment_file_path_tsv) as f:
    for r in csv.reader(f, delimiter="\t"):
        motif_acc, match, reducing_end_match = r

        if motif_acc not in motif_alignment:
            motif_alignment[motif_acc] = copy.deepcopy(motif_alignment_single)

        if reducing_end_match == "True":
            motif_alignment[motif_acc]["red_only"].append(match)
        elif reducing_end_match == "False":
            motif_alignment[motif_acc]["other"].append(match)


rdf_template = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE rdf:RDF[
    <!ENTITY rdf 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'>
    <!ENTITY rdfs 'http://www.w3.org/2000/01/rdf-schema#'>
    <!ENTITY swivt 'http://semantic-mediawiki.org/swivt/1.0#'>
    <!ENTITY glycomotif 'http://glycandata.glygen.org/glycomotif#'>
]>

<rdf:RDF
    xmlns:rdf="&rdf;"
    xmlns:rdfs="&rdfs;"
    xmlns:swivt="&swivt;"
    xmlns:glycomotif="&glycomotif;"
    xmlns:skos="http://www.w3.org/2004/02/skos/core#">
%s
</rdf:RDF>"""
rdf_s = ""

for acc in sorted(motif_alignment.keys()):
    print acc
    red = sorted(motif_alignment[acc]["other"])
    any = sorted(motif_alignment[acc]["red_only"] + motif_alignment[acc]["other"])

    rdf_per_motif = '\t<swivt:Subject rdf:about="http://glycandata.glygen.org/glycomotif/Special:URIResolver/GM.' + acc + '">\n%s\n\t</swivt:Subject>\n'
    rdf_normal_align = '\t\t<glycomotif:has_alignment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</glycomotif:has_alignment>'
    rdf_red_align = '\t\t<glycomotif:has_redend_alignment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</glycomotif:has_redend_alignment>'

    normal_align = []
    red_align = []

    for i in any:
        normal_align.append(rdf_normal_align%i)
        if i in red:
            red_align.append(rdf_red_align%i)

    rdf_s += rdf_per_motif % ("\n".join(normal_align) + "\n" + "\n".join(red_align))

rdf_f = open(motif_alignment_file_path_rdf, "w")
rdf_f.write(rdf_template%rdf_s)
