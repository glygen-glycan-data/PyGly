#!/bin/env python27

import sys
import csv
import copy


motif_alignment_file_path_tsv = sys.argv[1]#"../data/motif_alignment.tsv"
prefix = sys.argv[2]


motif_alignment = {}
alignment_keys = ["core", "substructure", "whole", "nred"]
motif_alignment_single = {}
for i in alignment_keys:
    motif_alignment_single[i] = []


with open(motif_alignment_file_path_tsv) as f:

    for i, r in enumerate(csv.reader(f, delimiter="\t")):
        if i == 0:
            continue

        motif_acc, structure = r[0:2]
        flags = map(lambda x: True if x == "Y" else False, r[2:])
        core, substructure, whole, nred = flags

        if motif_acc not in motif_alignment:
            motif_alignment[motif_acc] = copy.deepcopy(motif_alignment_single)

        if core:
            motif_alignment[motif_acc]["core"].append(structure)

        if substructure:
            motif_alignment[motif_acc]["substructure"].append(structure)

        if whole:
            motif_alignment[motif_acc]["whole"].append(structure)

        if nred:
            motif_alignment[motif_acc]["nred"].append(structure)




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
    print >>sys.stderr, acc


    rdf_per_motif = '\t<swivt:Subject rdf:about="http://glycandata.glygen.org/%s/Special:URIResolver/GM.'%(prefix,) + acc + '">\n%s\n\t</swivt:Subject>\n'

    rdf_alignment_template = {
        "core": '\t\t<glycomotif:has_core_alignment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</glycomotif:has_core_alignment>',
        "substructure": '\t\t<glycomotif:has_substructure_alignment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</glycomotif:has_substructure_alignment>',
        "whole": '\t\t<glycomotif:has_whole_glycan_alignment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</glycomotif:has_whole_glycan_alignment>',
        "nred": '\t\t<glycomotif:has_non-reducing_end_alignment rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</glycomotif:has_non-reducing_end_alignment>',
    }
    rdf_alignments = []

    for key in alignment_keys:
        matched_structures = sorted(motif_alignment[acc][key])
        for gacc in matched_structures:
            rdf_alignments.append(rdf_alignment_template[key] % gacc)

    rdf_s += rdf_per_motif % ("\n".join(sorted(rdf_alignments)))

res = rdf_template%rdf_s
sys.stdout.write(res)


