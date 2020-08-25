import os
import sys
import csv
import copy
import hashlib
import rdflib
import rdflib.plugins.serializers.rdfxml


motif_tsv = sys.argv[1] # "../data/motif_alignment.tsv"
prefix = sys.argv[2]
rdf_file_path = sys.argv[3] # "../data/alignment.rdf"

alignments = {}
alignments_template = {
    "Core": [],
    "Substructure": [],
    "Whole": [],
    "Nonreducing-End": []
}
for i, line in enumerate(csv.reader(open(motif_tsv), dialect="excel-tab")):
    if i == 0:
        continue

    # Motif Structure Core Substructure Whole Nonreducing_end
    macc = line[0]
    gacc = line[1]
    flags = map(lambda x:True if x=="Y" else False, line[2:])

    if macc not in alignments:
        alignments[macc] = copy.deepcopy(alignments_template)

    if flags[0]:
        alignments[macc]["Core"].append(gacc)

    if flags[1]:
        alignments[macc]["Substructure"].append(gacc)

    if flags[2]:
        alignments[macc]["Whole"].append(gacc)

    if flags[3]:
        alignments[macc]["Nonreducing-End"].append(gacc)

"""
for macc, res0 in alignments.items():
    print macc
    for alignments_type, res00 in res0.items():
        print alignments_type[:4], res00[:10]
"""

rdfgraph = rdflib.Graph()

rdfns = rdflib.RDF
rdfsns = rdflib.RDFS
skosns = rdflib.Namespace('http://www.w3.org/2004/02/skos/core#')
swivtns = rdflib.Namespace('http://semantic-mediawiki.org/swivt/1.0#')
glycomotifns = rdflib.Namespace('http://glycandata.glygen.org/glycomotif#')


rdfgraph.bind("rdf", rdfns)
rdfgraph.bind("rdfs", rdfsns)
rdfgraph.bind("glycomotif", glycomotifns)
rdfgraph.bind("swivt", swivtns)
rdfgraph.bind("skos", skosns)

abbr_extension = {
    "Core": "Core Alignment",
    "Substructure": "Substructure Alignment",
    "Whole": "Whole-Glycan Alignment",
    "Nonreducing-End": "Nonreducing-End Alignment"
}

from getwiki import GlycoMotifWiki
w = GlycoMotifWiki()

for m in w.itermotif():
    pageid = m.get("id")
    motifacc = m.get("glytoucan")
    alignment_allowed_type = m.get("alignment")

    print >> sys.stderr, pageid

    if motifacc not in alignments:
        print >> sys.stderr, "No alignments found for %s(%s)" % (pageid, motifacc)
        continue

    alignments_per_motif = alignments[motifacc]

    motif_rdf_node = rdflib.URIRef("http://glycandata.glygen.org/%s/Special:URIResolver/%s" % (prefix, pageid))
    rdfgraph.add((motif_rdf_node, rdfns.type, swivtns["Subject"]))

    for alignment_type_abbr, structure_accs in alignments_per_motif.items():
        alignment_type = abbr_extension[alignment_type_abbr]

        if alignment_type_abbr not in alignment_allowed_type:
            print >> sys.stderr, "Skipping alignment for %s %s" % (pageid, alignment_type)
            continue
        print >> sys.stderr, "Loading alignment for %s %s" % (pageid, alignment_type)

        for acc in structure_accs:
            tmp = "_".join((alignment_type_abbr, motifacc, acc))
            # alignment_id = hashlib.sha256(tmp).hexdigest()
            alignment_id = "_".join((alignment_type_abbr, pageid, acc))

            matched_rdf_node = glycomotifns[alignment_id]
            rdfgraph.add((matched_rdf_node, glycomotifns["structure_accession"], rdflib.Literal(acc)))
            rdfgraph.add((matched_rdf_node, glycomotifns["alignment_type"], rdflib.Literal(alignment_type)))
            rdfgraph.add((motif_rdf_node, glycomotifns["has_alignment"], matched_rdf_node))


writer = rdflib.plugins.serializers.rdfxml.PrettyXMLSerializer(rdfgraph, max_depth=2)
writer.serialize(open(rdf_file_path, "w"))








