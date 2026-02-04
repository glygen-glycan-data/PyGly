#!/bin/env python3.12

import sys

from urllib.parse import quote
from getwiki import GlycanData
w = GlycanData()

if len(sys.argv) >= 2:
    database = sys.argv[1]
else:
    database = "glycandatadev"

head = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE rdf:RDF[
        <!ENTITY rdf 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'>
        <!ENTITY rdfs 'http://www.w3.org/2000/01/rdf-schema#'>
        <!ENTITY swivt 'http://semantic-mediawiki.org/swivt/1.0#'>
        <!ENTITY glycandata 'http://glyomics.org/glycandata#'>
]>

<rdf:RDF
        xmlns:rdf="&rdf;"
        xmlns:rdfs="&rdfs;"
        xmlns:swivt="&swivt;"
        xmlns:glycandata="&glycandata;"
        xmlns:skos="http://www.w3.org/2004/02/skos/core#">
"""
tail = """
</rdf:RDF>
"""
glycantmpl = """
        <swivt:Subject rdf:about="http://glyomics.org/%(database)s/Special:URIResolver/%(accession)s">
                <rdf:type rdf:resource="http://glyomics.org/glycandata#Glycan" />
                <rdfs:label>%(accession)s</rdfs:label>
                <glycandata:accession rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(accession)s</glycandata:accession>
                </swivt:Subject>"""
annotationtmpl = """
        <swivt:Subject rdf:about="http://glyomics.org/%(database)s/Special:URIResolver/%(accession)s-23ANN_%(type)s.%(property_safe)s.%(source)s.%(sourceid)s">
                <rdf:type rdf:resource="http://glyomics.org/glycandata#Annotation" />
                <glycandata:hasglycan rdf:resource="http://glyomics.org/%(database)s/Special:URIResolver/%(accession)s" />
                <glycandata:property rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(property)s</glycandata:property>
                <glycandata:source rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(source)s</glycandata:source>
%(sourceidline)s                <glycandata:type rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(type)s</glycandata:type>
%(values)s
                </swivt:Subject>"""
sourceidlinetmpl = """                <glycandata:sourceid rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(sourceid)s</glycandata:sourceid>
"""

valuetmpl = """                <glycandata:value rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(value)s</glycandata:value>"""

def escape(s):
    return s.replace('&','&amp;').replace('<','&lt;').replace('>','&gt;')

sys.stdout.write(head)
for m in w.iterglycan():
    acc = m.get('accession')
    print(acc,file=sys.stderr)
    sys.stdout.write(glycantmpl%dict(database=database,accession=acc))
    for an in m.annotations():
        anval = an.get('value')
        values = []
        if isinstance(anval,list):
            for v in anval:
                values.append(valuetmpl%dict(value=escape(v)))
        else:
            values.append(valuetmpl%dict(value=escape(anval)))
        sourceid = an.get('sourceid',"")
        if sourceid:
            sourceidline = sourceidlinetmpl%dict(sourceid=escape(sourceid))
        else:
            sourceidline = ""
        sys.stdout.write(annotationtmpl%dict(database=database,accession=acc,
                                  type=an.get('type'), property=an.get('property'),
                                  property_safe=an.get('property').replace(' ','_'),
                                  source=an.get('source'), sourceid=escape(quote(sourceid,safe="")),
                                  sourceidline=sourceidline, values = '\n'.join(values)))
sys.stdout.write(tail)
