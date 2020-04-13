#!/bin/env python27

from getwiki import GlycanDataDiskCache

import sys

def xmlescape(s):
    s = s.replace('&','&amp;')
    s = s.replace('<','&lt;')
    s = s.replace('>','&gt;')
    s = s.replace('\'','&apos;')
    s = s.replace('"','&quot;')
    return s

d = GlycanDataDiskCache(sys.argv[1])

prefix = sys.argv[2]

print """
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE rdf:RDF[
        <!ENTITY rdf 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'>
        <!ENTITY rdfs 'http://www.w3.org/2000/01/rdf-schema#'>
        <!ENTITY swivt 'http://semantic-mediawiki.org/swivt/1.0#'>
        <!ENTITY glycandata 'http://glycandata.glygen.org/glycandata#'>
]>

<rdf:RDF
        xmlns:rdf="&rdf;"
        xmlns:rdfs="&rdfs;"
        xmlns:swivt="&swivt;"
        xmlns:glycandata="&glycandata;"
        xmlns:skos="http://www.w3.org/2004/02/skos/core#">
""".strip('\n')

for g in d.iterglycan():
    glycan = g.get('id')
    print ("""
        <swivt:Subject rdf:about="http://glycandata.glygen.org/%(prefix)s/Special:URIResolver/%(glycan)s">
                <rdf:type rdf:resource="http://glycandata.glygen.org/glycandata#Glycan" />
                <rdfs:label>%(glycan)s</rdfs:label>
                <glycandata:accession rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(glycan)s</glycandata:accession>
                </swivt:Subject>
"""%dict(glycan=glycan,prefix=prefix)).strip('\n')
    for a in g.annotations():
	d = dict(glycan=glycan,prefix=prefix)
	d['property'] = a.get('property')
	d['type'] = a.get('type')
	d['source'] = a.get('source')
	if a.get('sourceid'):
	    d['sourceid'] = a.get('sourceid')
	    d['sourcewithsourceid'] = d['source']+"."+d['sourceid']
	else:
	    d['sourcewithsourceid'] = d['source']
	if isinstance(a.get('value'),basestring):
	    d['value'] = [ a.get('value') ]
	else:
	    d['value'] = a.get('value')
	for k in ('property','type','source','sourcewithsourceid'):
	    d[k+"_no_space"] = d[k].replace(' ','_')
	print ("""
        <swivt:Subject rdf:about="http://glycandata.glygen.org/%(prefix)s/Special:URIResolver/%(glycan)s.%(type)s.%(property_no_space)s.%(sourcewithsourceid)s">
                <rdf:type rdf:resource="http://glycandata.glygen.org/glycandata#Annotation" />
                <glycandata:hasglycan rdf:resource="http://glycandata.glygen.org/%(prefix)s/Special:URIResolver/%(glycan)s" />
                <glycandata:property rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(property)s</glycandata:property>
                <glycandata:type rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(type)s</glycandata:type>
                <glycandata:source rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(source)s</glycandata:source>
"""%d).strip('\n')
	if d.get('sourceid'):
	    print ("""
                <glycandata:sourceid rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%(sourceid)s</glycandata:sourceid>
"""%d).strip('\n')
	for val in d.get('value'):
	    print ("""
                <glycandata:value rdf:datatype="http://www.w3.org/2001/XMLSchema#string">%s</glycandata:value>
"""%(xmlescape(val),)).strip('\n')
	print ("""
                </swivt:Subject>
""").strip('\n')

print """
</rdf:RDF>
""".strip('\n')
