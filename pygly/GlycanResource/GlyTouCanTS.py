
from TripleStoreResource import TripleStoreResource
from GlycanResourceWrappers import partitioner, prefetcher

import hashlib
import base64
from collections import defaultdict

class GlyTouCanTS(TripleStoreResource):

    endpt = "http://ts.glytoucan.org/sparql"
    defns = "http://rdf.glycoinfo.org/glycan/"
    # verbose = True
    cachefile = ".gtc.cache"
    # usecache = True
    # prefetch = True
    # cachemode= 'c'
    
    sequence_formats = set(["wurcs", "glycoct", "iupac_extended", "iupac_condensed"])

    crossref_resources = set(['glycosciences_de', 'pubchem', 'kegg',
                              'unicarbkb', 'glyconnect', 'glycome-db',
                              'unicarb-db', 'carbbank', 'pdb', 'cfg',
                              'bcsdb','matrixdb','glycoepitope'])

    def __init__(self,**kw):
        if 'usecache' not in kw:
	    kw['usecache'] = True
	if 'prefetch' not in kw:
	    kw['prefetch'] = True
	if 'cachemode' not in kw:
	    kw['cachemode'] = 'c'
        super(GlyTouCanTS,self).__init__(**kw)
        for k in self.keys():
	    if k == "query_hashedseq":
                self.modify_method(k,partitioner(kwarg="hash",fmt="%%0%dx.*",values='hexidecimal'))
	    else:
                self.modify_method(k,partitioner())
	    if self._prefetch:
		if k == "query_hashedseq":
                    self.modify_method(k,prefetcher("hash",usecache=self._usecache))
		else:
                    self.modify_method(k,prefetcher(usecache=self._usecache))

    def getseq(self,accession,format='wurcs'):
        assert format in self.sequence_formats
        for row in self.query_sequence(accession=accession,format=format):
	    if row['format'] == 'wurcs' and not row['sequence'].startswith('WURCS'):
		continue
	    if row['format'] == 'glycoct' and not row['sequence'].startswith('RES'):
		continue
            return row['sequence']
        return None

    def allseq(self,format=None):
        assert format == None or format in self.sequence_formats
        for row in self.query_sequence(format=format):
	    if row['format'] == 'wurcs' and not row['sequence'].startswith('WURCS'):
                continue
	    if row['format'] == 'glycoct' and not row['sequence'].startswith('RES'):
                continue
            yield row['accession'], row['format'], row['sequence']

    def getmass(self,accession):
        for row in self.query_mass(accession=accession):
            try:
                return float(row['mass'])
            except ValueError:
                pass
        return None

    def allmass(self):
        for row in self.query_mass():
            try:
                yield row['accession'],float(row['mass'])
            except ValueError:
                pass

    def getmonocount(self,accession):
        for row in self.query_monocount(accession=accession):
            try:
                return int(row['count'])
            except ValueError:
                pass
        return None

    def allmonocount(self):
        for row in self.query_monocount():
            try:
                yield row['accession'],int(row['count'])
            except ValueError:
                pass

    def getrefs(self,accession):
        refs = set()
        for row in self.query_references(accession=accession):
            try:
                refs.add(int(row['ref']))
            except ValueError:
                continue
        return sorted(refs)

    def allrefs(self):
        for row in self.query_references():
            try:
                yield row['accession'],int(row['ref'])
            except ValueError:
                continue

    def getcrossrefs(self,accession,resource=None):
        assert resource == None or resource in self.crossref_resources
        for row in self.query_crossrefs(accession=accession,resource=resource):
            if row['resource'] in self.crossref_resources:
                yield row['resource'],row['entry']

    def allcrossrefs(self,resource=None):
        assert resource == None or resource in self.crossref_resources
        for row in self.query_crossrefs(resource=resource):
            if row['resource'] in self.crossref_resources:
                yield row['accession'],row['resource'],row['entry']

    # Named for consistency with GlyTouCan class...
    def getmotif(self,accession):
        return [ row['motif'] for row in self.query_motifs(accession=accession) ]

    def allmotifaligns(self):
        for row in self.query_motifs():
            yield row['accession'],row['motif']

    def allmotifs(self):
        for row in self.query_allmotif():
            yield row['accession'],row['label'],row['redend']

    def exists(self,accession):
        for row in self.query_exists(accession=accession):
            return True
        return False

    def allvalid(self):
	seen = set()
	for row in self.query_validacc1():
	    if row['validseqhash'] == row['seqhash']:
		if row['accession'] in seen:
		    continue
		seen.add(row['accession'])
		if self.invalid(row['accession']):
		    continue
	        yield row['accession']

    def allinvalid(self):
	seen = set()
	for row in self.query_validacc1():
	    if row['validseqhash'] != row['seqhash']:
		if row['accession'] in seen:
		    continue
		seen.add(row['accession'])
	        yield row['accession']

    def replace(self):
	vsh2acc = defaultdict(set)
	acc2vsh = dict()
	for row in self.query_validacc1():
	    vsh2acc[row['validseqhash']].add(row['accession'])
	    acc2vsh[row['accession']] = row['validseqhash']
	for acc,vsh in acc2vsh.items():
	    if len(vsh2acc[vsh]) > 1 and self.invalid(acc):
	        any = False
		for acc2 in vsh2acc[vsh]:
		    if acc == acc2:
			continue
		    if self.invalid(acc2):
			continue
		    assert any==False
		    any = True
		    yield acc,acc2
		if not any:
		    yield acc,None
		    # if acc != min(vsh2acc[vsh]):
		    #     yield acc,min(vsh2acc[vsh]),"*"

    def invalid(self,accession):
	for row in self.query_validacc1(accession=accession):
	    if row['validseqhash'] != row['seqhash']:
	        return True
	return False

    def allaccessions(self):
        for row in self.query_exists():
            yield row['accession']

    def gethash(self,accession):
	for row in self.query_hash(accession=accession):
	    return row['hash']

    def allhash(self):
	for row in self.query_hash():
	    yield row['accession'],row['hash']

    def allhashedseq(self):
	for row in self.query_hashedseq():
	    yield row['hash'],row['seq'],row['accession']

    def gethashedseq(self,hash=None,seq=None):
        # we can lookup with hash or seq
	if seq != None:
	    thehash = hashlib.sha256(seq.strip()).hexdigest().lower()
	if hash != None:
	    thehash=hash
	for row in self.query_hashedseq(hash=thehash):
	    if hash == None:
	        return row['hash'],row['accession']
	    if seq == None:
		return row['hash'],row['accession']
	return None,None

    def gettaxa(self,accession):
	for row in self.query_taxonomy(accession=accession):
	    yield row['taxon']

    def bytaxa(self,taxon):
	for row in self.query_taxonomy(taxon=taxon):
	    yield row['accession']

    def alltaxa(self):
	for row in self.query_taxonomy():
	    yield row['accession'],row['taxon']

    def gettopo(self,accession):
	for row in self.query_topology(accession=accession):
	    return row['topology']
	return None

    def alltopo(self):
	for row in self.query_topology():
	    yield row['accession'],row['topology']

    def getcomp(self,accession):
	for acc in set([self.gettopo(accession),accession]):
	    for row in self.query_composition(accession=acc):
	        return row['composition']

    def allcomp(self):
        topo = defaultdict(set)
        for s, t in self.alltopo():
            topo[t].add(s)
        seen = set()
        for row in self.query_composition():
	    t,c = row['accession'],row['composition']
            for s in topo[t]:
                if (s, c) not in seen:
                    seen.add((s, c))
                    yield s, c
                if (c, c) not in seen:
                    seen.add((c, c))
                    yield c, c

    def getbasecomp(self,accession):
	for acc in set([self.getcomp(accession),accession]):
	    for row in self.query_basecomposition(accession=acc):
	        return row['basecomposition']

    def allbasecomp(self):
	comp = defaultdict(set)
	for s, c in self.allcomp():
	    comp[c].add(s)
	seen = set()
	for row in self.query_basecomposition():
	    c,bc = row['accession'],row['basecomposition']
	    for s in comp[c]:
		if (s,bc) not in seen:
		    seen.add((s,bc))
		    yield s, bc
		if (bc,bc) not in seen:
		    seen.add((bc,bc))
		    yield bc,bc

    def _query_date_helper(self,**kwargs):
        lastacc = None; lastaccdates = set()
	for row in sorted(self.query_date(**kwargs),key=lambda r: r['accession']):
	    row['date'] = dateutil_parser.parse(row['date']).date()
	    if row['accession'] != lastacc:
		if len(lastaccdates) > 0:
		    yield lastacc,min(lastaccdates).isoformat(),max(lastaccdates).isoformat()
		lastacc = row['accession']
		lastaccdates = set([row['date']])
	    else:
		lastaccdates.add(row['date'])
	if len(lastaccdates) > 0:
	    yield lastacc,min(lastaccdates).isoformat(),max(lastaccdates).isoformat()

    def getdate(self,accession):
	for row in self._query_date_helper(accession=accession):
	    yield row[1],row[2]

    def alldate(self):
	for row in self._query_date_helper():
	    yield row

    def getimage(self,accession,notation='snfg',style='extended',format='svg'):
        assert notation in ('snfg','cfg')
        assert style in ('extended',)
        assert format in ('png','svg')
        for row in self.query_image(accession=accession,notation=notation,style=style,format=format):
            if format == "png":
                return base64.standard_b64decode(row['imagedata'])
            else:
                return row['imagedata']

    def allimage(self,notation='snfg',style='extended',format='svg'):
	assert notation in ('snfg','cfg')
	assert style in ('extended',)
	assert format in ('png','svg')
	for row in self.query_image(notation=notation,style=style,format=format):
	    if format == "png":
		yield row['accession'],base64.standard_b64decode(row['imagedata'])
	    else:
		yield row['accession'],row['imagedata']
