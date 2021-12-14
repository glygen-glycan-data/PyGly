#!/bin/env python2

import sys,traceback,time

from getwiki import GlycoMotifWiki, Publication
w = GlycoMotifWiki()

refs = []
if len(sys.argv) > 1:
  h = open(sys.argv[1])
else:
  h = sys.stdin
for l in h:
    sl = l.split()
    try:
        refid = int(sl[0])
    except:
        refid = sl[0]
    if refid == "-":
	refid = None
    pmid = sl[1]
    refs.append((refid,pmid))

currentpmids = set()
for pub in w.iterpages(include_categories=['Publication']):
    pmid = w.get(pub.name).get('pmid')
    if pmid != None: 
	currentpmids.add(pmid)

from Bio import Entrez
Entrez.email = "nje5@georgetown.edu"

for ref,pmid in refs:

    if pmid in currentpmids:
	continue

    handle = Entrez.efetch(db="pubmed", id=pmid, retmode='xml')
    record = Entrez.read(handle)                                                                                         
    handle.close()                                                                                                       
    for article in record['PubmedArticle']:                                                                              
        pmid = None                                                                                                      
        doi = None                                                                                                       
        for artid in article['PubmedData']['ArticleIdList']:                                                             
            if artid.attributes["IdType"] == "pubmed":                                                                   
                pmid = str(artid)                                                                                        
            if artid.attributes["IdType"] == "doi":                                                                      
                doi = str(artid)                                                                                         
        theart = article['MedlineCitation']['Article']                                                                   
        vol = theart["Journal"]["JournalIssue"]["Volume"]                                                                
        issue = theart["Journal"]["JournalIssue"].get("Issue")                                                           
        pubdate = theart["Journal"]["JournalIssue"]["PubDate"]
	if "Year" in pubdate:
	    year = pubdate["Year"]
	elif "MedlineDate" in pubdate:
            year = pubdate["MedlineDate"].split(None,1)[0]
	else:
	    print theart["Journal"]["JournalIssue"]
	    raise RuntimeError("Can't figure out the year")
        journal = theart["Journal"]["ISOAbbreviation"].rstrip('.')
        authors = []                                                                                                     
        for au in theart['AuthorList']:                                                                                  
            inits = au['Initials']                                                                                       
            lastname = au['LastName']                                                                                    
            authors.append("%s %s"%(lastname, inits))                                                                     
        authors = ", ".join(authors)                                                                                     
        title = theart['ArticleTitle'].rstrip('.')
        page = theart['Pagination']['MedlinePgn']                                                                        
	if issue:
	    volume = "%s(%s)"%(vol,issue)
	else:
	    volume = vol

    if isinstance(ref,int):
	citedby = ('Cummings_(2009)_The_repertoire_of_glycan_determinants_in_the_human_glycome',ref)
    elif isinstance(ref,basestring) and ref.startswith('EP'):
	citedby = ('GE',ref)
    else:
	citedby = None

    pub = Publication(authors=authors,title=title,journal=journal,year=year,volume=volume,pages=page,
                      pmid=pmid,doi=doi,citedby=citedby)
    print Publication.pagename(**pub.data)
    print pub.astemplate()
    if w.put(pub):
        print >>sys.stderr, pmid

    time.sleep(5)
