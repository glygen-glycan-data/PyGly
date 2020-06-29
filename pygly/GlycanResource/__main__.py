#!/bin/env python27

import sys, re, os.path

sys.path.append(os.path.dirname(os.path.realpath(os.path.join(os.path.dirname(__file__),".."))))

from pygly.GlycanResource import *

def tostr(s):
    try:
	return str(s)
    except ValueError:
	pass
    try:
	return s.value
    except ValueError:
	pass
    return s

cls = sys.argv[1]
resource = eval(cls+"()")
query = sys.argv[2]
method = getattr(resource,query)
headers = None
args = sys.argv[3:]
kwargs = {}
for i in range(len(args)-1,-1,-1):
    if re.search(r'^[a-z]+=',args[i]):
        k,v = args[i].split('=',1)
        kwargs[k] = v
        del args[i]
result = method(*args,**kwargs)
if isinstance(result,basestring) or not hasattr(result,'next'):
    if not isinstance(result,list) and not isinstance(result,tuple):
	result = [ result ]
    print "\t".join(map(tostr,result))
else:
    for r in result:
        if isinstance(r,basestring):
            print r
        elif isinstance(r,dict):
            if headers == None:
                headers = sorted(r.keys())
                if 'accession' in headers:
                    headers.remove('accession')
                    headers = ['accession'] + headers
                print "\t".join(headers)
            print "\t".join(map(tostr,map(r.get,headers)))
        else:
	    # print r
            # print map(tostr,r)
            print "\t".join(map(tostr,r))

