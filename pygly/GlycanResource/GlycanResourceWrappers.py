
# Defaults ensure a 10-way partition of GlyTouCan accessions
def partitioner(kwarg="accession",fmt="G%%0%dd.*",digits=1,values='decimal'):
    fmtstr = fmt%(digits,)
    def partition(fn):
        def wrapper(self,*args,**kw):
            if kwarg not in kw:
		if values == "decimal":
                    for i in range(0,10**digits):
			kw[kwarg] = fmtstr%(i,)
                        for row in fn(self,*args,**kw):
                            yield row
		elif values == "hexidecimal":
		    for i in range(0,16**digits):
			kw[kwarg] = fmtstr%(i,)
                        for row in fn(self,*args,**kw):
                            yield row
            else:
                for row in fn(self,*args,**kw):
                    yield row
        return wrapper
    return partition

# Defaults presume acccession is the (only) kwarg to prefetch for...
def prefetcher(kwarg='accession',**kw):
    usecache = kw.get('usecache',False)
    def prefetch(fn):
        def wrapper(self,**kw):
            # print >>sys.stderr, kw
            kw1 = dict((k,v) for k,v in kw.items() if k != kwarg)
            key = fn.__name__+":"+":".join("%s=%s"%(k,v) for k,v in sorted(kw1.items()))
            # print >>sys.stderr, "cache key:",key
            if key not in self._cache:
		if not usecache or not self._cacheondisk.has_key(key):
                    # print >>sys.stderr, "fill cache:",key
		    self._cache[key] = {}
                    if usecache:
                        self._cachedirty[key] = True
                    for row in fn(self,**kw1):
                        if row[kwarg] not in self._cache[key]:
                            self._cache[key][row[kwarg]] = []
                        self._cache[key][row[kwarg]].append(row)
                    if usecache:
                        self.writecache()
		else:
		    # print >>sys.stderr, key
                    self._cache[key] = self._cacheondisk[key]
		    self._cachedirty[key] = False

            if kwarg in kw:
                for row in self._cache[key].get(kw[kwarg],[]):
                    yield row
	    else:
                for acc in self._cache[key]:
                    for row in self._cache[key][acc]:
                        yield row
        return wrapper
    return prefetch
