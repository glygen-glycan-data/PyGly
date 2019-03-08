
def memoize(key=None):
    thecache = {}
    def wrap(f):        
        def decorated_function(*args,**kw):
            if key != None:
                _key = key(args,kw)
            else:
                _key = args
	    _key = (f,args)
            if _key not in thecache:
                thecache[_key] = f(*args,**kw)
            return thecache[_key]
        return decorated_function                                                                                 
    return wrap
