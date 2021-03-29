
__all__ = [ "GPTWiki", "Glycan", "Protein", "Peptide", "TransitionGroup", "Transition", "Acquisition", "Alignment" ]

import sys, re, os, os.path, json
from operator import itemgetter
import findpygly
from pygly.lockfile import FileLock, AlreadyLocked, LockTimeout

from smw import SMW

class Glycan(SMW.SMWClass):
    template = 'Glycan'
    classes = [ 'N-Glycan', 'O-Glycan' ]

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('accession')
        return kwargs.get('accession')
    
    def asclass(self,value):
        if value not in self.classes:
	    raise ValueError("Bad class value: %r"%(value,))
	return value

    def toPython(self,data):
	data = super(Glycan,self).toPython(data)

        # name is newline separated
        if isinstance(data.get('name'),basestring):
            data['name'] = map(lambda s: s.strip(),data.get('name').split('\n'))

        # topo is comma separated
        if isinstance(data.get('topo'),basestring):
            data['topo'] = map(lambda s: s.strip(),data.get('topo').split(','))

        # mw is a float
        if isinstance(data.get('mw'),basestring):
            data['mw'] = float(data.get('mw'))

        # class is comma separated, with specific possible values, sorted, so behaves as set
        if isinstance(data.get('class'),basestring):
            data['class'] = sorted(map(self.asclass,map(lambda s: s.strip(),data.get('class').split(','))))
	elif data.get('class') != None:
	    data['class'] = sorted(map(self.asclass,data.get('class')))
	
	return data

    def toTemplate(self,data):
	data = super(Glycan,self).toTemplate(data)

        if 'name' in data:
	    data['name'] = "\n".join(data['name'])

        if 'topo' in data:
            data['topo'] = ",".join(data['topo'])

        if 'mw' in data:
            data['mw'] = str(data.get('mw'))
        
	if 'class' in data:
	    data['class'] = ",".join(sorted(data['class']))

	return data
	
class Protein(SMW.SMWClass):
    template = 'Protein'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('accession')
        return kwargs.get('accession')
    
class Transition(SMW.SMWClass):
    template = 'Transition'

    @staticmethod
    def key(peptide,z1,label,z2):
	return "^".join(map(str,[peptide,z1,label,z2]))

    def toPython(self,data):
	data = super(Transition,self).toPython(data)

        for k in ('nrt','prt','rt','mz1','mz2'):
            if isinstance(data.get(k),basestring):
                data[k] = float(data.get(k))
                
        for k in ('z1','z2'):
            if isinstance(data.get(k),basestring):
                data[k] = int(data.get(k))

	return data

    def toTemplate(self,data):
	data = super(Transition,self).toTemplate(data)

        for k in ('nrt','prt','rt','mz1','mz2'):
            if k in data:
                data[k] = "%.3f"%(data[k])

        for k in ('z1','z2'):
            if k in data:
                data[k] = "%d"%(data[k])

	return data

class TransitionGroup(SMW.SMWClass):
    template = 'TransitionGroup'

    @staticmethod
    def key(peptide,z1,spectra):
	return "^".join(map(str,[peptide,z1,spectra]))

    def astrans(self,tr):
        str = tr.split(';')
        return (str[0],float(str[1]))

    def asscans(self,sc):
        str = sc.split(';')
	if len(str) > 0:
	    str[0] = int(str[0])
	if len(str) > 2:
	    str[2] = float(str[2])
        if len(str) > 3:
	    try:
	        str[3] = int(str[3])
	    except (ValueError,TypeError):
	        str[3] = None
	if len(str) > 4:
	    str[4] = float(str[4])
        return tuple(str)

    def toPython(self,data):
	data = super(TransitionGroup,self).toPython(data)

        # name is newline separated
        if isinstance(data.get('transitions'),basestring):
            data['transitions'] = map(lambda t: self.astrans(t),data.get('transitions').split(','))

        if isinstance(data.get('scans'),basestring):
            data['scans'] = map(lambda t: self.asscans(t),data.get('scans').split(','))

        for k in ('rt','prt','nrt','mz1','score','intensity','fdr'):
            if isinstance(data.get(k),basestring):
                data[k] = float(data.get(k))
                
        for k in ('z1','ntransition'):
            if isinstance(data.get(k),basestring):
                data[k] = int(data.get(k))

	if 'method' in data:
	    del data['method']

	return data

    @staticmethod
    def scan2str(scan):
        if scan[3] in (None,"None",""):
	    scan[3] = ""
	if len(scan) < 6:
	    scan = (list(scan) + [""])
        return "%d;%s;%.3f;%s;%.3f;%s"%(int(scan[0]),scan[1],float(scan[2]),scan[3],float(scan[4]),scan[5])

    def toTemplate(self,data):
	data = super(TransitionGroup,self).toTemplate(data)

        if 'transitions' in data:
            data['transitions'] = ",".join(map(lambda t: "%s;%.1f"%t,data['transitions']))

        if 'scans' in data:
            data['scans'] = ",".join(map(self.scan2str,data['scans']))

        for k in ('rt','prt','nrt','mz1','score','intensity','fdr'):
            if k in data:
                data[k] = "%.3f"%(data[k])

        for k in ('z1','ntransition'):
            if k in data:
                data[k] = "%d"%(data[k])

	return data


class Peptide(SMW.SMWClass):
    template = 'Peptide'

    @staticmethod
    def key(sequence=None,glycan=[],mod=[]):
        l = [ sequence ]
	for t in sorted(glycan,key=lambda t: int(t[1][1:])):
	    l.append(":".join(map(str,t)))
	for t in sorted(mod,key=lambda t: (int(t[1][1:]),"%+.3f"%t[0])):
	    l.append("%+.3f"%(t[0],)+":"+":".join(t[1:]))
  	return "^".join(l)

    @staticmethod
    def asglycan(g):
        sg = g.split(';')
        return tuple(sg)

    @staticmethod
    def asmod(m):
        sm = m.split(';')
        return tuple([float(sm[0])] + sm[1:])

    @staticmethod
    def asalign(a):
        sa = a.split(';')
        return (sa[0],int(sa[1]),int(sa[2]),sa[3],sa[4])

    def toPython(self,data):
	data = super(Peptide,self).toPython(data)

        # mw is a float
        if isinstance(data.get('mw'),basestring):
            data['mw'] = float(data.get('mw'))
            
        if isinstance(data.get('nrt'),basestring):
            data['nrt'] = float(data.get('nrt'))
            
        if isinstance(data.get('glycan'),basestring):
            data['glycan'] = map(self.asglycan,data.get('glycan').split(','))

        if isinstance(data.get('mod'),basestring):
            data['mod'] = map(self.asmod,data.get('mod').split(','))

        # if isinstance(data.get('alignment'),basestring):
        #     data['alignment'] = sorted(map(self.asalign,data.get('alignment').split(',')))
        #     # data['protein'] = sorted(set(map(itemgetter(0),data.get('alignment'))))

	alignments = []
	for so in data.get('_subobjs',[]):
	    if so.get('template') == 'Alignment':
		alignments.append(so)
	if len(alignments) > 0:
	    alignments = sorted(alignments,key=Alignment.sortkey)
	    data['alignments'] = alignments

        if isinstance(data.get('transgroup'),basestring):
            data['transgroup'] = sorted(data.get('transgroup').split(','))

	if 'sample' in data:
	    del data['sample']

	return data

    def toTemplate(self,data):
	data = super(Peptide,self).toTemplate(data)

        if 'mw' in data:
            data['mw'] = "%.3f"%(data.get('mw'),)

        if 'nrt' in data:
            data['nrt'] = "%.3f"%(data.get('nrt'),)

        if 'glycan' in data:
            data['glycan'] = ",".join(map(lambda t: "%s;%s"%t,sorted(data['glycan'],key=lambda t: int(t[1][1:]))))            
        if 'mod' in data:
            data['mod'] = ",".join(map(lambda t: "%+.3f;%s"%(t[0],";".join(t[1:])),
                                       sorted(data['mod'],key=lambda t: int(t[1][1:]))))

        # if 'alignment' in data:
        #     data['alignment'] = ",".join(map(lambda t: "%s;%d;%d;%s;%s"%t,sorted(data['alignment'])))
        #     # data['protein'] = ",".join(set(map(lambda t: t.split(';')[0],sorted(data['alignment'].split(',')))))

	if 'alignments' in data:
	    data['alignments'] = sorted(data['alignments'],key=Alignment.sortkey)
	    data['_subobjs'] = []
	    for al in data['alignments']:
	        data['_subobjs'].append(al)
	    del data['alignments']

        if 'transgroup' in data:
            data['transgroup'] = ",".join(data.get('transgroup'))

	return data    

class Alignment(SMW.SMWClass):
    template = 'Alignment'
    def toPython(self,data):
	data = super(Alignment,self).toPython(data)

        for k in ('start','end'):
            if isinstance(data.get(k),basestring):
                data[k] = int(data.get(k))
	
	return data

    def toTemplate(self,data):
	data = super(Alignment,self).toTemplate(data)
	return data

    @staticmethod
    def sortkey(al):
	return al.get('protein'),al.get('start')
	
class Acquisition(SMW.SMWClass):
    template = 'Acquisition'

    @staticmethod
    def pagename(**kwargs):
        assert kwargs.get('name')
        return kwargs.get('name')

    def toPython(self,data):
	data = super(Acquisition,self).toPython(data)

        for k in ('nrtslope','nrtintercept'):
            if isinstance(data.get(k),basestring):
                data[k] = float(data.get(k))

	return data

    def toTemplate(self,data):
	data = super(Acquisition,self).toTemplate(data)

        for k in ('nrtslope','nrtintercept'):
            if k in data:
		# no trailing zeros...
                data[k] = "%s"%(data[k],)

	return data

class GPTWiki(SMW.SMWSite):
    _name = 'gptwiki'

    def __init__(self,*args,**kwargs):
        super(GPTWiki,self).__init__(*args,**kwargs)
        self._cache = dict()
	self._cache_locks = dict()

    def __del__(self):
	for clsname in ("peptides","tgs","trans"):
	    if self.has_lock(clsname):
		if self._cache.get(clsname) is not None:
                    self.write_cache(clsname)
		self.release_lock(clsname)

    def cache_file(self,clsname):
	return ".%s.%s.cache"%(self.prefix,clsname,)

    def has_lock(self,clsname):
	if clsname in self._cache_locks and self._cache_locks[clsname].i_am_locking():
	    return True
        return False

    def lock_cache(self,clsname):
	if self.has_lock(clsname):
	    return
	filename = self.cache_file(clsname)
	lock = FileLock(filename)
	try:
	    lock.acquire(10)
	except LockTimeout:
	    if self._cache_locks.get(clsname):
		self._cache_locks[clsname].release()
	    raise RuntimeError("Can't obtain lock for %s cache file"%(clsname,))
	self._cache_locks[clsname] = lock

    def release_cache(self,clsname):
	if clsname not in self._cache_locks or not self._cache_locks[clsname].i_am_locking():
	    return
	self._cache_locks[clsname].release()
	del self._cache_locks[clsname]

    def read_cache(self,clsname):
	filename = self.cache_file(clsname)
	if os.path.exists(filename):
	    self._cache[clsname] = json.loads(open(filename).read())
	else:
	    if clsname in self._cache:
	        del self._cache[clsname]

    def write_cache(self,clsname):
	filename = self.cache_file(clsname)
	wh = open(filename,'w')
	wh.write(json.dumps(self._cache[clsname]))
	wh.close()

    def has_cache(self,clsname):
	return self._cache.get(clsname) != None

    def cache_ids(self,clsname):
	return self._cache.get(clsname,{}).values()

    def verify_ids(self,clsname,iterids):
	return set(self.cache_ids(clsname)) == set(iterids)

    def add_to_cache(self,clsname,key,value):
	self.lock_cache(clsname)
	if self._cache.get(clsname) is None:
	    self._cache[clsname] = {}
	self._cache[clsname][key] = value

    def get_from_cache(self,clsname,key):
	if self._cache.get(clsname) is not None:
	    return self._cache[clsname].get(key)

    def is_in_cache(self,clsname,key):
	return (self.get_from_cache(clsname,key) is not None)

    def peptideindex(self):
	self.lock_cache("peptides")
        self.read_cache("peptides")
	if self.has_cache("peptides"):
	    if self.verify_ids("peptides",self.iterpeptideids()):
		self.release_cache("peptides")
		return
	print >>sys.stderr, "(Re-)computing peptide index from site..."
        for p in self.iterpeptides():
            key = Peptide.key(p.get('sequence'),p.get('glycan',[]),p.get('mod',[]))
	    self.add_to_cache("peptides",key,p.get('id'))
	self.write_cache("peptides")
	self.release_cache("peptides")
	print >>sys.stderr, "(Re-)computing peptide index from site... done."

    def tgindex(self):
	self.lock_cache("tgs")
	self.read_cache("tgs")
	if self.has_cache("tgs"):
	    if self.verify_ids("tgs",self.itertransgroupids()):
		self.release_cache("tgs")
		return
	print >>sys.stderr, "(Re-)computing transition group index..."
        for tg in self.itertransgroups():
            key = TransitionGroup.key(tg.get('peptide'),tg.get('z1'),tg.get('spectra'))
	    self.add_to_cache("tgs",key,tg.get('id'))
	self.write_cache("tgs")
	self.release_cache("tgs")
	print >>sys.stderr, "(Re-)computing transition group index... done."

    def transindex(self):
	self.lock_cache("trans")
	self.read_cache("trans")
	if self.has_cache("trans"):
	    if self.verify_ids("trans",self.itertransitionids()):
		self.release_cache("trans")
		return
	print >>sys.stderr, "(Re-)computing transition index..."
        for t in self.itertransitions():
            key = Transition.key(t.get('peptide'),t.get('z1'),t.get('label'),t.get('z2'))
	    self.add_to_cache("trans",key,t.get('id'))
	self.write_cache("trans")
	self.release_cache("trans")
	print >>sys.stderr, "(Re-)computing transition index... done."

    def addprotein(self,accession,**kw):
        p = self.get(accession)
        if p:
            p.update(kw)
        else:
            p = Protein(accession=accession,**kw)
        self.put(p)
        
    def addglycan(self,accession,**kw):
        p = self.get(accession)
        if p:
            p.update(kw)
        else:
            p = Glycan(accession=accession,**kw)
        self.put(p)

    def addacquisition(self,name,**kw):
	a = self.get(name)
	if a:
	    a.update(**kw)
	else:
	    a = Acquisition(name=name,**kw)
	self.put(a)

    def newtransid(self):
        maxid = 0
        for i,id in enumerate(sorted(map(lambda id: int(id[2:]),self.cache_ids("trans")))):
            if (i+1) < id:
                return "TR%06d"%(i+1,)
            maxid = id
        return "TR%06d"%(maxid+1,)

    def addtransition(self,peptide,z1,label,z2,**kw):
	if not self.has_cache("trans"):
	    self.transindex()
	key = Transition.key(peptide,z1,label,z2)
        id = self.get_from_cache("trans",key)
        if id:
            t = self.get(id)
            t.update(**kw)
        else:
            id = self.newtransid()
            t = Transition(id=id,peptide=peptide,z1=z1,label=label,z2=z2,**kw)
	    self.add_to_cache("trans",key,id)
        return t,self.put(t)

    def newtransgrpid(self):
        maxid = 0
        for i,id in enumerate(sorted(map(lambda id: int(id[2:]),self.cache_ids("tgs")))):
            if (i+1) < id:
                return "TG%06d"%(i+1,)
            maxid = id
        return "TG%06d"%(maxid+1,)

    def addtransgroup(self,peptide,z1,spectra,**kw):
	if not self.has_cache("tgs"):
	    self.tgindex()
	key = TransitionGroup.key(peptide,z1,spectra)
        id = self.get_from_cache("tgs",key)
        if id:
            tg = self.get(id)
            tg.update(**kw)
        else:
            id = self.newtransgrpid()
            tg = TransitionGroup(id=id,peptide=peptide,z1=z1,spectra=spectra,**kw)
	    self.add_to_cache("tgs",key,id)
        return tg,self.put(tg)

    # def iterpeptides(self):
    #	for p in map(self.get,self.itercat('Peptide')):
    #        yield p

    def iterproteins(self):
	for p in map(self.get,self.itercat('Protein')):
            yield p

    def iterprefix(self,prefix):
	for pagename in self.site.allpages(prefix=prefix,generator=False):
	    yield self.get(pagename)

    def iteridbyregex(self,regex,prefix=None):
	regex = re.compile(regex)
        for pagename in self.site.allpages(prefix=prefix,generator=False):
            if regex.search(pagename):
		yield pagename

    def iterregex(self,regex,prefix=None):
	regex = re.compile(regex)
        for page in self.site.allpages(prefix=prefix,generator=True):
            if regex.search(page.name):
                yield self.parse_pagetext(page.name,page.text())

    def iterglycans(self):
	for g in self.iterregex(r'^G\d{5}[A-Z]{2}$',prefix="G"):
	    yield g

    def itertransgroupids(self):
	for p in self.iteridbyregex(r'^TG\d{6}$',prefix="TG"):
	    yield p

    def itertransgroups(self):
	for tg in self.iterregex(r'^TG\d{6}$',prefix="TG"):
	    yield tg

    def itertransitionids(self):
	for p in self.iteridbyregex(r'^TR\d{6}$',prefix="TR"):
	    yield p

    def itertransitions(self):
	for tg in self.iterregex(r'^TR\d{6}$',prefix="TR"):
	    yield tg

    def iterask(self,query):
	for it in self.site.ask(query):
	    if isinstance(it,basestring):
	        yield self.get(it)
	    else:
	        yield self.get(it['fulltext'])

    def iterspec(self,**kw):
	query = """
	[[Category:Acquisition]]
        """.strip()
	if kw.get('method'):
	    query += "\n[[gptwiki:formethod::%(method)s]]"%kw
	if kw.get('sample'):
	    query += "\n[[gptwiki:forsample::%(sample)s]]"%kw
	if kw.get('anfrac'):
	    query += "\n[[gptwiki:analfraction::%(anfrac)s]]"%kw
	if kw.get('acqtype'):
	    query += "\n[[gptwiki:type::%(acqtype)s]]"%kw
	if kw.get('inst'):
	    query += "\n[[gptwiki:spectra.gptwiki:instrument::%(inst)s]]"%kw
	for sp in self.iterask(query):
	    yield sp

    def itertgs(self,**kw):
	query = """
	[[Category:TransitionGroup]]
        """
	if kw.get('spectra'):
	    query += "\n[[gptwiki:spectra::%(spectra)s]]"%kw
	if kw.get('method'):
	    query += "\n[[gptwiki:spectra.gptwiki:formethod::%(method)s]]"%kw
	if kw.get('sample'):
	    query += "\n[[gptwiki:spectra.gptwiki:forsample::%(sample)s]]"%kw
	if kw.get('anfrac'):
	    query += "\n[[gptwiki:spectra.gptwiki:analfraction::%(anfrac)s]]"%kw
	if kw.get('acqtype'):
	    query += "\n[[gptwiki:spectra.gptwiki:type::%(acqtype)s]]"%kw
	if kw.get('inst'):
	    query += "\n[[gptwiki:spectra.gptwiki:instrument::%(inst)s]]"%kw
	for tg in self.iterask(query):
	    yield tg

    def iterpep(self,**kw):
	query = """
	[[Category:Peptide]]
        """
	for pep in self.iterask(query):
	    yield pep

    def iterpeptides(self):
	for p in self.iterregex(r'^PE\d{6}$',prefix="PE"):
	    yield p

    def iterpeptideids(self):
	for p in self.iteridbyregex(r'^PE\d{6}$',prefix="PE"):
	    yield p

    def newpeptideid(self):
        maxid = 0
        for i,id in enumerate(sorted(map(lambda id: int(id[2:]),self.cache_ids("peptides")))):
            if (i+1) < id:
                return "PE%06d"%(i+1,)
            maxid = id
        return "PE%06d"%(maxid+1,)

    def findpeptide(self,sequence,glycans=[],mods=[]):
	if not self.has_cache("peptides"):
	    self.peptideindex()
        key = Peptide.key(sequence,glycan=glycans,mod=mods)
        return key,self.get_from_cache("peptides",key)

    def addpeptide(self,sequence,glycans=[],mods=[],**kw):
	if not self.has_cache("peptides"):
	    self.peptideindex()

        if kw.get('id'):
            key = Peptide.key(sequence,glycan=glycans,mod=mods)
            p = Peptide(sequence=sequence,glycan=glycans,mod=mods,**kw)
	    assert self.get_from_cache("peptides",key) in (None,kw.get('id'))
	    if self.get_from_cache("peptides",key) is None:
		self.add_to_cache("peptides",key,kw.get('id'))
            return p,self.put(p)
        
        key,id = self.findpeptide(sequence,glycans=glycans,mods=mods)
        if id:
            p = self.get(id)
            p.update(**kw)
        else:
            id = self.newpeptideid()
            p = Peptide(id=id,sequence=sequence,glycan=glycans,mod=mods,**kw)
	    self.add_to_cache("peptides",key,id)
	return p,self.put(p)

