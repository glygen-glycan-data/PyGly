
import sqlite3, re
import cPickle as pickle
import zlib, os, os.path, sys, math

from GlyMWFilter import GlyMWFilter
from GlyCompFilter import GlyCompFilter
from GlyNLinkedFilter import GlyNLinkedFilter
from GlyLactosamineFilter import GlyLactosamineFilter
from GlyHighMannoseFilter import GlyHighMannoseFilter
from GlyOxfordFilter import GlyOxfordFilter
from GlyLinCodeFilter import GlyLinCodeFilter
from GlyIUPACFilter import GlyIUPACFilter
from CompositionTable import ResidueCompositionTable, PermethylCompositionTable
from ElementMass import MonoisotopicElementMass
from GlyRecord import GlyRecord
from subtree import linearcodeSubtreeEquals, iupacSubtreeEquals

class GlyDbIndex:
    createTable = """
        create table theindex (
          id integer primary key autoincrement, 
          source varchar not null,
          accession varchar not null,
          mw float not null,
	  intmw integer not null,
          pmw float not null,
	  intpmw integer not null,
          nlinked bool not null,
          lactosamine bool not null,
          highmannose bool not null,
	  oxford varchar,
          lincode bool not null,
	  toporepr int,
          topogrp varchar,
	  lincoderepr int,
          Hex int not null,
          HexNAc int not null,
          NeuAc int not null,
          NeuGc int not null,
          Fuc int not null,
          Xyl int not null,
          Xxx int not null,
          data blob not null
        );
    """
    createTable1 = """
       create table accindex (
          glyid integer, 
          source varchar not null,
          accession varchar not null
        );
    """
    valid = set(['Hex','HexNAc','NeuAc','NeuGc','Fuc','Xyl'])
    insert = """
        insert into theindex values (
          ?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?
        )
    """
    insert1 = """
        insert into accindex values (
          ?,?,?
        )
    """
    select = """
        select * from theindex
        where mw >= ? and mw <= ?
    """
    select1 = """
        select * from theindex
        order by abs(mw-?) asc limit 1
    """
    select2 = """
        select id,data from theindex
        order by mw asc
    """
    index = """
        create index mwindex on theindex (mw,nlinked)
    """
    index1 = """
        create index nlindex on theindex (nlinked,mw)
    """
    index2 = """
        create unique index theaccindex on theindex (accession)
    """
    index3 = """
        create index pmwindex on theindex (pmw,nlinked)
    """
    index4 = """
        create index theaccindex1 on accindex (accession)
    """
    index5 = """
        create index glyidindex on accindex (glyid)
    """
    index6 = """
	create index index6 on theindex (toporepr,nlinked,intmw)
    """
    index7 = """
	create index index7 on theindex (toporepr,nlinked,intpmw)
    """
    minsel = """
        select min(mw) from theindex
    """
    maxsel = """
        select max(mw) from theindex
    """
    minpmwsel = """
        select min(pmw) from theindex where pmw > 0
    """
    maxpmwsel = """
        select max(pmw) from theindex
    """
    cntsel = """
        select count(*) from theindex
    """

    def __init__(self,glydb=None,force=False,filters=None,indexfile=None):
	if isinstance(glydb, basestring):
	    if glydb.endswith('.index'):
	        self.indexfile = glydb
	    else:
	        self.indexfile = glydb+'.index'
            self.glydb = None
	else:
            self.glydb = glydb
	    if not indexfile:
                self.indexfile = glydb.filename+'.index'
	    else:
                self.indexfile = indexfile
        self.conn = None
        self.filters = []
        if filters:
            self.filters.extend(filters)
	if self.glydb:
            self.build(force=force)

    def clean(self):
        try:
            os.unlink(self.indexfile)
        except OSError:
            pass
        try:
            os.unlink(self.indexfile+"-journal")
        except OSError:
            pass

    def exists(self):
	if not os.path.exists(self.indexfile):
	    return False
	if not os.path.exists(self.indexfile[:-6]):
	    return True
	if os.path.getmtime(self.indexfile) < \
	   os.path.getmtime(self.indexfile[:-6]):
	    return False
        # Truncated?
        if os.path.getsize(self.indexfile) < os.path.getsize(self.indexfile[:-6]):
            return False
        return True

    def build(self,force=False):
        if self.exists() and not force:
            return
        glyind = 0;
        self.clean()
        self.conn = sqlite3.connect(self.indexfile,isolation_level='DEFERRED')
        self.conn.execute(self.createTable)
        self.conn.execute(self.createTable1)
        self.conn.execute(self.index)
        self.conn.execute(self.index1)
        self.conn.execute(self.index2)
        self.conn.execute(self.index3)
        self.conn.execute(self.index4)
        self.conn.execute(self.index5)
	gdb = self.glydb
        for f in self.filters:
            gdb = f(gdb)
        gdb = GlyMWFilter(gdb,
                          ResidueCompositionTable(),
                          PermethylCompositionTable(),
                          MonoisotopicElementMass())
        gdb = GlyCompFilter(gdb)
        gdb = GlyNLinkedFilter(gdb)
        gdb = GlyLactosamineFilter(gdb)
        gdb = GlyHighMannoseFilter(gdb)
	gdb = GlyOxfordFilter(gdb)
        gdb = GlyLinCodeFilter(gdb)
        gdb = GlyIUPACFilter(gdb)
        for gr in gdb:
            if 'molecular_weight' not in gr:
                continue
            if 'composition' not in gr:
                continue
	    if 'high-mannose' not in gr:
		continue
            cmp = gr['composition']
            remaining_residues = sum(v for k,v in cmp.iteritems()
                                     if k not in GlyDbIndex.valid)
            glyind += 1
            self.conn.execute(self.insert,(glyind,
                                           gr.source,
                                           gr.accession,
                                           gr['molecular_weight'],
					   int(gr['molecular_weight']),
                                           gr.get('permethylated_molecular_weight',0.0),
                                           int(gr.get('permethylated_molecular_weight',0.0)),
                                           gr['mini-nlinked'],
                                           gr['lactosamine'],
                                           gr['high-mannose'],
					   gr.get('oxford',""),
                                           gr['lincode'],
					   None,
					   None,
                                           None,
                                           cmp['Hex'],
                                           cmp['HexNAc'],
                                           cmp['NeuAc'],
                                           cmp['NeuGc'],
                                           cmp['Fuc'],
                                           cmp['Xyl'],
                                           remaining_residues,
                                           self.tostr(gr)))
            self.conn.execute(self.insert1,(glyind,gr.source,gr.accession))
            for res in gr.get('resource',[]):
                source,accession = res.split(':')
                self.conn.execute(self.insert1,(glyind,source,accession))
            for res in gr.get('taxon',[]):
                source,accession = res.split(':')
                self.conn.execute(self.insert1,(glyind,source+"taxa",accession))
        self.conn.commit()
	lastmwstr = None
	mwind = 0
	group = []
        for r in list(self.getallbymw()):
	    mwstr = "%.3f"%r['molecular_weight']
	    if mwstr != lastmwstr:
                if len(group) != 0:
                    self.processgroup(self.conn,mwind,group)
		mwind += 1
                lastmwstr = mwstr
                group = []
	    group.append(r)
        if len(group) != 0:
            self.processgroup(self.conn,mwind,group)
        self.conn.commit()
        self.conn = None

    @staticmethod
    def gdbsortkey(acc):
	if '-' in acc:
	    return (1,int(acc.split('-')[0][3:]))
	else:
	    return (0,int(acc[3:]))

    @staticmethod
    def processgroup(conn,grpind,group):
        try:
            group.sort(key=lambda gr: GlyDbIndex.gdbsortkey(gr.accession))
        except:
            group.sort(key=lambda gr: gr.accession)
        # print "++++++++++++++++++++++++++"
        # print grpind,map(lambda gr: gr.accession,group)
        # for gr in group:
        #     try:
        #         print gr
        #     except:
        #         pass
        glycmp = iupacSubtreeEquals()
        lc = {}; lcre = {}
        for i,gr in enumerate(group):
            if 'LinearCode' in gr:
                lc[i] = gr['LinearCode']
                lcre[i] = re.compile(lc[i].replace('?','.').replace('(','\(').replace(')','\)'))
        lcrepr = {}
        toporepr = {}
        for i1,m1 in enumerate(group):
            lcrepr[i1] = set([i1])
        toremove = set()
        for i1,m1 in enumerate(group):
            if i1 not in lc:
                continue
            lckeep = True
            for i2,m2 in enumerate(group):
                if i2 not in lc:
                    continue
                if i1 == i2:
                    continue
                if lc[i1] == lc[i2]:
                    if i2 < i1:
                        lckeep = False
                        lcrepr[i2].add(i1)
                    continue
                if lcre[i1].search(lc[i2]):
                    lckeep = False
                    lcrepr[i2].add(i1)
                    continue
            if not lckeep:
                toremove.add(i1)
        for i1 in toremove:
            del lcrepr[i1]
        for i1 in lcrepr.keys():
            toporepr[i1] = set(lcrepr[i1])
        toremove = set()
        for i1 in lcrepr.keys():
            m1 = group[i1]
            topokeep = True
            for i2 in lcrepr.keys():
                m2 = group[i2]
                if i1 > i2 and glycmp.compare(m1.glycan,m2.glycan):
                    topokeep = False
                    toporepr[i2].update(toporepr[i1])
            if not topokeep:
                toremove.add(i1)
        for i1 in toremove:
            del toporepr[i1]
        # print lcrepr
        # print toporepr
        lcreprmap = {}
        for repr,others in lcrepr.items():
            for o in others:
                lcreprmap[group[o]['id']] = group[repr]['id']
            lcreprmap[group[repr]['id']] = 0
        toporeprmap = {}
        topogrp = {}
        i = 1
        for repr,others in toporepr.items():
            for o in others:
                toporeprmap[group[o]['id']] = group[repr]['id']
                topogrp[group[o]['id']] = "%d.%d"%(grpind,i)
            toporeprmap[group[repr]['id']] = 0
            i += 1
        # print lcreprmap
        # print toporeprmap
        # print topogrp
        assert len(toporeprmap) == len(group)
        assert len(lcreprmap) == len(group)
        assert len(topogrp) == len(group)
        for id in topogrp:
            conn.execute("update theindex set toporepr = ?, lincoderepr = ?, topogrp = ? where id = ?",
                         (toporeprmap[id],lcreprmap[id],topogrp[id],id))
                
    @staticmethod
    def tostr1(o):
        return sqlite3.Binary(zlib.compress(o))
    @staticmethod
    def toobj1(s):
        return zlib.decompress(s)
    @staticmethod
    def tostr(o):
        return sqlite3.Binary(zlib.compress(pickle.dumps(o)))
    @staticmethod
    def toobj(s):
        return pickle.loads(zlib.decompress(s))
    def get(self,**kw):
        if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
        c = self.conn.cursor()
	count = False
	if 'count' in kw:
	    if kw['count']:
                select = "select count(*) from theindex"
		count = True
	    else:
                select = "select data from theindex"
	    del kw['count']
	else:
            select = "select data from theindex"
        values = []
        for i,(k,v) in enumerate(kw.items()):
            if i != 0:
                select += " and "
            else:
                select += " where "
            m = re.search(r'((min|max)?)(.*)$',k)
            sqlcol = m.group(3)
            minmaxeq = m.group(1)
            if minmaxeq == "min":
                select += "%s >= ?\n"%(sqlcol,)
            elif minmaxeq == "max":
                select += "%s <= ?\n"%(sqlcol,)
            elif minmaxeq == "":
                select += "%s == ?\n"%(sqlcol,)
            values.append(v)
	if 'minmw' in kw and 'maxmw' in kw:
	    select += 'and intmw in (%s)'%(','.join(map(str,range(int(math.floor(kw['minmw'])),int(math.ceil(kw['maxmw']))+1))),)
	if 'minpmw' in kw and 'maxpmw' in kw:
	    select += 'and intpmw in (%s)'%(','.join(map(str,range(int(math.floor(kw['minpmw'])),int(math.ceil(kw['maxpmw']))+1))),)
	if count:
	    yield c.execute(select,values).next()[0]
	    return
        for r in c.execute(select,values):
            yield self.toobj(r[0])
        c.close()
    def count(self,**kwargs):
	kwargs['count'] = True
	return self.get(**kwargs).next()
    def __iter__(self):
	return self.getall()
    def getallbymw(self):
	if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
        c = self.conn.cursor()
        for r in c.execute(self.select2):
            o = self.toobj(r[1])
            o['id'] = r[0]
	    yield o
    def getall(self):
        for r in self.get():
            yield r
    def minmw(self):
        if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
        c = self.conn.cursor()
	return c.execute(self.minsel).next()[0]
    def maxmw(self):
        if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
        c = self.conn.cursor()
	return c.execute(self.maxsel).next()[0]
    def minpmw(self):
        if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
        c = self.conn.cursor()
	return c.execute(self.minpmwsel).next()[0]
    def maxpmw(self):
        if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
        c = self.conn.cursor()
	return c.execute(self.maxpmwsel).next()[0]
    def toporepr(self,acc):
	if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
	c = self.conn.cursor()
	toporepr = c.execute("select toporepr from theindex where accession = ?",(acc,)).next()[0]
	for r in self.get(id=toporepr):
            return r
	return None
    def topoequiv(self,acc):
	if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
	c = self.conn.cursor()
	toporepr = c.execute("select id from theindex where accession = ?",(acc,)).next()[0]
	for r in self.get(toporepr=toporepr):
            yield r
    def lincodeequiv(self,acc):
	if not self.conn:
            self.conn = sqlite3.connect(self.indexfile,isolation_level=None)
	c = self.conn.cursor()
	lincoderepr = c.execute("select id from theindex where accession = ?",(acc,)).next()[0]
	for r in self.get(lincoderepr=lincoderepr):
            yield r
    def allaccessions(self,acc,mergelc=False,mergetopo=False):
	allacc = set([acc])
	if mergetopo:
	    for gr in self.topoequiv(acc):
		allacc.add(gr.accession)
	if mergelc:
	    for gr in self.lincodeequiv(acc):
		allacc.add(gr.accession)
	return allacc

if __name__ == '__main__':
    from GlycoCTDatabase import GlycoCTDatabase, GlycomeDBDatabase
    gdb = GlycoCTDatabase(sys.argv[1])
    index = GlyDbIndex(gdb,force=True)
    print index.count()

	
	
