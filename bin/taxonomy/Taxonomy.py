
import tempfile, os, os.path, sys, gzip, sqlite3, urllib, string, zipfile
from operator import itemgetter

class NCBITaxonomy(object):
    createTable1 = """
        create table names (
            tax_id INT,
            name TEXT,
            name_class TEXT,
            type INT
        );
    """
    createTable2 = """
        create table nodes (
            tax_id INTEGER PRIMARY KEY,
            rank TEXT,
            parent_id INT
        );
    """
    createTable3 = """
        create table nucgimap (
            tax_id INT,
            gi INT
        );
    """
    createTable4 = """
        create table protgimap (
            tax_id INT,
            gi INT
        );
    """
    createTable5 = """
        create table merged (
            old_id INT,
            tax_id INT
        );
    """
    indexes = """
        create index index1 on names (tax_id,type);
        create index index2 on names (name);
        create index index3 on names (name_class, name);
        create index index4 on nucgimap (gi);
        create index index5 on protgimap (gi);
        create index index6 on nucgimap (tax_id);
        create index index7 on protgimap (tax_id);
        create index index8 on merged (old_id);
    """
    insert1 = """
        insert into nodes values (?,?,?);
    """
    insert2 = """
        insert into names values (?,?,?,?);
    """
    insert3 = """
        insert into nucgimap values (?,?);
    """
    insert4 = """
        insert into protgimap values (?,?);
    """
    insert5 = """
        insert into merged values (?,?);
    """
    select1 = """
        select * from nodes
        where tax_id = ? limit 1
    """    
    select2 = """
        select * from merged
        where old_id = ? limit 1
    """    
    select3 = """
        select * from names
        where name = ?
    """
    select4 = """
        select * from names
        where tax_id = ? and type = 1 and name_class = 'scientific name' limit 1
    """
    select5 = """
        select * from names
        where tax_id = ? and type = 1
    """
    select6 = """
        select * from nodes
        where parent_id = ?
    """
    select7 = """
        select * from nucgimap
        where gi = ? limit 1
    """
    select8 = """
        select * from protgimap
        where gi = ? limit 1
    """
    select9 = """
	select distinct tax_id from names
	where name LIKE ?        
    """
    basefilename = 'ncbitaxonomy'
    def __init__(self,datadir=None):
        pkgdir = __file__
        while not os.path.isdir(pkgdir):
            pkgdir,dummy = os.path.split(pkgdir)
        self.data = pkgdir
	if datadir:
	    self.data = datadir
        self.filename = self.basefilename+".sqlite3"
        self.filename = os.path.join(pkgdir,self.filename)
        if self.exists():
            self.init()

    def exists(self):
	return os.path.exists(self.filename)

    def build(self,force=False):
        return NCBITaxonomy._build(self,force=force,totaxid=True)

    def _build(self,force=False,totaxid=True):
	if self.exists() and not force:
	    return
	self.clean()
        conn = sqlite3.connect(self.filename,
                               isolation_level='EXCLUSIVE')
	conn.text_factory = str
	try:
            self.parse(self.data,conn,verbose=True,totaxid=totaxid)
	except:
	    self.clean()
	    raise

    def clean(self):
        try:
            os.unlink(self.filename)
        except OSError:
            pass
        try:
            os.unlink(self.filename+"-journal")
        except OSError:
            pass

    def parse(self,datadir,conn,verbose=False,totaxid=True):
        if verbose:
            print >>sys.stderr, "Progress: Creating tables"        

        conn.execute(self.createTable1)
        conn.execute(self.createTable2)
        if totaxid:
            conn.execute(self.createTable3)
            conn.execute(self.createTable4)
        conn.execute(self.createTable5)

        if verbose:
            print >>sys.stderr, "Progress: Loading nodes"        

        zf = zipfile.ZipFile(os.path.join(datadir,'taxdmp.zip'),'r')
        for l in zf.read('nodes.dmp').splitlines():
            l = l.rstrip('\t|\n')
            sl = l.split('\t|\t',3)
            taxid = int(sl[0])
            partaxid = int(sl[1])
            rank = sl[2].strip()
            if rank == 'no rank':
                rank = ""
            conn.execute(self.insert1,(taxid,rank,partaxid))

        if verbose:
            print >>sys.stderr, "Progress: Loading names"        

        for i,l in enumerate(zf.read('names.dmp').splitlines()):
            l = l.rstrip('\t|\n')
            sl = map(str.strip,l.split('\t|\t'))
            taxid = int(sl[0])
            conn.execute(self.insert2,(taxid,sl[1],sl[3],1));
            if self.normalize(sl[1]) != sl[1]:
                conn.execute(self.insert2,(taxid,self.normalize(sl[1]),sl[3] + ' normalized',2));
            if sl[2]:
                conn.execute(self.insert2,(taxid,sl[2],sl[3] + ' unique',3));
                if self.normalize(sl[2]) != sl[2]:
                    conn.execute(self.insert2,(taxid,self.normalize(sl[2]),sl[3] + ' normalized unique',4));
            if verbose and i % 10000 == 0 and i > 0:
                print >>sys.stderr, "Inserted record %d"%i
            
        if verbose:
            print >>sys.stderr, "Progress: Loading merged taxids"        

        for l in zf.read('merged.dmp').splitlines():
            l = l.rstrip('\t|\n')
            sl = map(int,l.split('\t|\t'))
            conn.execute(self.insert5,(sl[0],sl[1]));

        zf.close()

        if totaxid:
            if verbose:
                print >>sys.stderr, "Progress: Loading nucleotide gimap"        
            
            for i,l in enumerate(gzip.open(os.path.join(datadir,'gi_taxid_nucl.dmp.gz'),'r')):
                sl = map(int,l.split())
                conn.execute(self.insert3,(sl[1],sl[0]))
                if verbose and i % 10000 == 0 and i > 0:
                    print >>sys.stderr, "Inserted record %d"%i

            if verbose:
                print >>sys.stderr, "Progress: Loading protein gimap"        
            
            for i,l in enumerate(gzip.open(os.path.join(datadir,'gi_taxid_prot.dmp.gz'),'r')):
                sl = map(int,l.split())
                conn.execute(self.insert4,(sl[1],sl[0]))
                if verbose and i % 10000 == 0 and i > 0:
                    print >>sys.stderr, "Inserted record %d"%i
            
        if verbose:
            print >>sys.stderr, "Progress: Making indicies"

        for l in self.indexes.splitlines():
            if 'gimap' not in l or totaxid:
                conn.execute(l)

        if verbose:
            print >>sys.stderr, "Progress: Commit"

        conn.commit()
        conn.close()

        if verbose:
            print >>sys.stderr, "Progress: Done"

    def init(self):
	assert self.exists(), "Taxonomy database missing"
        self.conn = sqlite3.connect(self.filename)
	self.conn.text_factory = str

    def normalize(self,n):
	if isinstance(n,int):
	    return str(int)
        return ' '.join(n.split()).lower()
    
    def get_taxid(self,name):
	if not name:
	    return None
        try:
            name = int(name)
        except:
            pass
        if type(name) == int:
            # print >>sys.stderr, self.select1, name
            for r in self.conn.execute(self.select1,(name,)):
                return r[0]
            # print >>sys.stderr, self.select2, name
            for r in self.conn.execute(self.select2,(name,)):
                return r[1]
        name = self.normalize(name)
        # print >>sys.stderr, self.select3, name
        for r in self.conn.execute(self.select3,(name,)):
            return r[0]
        return None

    def search(self,term):
	return map(self.get_taxid,map(itemgetter(0),self.conn.execute(self.select9,("%"+term+"%",))))

    def istaxid(self,tid):
        for r in self.conn.execute(self.select1,(tid,)):
            return True
        return False

    def get_scientific_name_(self,tid):
        # print >>sys.stderr, self.select4, tid
        for r in self.conn.execute(self.select4,(tid,)):
            return r[1]
        return None

    def get_scientific_name(self,name):
        tid = self.get_taxid(name)
        if tid is not None:
            return self.get_scientific_name_(tid)
        return None

    def get_synonyms(self,name):
        tid = self.get_taxid(name)
        if tid is not None:
            return map(itemgetter(1),self.conn.execute(self.select5,(tid,)))
        return None

    def get_rank(self,name):
        tid = self.get_taxid(name)
        if tid is not None:
            return self.get_node_(tid)[0]
        return None

    def get_parent(self,name):
        tid = self.get_taxid(name)
        if tid is not None:
            return self.get_node_(tid)[1]
        return None

    def get_node_(self,tid):
        for r in self.conn.execute(self.select1,(tid,)):
            return r[1],r[2]
        return None,None

    def get_taxonomy(self,name):
        tid = self.get_taxid(name)
        if tid is not None:
            tax = []
            while tid != 1:
                rank,parent = self.get_node_(tid)
                if not parent:
                    break
                if rank:
                    tax.append(self.get_scientific_name_(tid))
                tid = parent
            return ';'.join(tax)
        return None

    def get_ancestor_taxid(self,tid,max_rank=None,filter=True):
        tid = self.get_taxid(tid)
        retval = []
	rank = None
        while tid != 1:
            rank,parent = self.get_node_(tid)
            if not parent:
                break
            if rank or not filter:
                retval.append(tid)
                if rank == max_rank:
                    break
            tid = parent
        if rank == max_rank:
            return (retval,True)
        return (retval,False)

    def get_children_(self,tid):
        return map(itemgetter(0),self.conn.execute(self.select6,(tid,)))

    def get_children(self,tid):
        tid = self.get_taxid(tid)
        if tid != None:
            return self.get_children_(tid)
        return []

    def get_descendents(self,tid):
	retval = []
	for ch in self.get_children_(tid):
	    retval.append(ch)
	    retval.extend(self.get_descendents_(ch))
        return retval

    def get_descendents(self,tid):
        tid = self.get_taxid(tid)
        if tid:
            return self.get_descendants_(tid)
        return []

    def nucl2taxid(self,gi):
        for r in self.conn.execute(self.select7,(gi,)):
            return r[0]

    def prot2taxid(self,gi):
        for r in self.conn.execute(self.select8,(gi,)):
            return r[0]


class NCBITaxonomyTree(NCBITaxonomy):
    basefilename = 'ncbitaxonomytree'

    def build(self,force=False):
        return NCBITaxonomy._build(self,force=force,totaxid=False)

