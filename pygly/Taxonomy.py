
import tempfile, os, os.path, sys, gzip, sqlite3, urllib, string, zipfile
from operator import itemgetter
from urllib.request import urlopen

class NCBITaxonomy(object):
    ncbitaxdmpzip = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip'
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
        create index index8 on merged (old_id);
    """
    insert1 = """
        insert into nodes values (?,?,?);
    """
    insert2 = """
        insert into names values (?,?,?,?);
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
    select9 = """
        select distinct tax_id from names
        where name LIKE ?        
    """
    def __init__(self,**kw):
        if kw.get('memory'):
             self.dbfile = None
        elif kw.get('filename'):
             self.dbfile = kw.get('filename')
        else:
             self.dbfile = '.ncbitaxa.db3'
        self.init(**kw)

    def init(self,verbose=True):
        if self.dbfile and os.path.exists(self.dbfile):
            if verbose:
                print("Progress: Connecting to cached NCBI taxonomy",file=sys.stderr)
            self.connect()
        else:
            tmpfile = self.download(verbose=verbose)
            self.build(tmpfile,verbose=verbose)
            os.unlink(tmpfile)

    def connect(self):
        if not self.dbfile:
            self.conn = sqlite3.connect(':memory:', isolation_level='EXCLUSIVE')
        else:
            self.conn = sqlite3.connect(self.dbfile, isolation_level='EXCLUSIVE')
        self.conn.text_factory = str

    def build(self,datafile,verbose=False):
        self.connect()
        try:
            self.parse(datafile,self.conn,verbose=verbose)
        except:
            raise

    def download(self,verbose=True):
        if verbose:
            print("Progress: Download taxonomy file",file=sys.stderr)

        filename = tempfile.NamedTemporaryFile(suffix='.zip',prefix='taxdmp_',delete=False).name
        h = urlopen(self.ncbitaxdmpzip)
        wh = open(filename,'wb')
        wh.write(h.read())
        wh.close()
        return filename

    def parse(self,filename,conn,verbose=False):
        if verbose:
            print("Progress: Creating tables",file=sys.stderr)

        conn.execute(self.createTable1)
        conn.execute(self.createTable2)
        conn.execute(self.createTable5)

        if verbose:
            print("Progress: Loading nodes",file=sys.stderr)

        zf = zipfile.ZipFile(filename,'r')
        for l in zf.read('nodes.dmp').splitlines():
            l = l.decode()
            l = l.rstrip('\t|\n')
            sl = l.split('\t|\t',3)
            taxid = int(sl[0])
            partaxid = int(sl[1])
            rank = sl[2].strip()
            if rank == 'no rank':
                rank = ""
            conn.execute(self.insert1,(taxid,rank,partaxid))

        if verbose:
            print("Progress: Loading names",file=sys.stderr)

        for i,l in enumerate(zf.read('names.dmp').splitlines()):
            l = l.decode() 
            l = l.rstrip('\t|\n')
            sl = list(map(str.strip,l.split('\t|\t')))
            taxid = int(sl[0])
            conn.execute(self.insert2,(taxid,sl[1],sl[3],1));
            if self.normalize(sl[1]) != sl[1]:
                conn.execute(self.insert2,(taxid,self.normalize(sl[1]),sl[3] + ' normalized',2));
            if sl[2]:
                conn.execute(self.insert2,(taxid,sl[2],sl[3] + ' unique',3));
                if self.normalize(sl[2]) != sl[2]:
                    conn.execute(self.insert2,(taxid,self.normalize(sl[2]),sl[3] + ' normalized unique',4));
            if verbose and i % 10000 == 0 and i > 0:
                print("Inserted record %d"%i,file=sys.stderr)
            
        if verbose:
            print("Progress: Loading merged taxids",file=sys.stderr)

        for l in zf.read('merged.dmp').splitlines():
            l = l.decode()
            l = l.rstrip('\t|\n')
            sl = list(map(int,l.split('\t|\t')))
            conn.execute(self.insert5,(sl[0],sl[1]));

        zf.close()

        if verbose:
            print("Progress: Making indicies",file=sys.stderr)

        for l in self.indexes.splitlines():
            conn.execute(l)

        if verbose:
            print("Progress: Commit",file=sys.stderr)

        conn.commit()

        if verbose:
            print("Progress: Done",file=sys.stderr)

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
            # print(self.select1, name,file=sys.stderr)
            for r in self.conn.execute(self.select1,(name,)):
                return r[0]
            # print(self.select2, name,file=sys.stderr)
            for r in self.conn.execute(self.select2,(name,)):
                return r[1]
        name = self.normalize(name)
        # print(self.select3, name,file=sys.stderr)
        for r in self.conn.execute(self.select3,(name,)):
            return r[0]
        return None

    def search(self,term):
        return list(map(self.get_taxid,map(itemgetter(0),self.conn.execute(self.select9,("%"+term+"%",)))))

    def istaxid(self,tid):
        for r in self.conn.execute(self.select1,(tid,)):
            return True
        return False

    def get_scientific_name_(self,tid):
        # print(self.select4, tid,file=sys.stderr)
        for r in self.conn.execute(self.select4,(tid,)):
            return r[1]
        return None

    def get_scientific_name(self,name):
        tid = self.get_taxid(name)
        if tid is not None:
            return self.get_scientific_name_(tid)
        return None

    def get_common_name(self,name,default=None):
        return self.get_synonym(name,'common name',default)

    def get_synonym(self,name,nameclass,default=None):
        names = self.get_synonyms(name)
        if names is not None:
            return names.get(nameclass,default)
        return None

    def get_synonyms(self,name):
        tid = self.get_taxid(name)
        if tid is not None:
            retval = {}
            for r in self.conn.execute(self.select5,(tid,)):
                retval[r[2]] = r[1]
            return retval
        return None

    def get_rank(self,name,default=None):
        tid = self.get_taxid(name)
        if tid is not None:
            r = self.get_node_(tid)[0]
            return (r if r else default)
        return default

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

    def get_ranked_ancestors(self,tid):
        tid = self.get_taxid(tid)
        retval = {}
        while tid != 1 and tid:
            rank,parent = self.get_node_(tid)
            if not parent:
                break
            if rank:
                retval[rank] = tid
            tid = parent
        return retval

    def get_genus(self,tid):
        return self.get_ranked_ancestors(tid).get('genus')

    def get_species(self,tid):
        return self.get_ranked_ancestors(tid).get('species')

    def get_family(self,tid):
        return self.get_ranked_ancestors(tid).get('family')

    def get_children_(self,tid):
        return list(map(itemgetter(0),self.conn.execute(self.select6,(tid,))))

    def get_children(self,tid):
        tid = self.get_taxid(tid)
        if tid != None:
            return self.get_children_(tid)
        return []

    def get_descendents_(self,tid):
        retval = []
        for ch in self.get_children_(tid):
            retval.append(ch)
            retval.extend(self.get_descendents_(ch))
        return retval

    def get_descendents(self,tid):
        tid = self.get_taxid(tid)
        if tid:
            return self.get_descendents_(tid)
        return []

if __name__ == "__main__":
    t = NCBITaxonomy(verbose=True)
    print(t.get_descendents(9606))
    for tid in t.search('sapiens'):
        print(tid,t.get_scientific_name(tid))
    for k,v in t.get_ranked_ancestors(9606).items():
        print(k,v)
