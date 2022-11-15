
import sys, os, os.path, sqlite3, tempfile, itertools
from operator import itemgetter
from pickle import dumps, loads

class ProteinDatabase:
    prtable = """
    create table protein (
      id         integer primary key autoincrement,
      accession  varchar,
      dbindex    integer,
      data       blob
    );
    """
    insprot = """
    insert into protein (accession,dbindex) values (?,?);
    """
    peptable = """
    create table peptide (
      id integer primary key autoincrement,
      peptide varchar
    );
    """
    inspep = """
    insert into peptide (peptide) values (?);
    """
    getpep = """
    select id from peptide
    where peptide = ?;
    """
    getpep1 = """
    select id,peptide from peptide
    order by id limit ?,?
    """
    getpep2 = """
    select id,peptide from peptide
    order by id
    """
    pepcnt = """
    select count(*) from peptide
    """
    getpro0 = """
    select dbindex,accession,data from protein
    where id = ?
    """
    getpro = """
    select id from protein
    where dbindex = ? and accession = ?
    """
    getpro1 = """
    select id,accession from protein
    where dbindex = ?
    order by id limit ?,?
    """
    getpro2 = """
    select id,accession from protein
    where dbindex = ?
    order by id
    """
    getpro3 = """
    select id,data from protein
    where dbindex = ?
    and accession = ?
    order by id
    """
    getpro4 = """
    select id from protein
    """
    setpro = """
    update protein set data = ?
    where id = ?
    """
    hittable = """
    create table hit (
      protein_id integer,
      peptide_id integer,
      start integer,
      end integer,
      leftaa char,
      rightaa char,
      deltamass float,
      substitutions varchar
    )
    """
    inshit = """
    insert into hit (protein_id,peptide_id,start,end,leftaa,rightaa,deltamass,substitutions) values (?,?,?,?,?,?,?,?);
    """
    gethit = """
    select peptide.peptide, hit.start, hit.end, hit.leftaa, hit.rightaa, hit.deltamass, hit.substitutions
    from hit join peptide on (hit.peptide_id = peptide.id)
    where hit.protein_id = ? and
          peptide.id = ?
    """
    gethits = """
    select hit.protein_id,hit.start,hit.end,hit.leftaa,hit.rightaa, hit.deltamass, hit.substitutions
    from hit
    where
      hit.peptide_id = ?
    """
    getprothits = """
    select peptide.peptide, hit.start, hit.end, hit.leftaa, hit.rightaa, hit.deltamass, hit.substitutions
    from hit, peptide
    where
      hit.protein_id = ? and
      hit.peptide_id = peptide.id
    """
    getprothits1 = """
    select peptide.peptide, hit.start, hit.end, hit.leftaa, hit.rightaa, hit.deltamass, hit.substitutions
    from hit, peptide
    where
      hit.protein_id = ? and
      hit.peptide_id = peptide.id
    order by hit.start, hit.end
    """
    getprothits2 = """
    select peptide.peptide, hit.start, hit.end, hit.leftaa, hit.rightaa, hit.deltamass, hit.substitutions
    from hit, peptide
    where
      hit.protein_id = ? and
      hit.peptide_id = peptide.id
    order by hit.end, hit.start
    """
    indices = [_f for _f in map(str.strip,"""
    create %s index index1 on protein (accession,dbindex);
    create %s index index2 on peptide (peptide);
    create %s index index3 on hit (peptide_id,protein_id,start,end);
    create index index4 on hit (protein_id);
    """.split('\n')) if _f]
    def __init__(self,filename=None,flag='r'):
        self.delfilename = None
        if not filename:
            flag = 'n'
            fp,fn = tempfile.mkstemp(prefix='protdb',suffix='.db3')
            os.close(fp)
            self.delfilename = fn
            filename=fn
        assert(filename)
        if 'n' in flag:
            try:
                os.unlink(filename)
            except:
                pass
        self.db = sqlite3.connect(filename,
                                  isolation_level='EXCLUSIVE')
        self.db.text_factory = str
        cur = self.db.cursor()
        if 'n' in flag:
            cur.execute(self.prtable)
            cur.execute(self.peptable)
            cur.execute(self.hittable)
            self.db.commit()
            self.create_protein_indices(True)
            self.create_peptide_indices(True)
            self.create_hit_indices(True)
    def __del__(self):
        if self.delfilename:
            try:
                os.unlink(self.delfilename)
            except OSError:
                pass
    def drop_protein_indices(self):
        self.db.execute('drop index if exists index1;')
        self.db.commit()
    def create_protein_indices(self,unique=False):
        self.db.execute(self.indices[0]%("unique" if unique else ""))
        self.db.commit()
    def drop_peptide_indices(self):
        self.db.execute('drop index if exists index2;')
        self.db.commit()
    def create_peptide_indices(self,unique=False):
        self.db.execute(self.indices[1]%("unique" if unique else ""))
        self.db.commit()
    def drop_hit_indices(self):
        self.db.execute('drop index if exists index3;')
        self.db.execute('drop index if exists index4;')
        self.db.commit()
    def create_hit_indices(self,unique=False):
        self.db.execute(self.indices[2]%("unique" if unique else ""))
        self.db.commit()
        self.db.execute(self.indices[3])
        self.db.commit()
    def load_peptides(self,peptides):
        self.drop_peptide_indices()
        cur = self.db.cursor()
        cur.executemany(self.inspep,map(lambda p: (p,),peptides))
        self.db.commit()
        self.create_peptide_indices(True)
    def clean_peptides(self):
        self.db.execute("""
        delete from peptide
        where rowid in (
          select p2.rowid
          from peptide as p1, peptide as p2
          where
            p1.peptide = p2.peptide and
            p2.rowid > p1.rowid)
        """)
        self.db.commit()
    def load_proteins(self,accessions,dbind=0):
        self.drop_protein_indices()
        cur = self.db.cursor()
        cur.executemany(self.insprot,[(a,dbind) for a in accessions])
        self.db.commit()
        self.create_protein_indices(True)
    def clean_proteins(self):
        self.db.execute("""
        delete from protein
        where rowid in (
          select p2.id from protein as p1, protein as p2
          where p1.accession = p2.accession and
                p1.dbindex = p2.dbindex and
                p2.rowid > p1.rowid)
        """)
        self.db.commit()
    def load_hits(self,hits,dbind=0):
        # peps = dict((p,self.get_peptide_id(p)) for p in sorted(set(map(itemgetter(0),hits))))
        # pros = dict((a,self.get_protein_id(a,dbind)) for a in sorted(set(map(itemgetter(1),hits))))
        # hits = ((pros.get(h[1]),h[0],h[2],h[3],h[4],h[5]) for h in hits)
        self.drop_hit_indices()
        cur = self.db.cursor()
        cur.executemany(self.inshit,hits)
        self.db.commit()
        self.create_hit_indices(False)
    def clean_hits(self):
        self.db.execute("""
        delete from hit
        where rowid in (
          select h2.rowid from hit as h1, hit as h2
          where h1.protein_id = h2.protein_id and
                h1.peptide_id = h2.peptide_id and
                h1.start      = h2.start      and
                h1.end        = h2.end        and
                h2.rowid > h1.rowid)
        """)
        self.db.commit()
    def add_peptide(self,peptide):
        cur = self.db.cursor()
        try:
            cur.execute(self.inspep,(peptide,))
        except sqlite3.IntegrityError:
            pass
    def peptide_count(self):
        cur = self.db.cursor()
        for r in cur.execute(self.pepcnt):
            return r[0]
    def get_peptides(self,count=None,offset=0):
        cur = self.db.cursor()
        if count == None:
            q = cur.execute(self.getpep2)
        else:
            q = cur.execute(self.getpep1,(offset,count))
        return q
    def has_peptide(self,pepseq):
        cur = self.db.cursor()
        for r in cur.execute(self.getpep,(pepseq,)):
            return True
        return False
    def get_peptide_id(self,pepseq):
        cur = self.db.cursor()
        for r in cur.execute(self.getpep,(pepseq,)):
            return r[0]
        return None
    def add_protein(self,acc,dbind=0):
        cur = self.db.cursor()
        try:
            cur.execute(self.insprot,(acc,dbind))
        except sqlite3.IntegrityError:
            pass
    def get_proteins(self,count=None,offset=0,dbind=0):
        cur = self.db.cursor()
        if count == None:
            q = cur.execute(self.getpro2,(dbind,))
        else:
            q = cur.execute(self.getpro1,(dbind,offset,count))
        return q
    def protein_idmap(self):
        return dict(map(itemgetter(1,0),self.get_proteins()))
    def get_protein_id(self,accession,dbind=0):
        cur = self.db.cursor()
        for r in cur.execute(self.getpro,(dbind,accession)):
            return r[0]
        return None
    def has_protein(self,accession,dbind=0):
        cur = self.db.cursor()
        for r in cur.execute(self.getpro,(dbind,accession)):
            return True
        return False
    def get_protein_data(self,id):
        cur = self.db.cursor()
        for r in cur.execute(self.getpro0,(id,)):
            if r[2]:
                d = loads(r[2])
            else:
                d = {}
            return (r[0],r[1],d)
        return None
    def set_protein_data(self,id,**kw):
        cur = self.db.cursor()
        cur.execute(self.setpro,(dumps(kw),id))
    def update_protein_data(self,id,**kw):
        try:
            dbind,acc,data = self.get_protein_data(id)
        except TypeError:
            return
        data.update(kw)
        cur = self.db.cursor()
        cur.execute(self.setpro,(dumps(data),id))
    def store_hit(self,pep,prot,st,ed,laa,raa,dbind=0):
        pepid = self.get_peptide_id(pep)
        proid = self.get_protein_id(prot,dbind)
        cur = self.db.cursor()
        try:
            cur.execute(self.inshit,(proid,pepid,st,ed,laa,raa))
        except sqlite3.IntegrityError:
            pass
    def store_protein_hit(self,proid,pep,st,ed,laa,raa,dbind=0):
        pepid = self.get_peptide_id(pep)
        cur = self.db.cursor()
        try:
            cur.execute(self.inshit,(proid,pepid,st,ed,laa,raa))
        except sqlite3.IntegrityError:
            pass
    def get_peptide_hits(self,pepseq,dbind=None):
        pepid = self.get_peptide_id(pepseq)
        cur = self.db.cursor()
        if dbind != None:
            return [r for r in cur.execute(self.gethits,(pepid,)) if r[0] == dbind]
        return cur.execute(self.gethits,(pepid,))

    def has_hit(self,protid,pepseq):
        return self.get_hit(protid,pepseq) != None

    def get_hit(self,protid,pepseq):
        cur = self.db.cursor()
        pepid = self.get_peptide_id(pepseq)
        q = cur.execute(self.gethit,(protid,pepid))
        for r in q:
            return Hit(*r)
        return None

    def get_protein_hits(self,protid):
        cur = self.db.cursor()
        q = cur.execute(self.getprothits,(protid,))
        return [Hit(*r) for r in q]

    def get_protein_hits_by_start(self,protid):
        cur = self.db.cursor()
        q = cur.execute(self.getprothits1,(protid,))
        return [Hit(*r) for r in q]

    def get_protein_hits_by_end(self,protid):
        cur = self.db.cursor()
        q = cur.execute(self.getprothits2,(protid,))
        return [Hit(*r) for r in q]

    def commit(self):
        self.db.commit()

    def get(self,accession):
        id = self.get_protein_id(accession)
        if id:
            return Protein(self,id)
        return None

    def __iter__(self):
        return next(self)

    def __next__(self):
        cur = self.db.cursor()
        for r in cur.execute(self.getpro4):
            yield Protein(self,r[0])

import re

class Protein:
    def __init__(self,pdb,id):
        self.pdb = pdb
        self.id = id
        dbind,self.accession,d = pdb.get_protein_data(id)
        for k,v in d.items():
            setattr(self,k,v)

    @staticmethod
    def getAcc(defline,accfn=None,accre=None,accgrp=None):
        defline = defline.strip()
        if defline.startswith('>'):
            defline = defline[1:]
        if accfn != None:
            return accfn(defline)
        if accre != None:
            if accgrp == None:
                accgrp = 1
            m = re.search(accre,defline)
            assert(m != None)
            return m.group(accgrp)
        try:
            acc,desc = defline.split(None,1)
        except:
            acc = defline
        return acc

    @staticmethod
    def fromFasta(defline,sequence,pdb,**kw):
        if defline.startswith('>'):
            defline = defline[1:]
        acc = Protein.getAcc(defline,**kw)
        try:
            dummy,desc = defline.split(None,1)
        except:
            desc = ''
        if not pdb.has_protein(acc):
            pdb.add_protein(acc)
            id = pdb.get_protein_id(acc)
            pdb.set_protein_data(id,
                                 defline=defline,
                                 length=len(sequence),
                                 seq=sequence,
                                 description=desc)
        else:
            id = pdb.get_protein_id(acc)
            pdb.update_protein_data(id,
                                    defline=defline,
                                    length=len(sequence),
                                    seq=sequence,
                                    description=desc)
        return Protein(pdb,id)

    def findPep(self,pep):
        assert(self.seq != None)
        p = self.seq.find(pep)
        if p < 0:
            return
        try:
            laa = self.seq[p-1]
        except IndexError:
            laa = '-'
        try:
            raa = self.seq[p+len(pep)]
        except IndexError:
            raa = '-'
        self.pdb.add_peptide(pep)
        self.pdb.store_protein_hit(self.id,pep,p,p+len,laa,raa)

    # def addPepHit(self,pephit):
    #     if not self.pephits.has_key(pephit.peptide):
    #         self.pephits[pephit.peptide] = pephit

    def hasPepHit(self,pep):
        return self.pdb.has_hit(self.id,pep)

    def getPepHit(self,pep):
        return self.pdb.get_hit(self.id,pep)

    def pepHits(self):
        for h in self.pdb.get_protein_hits(self.id):
            yield h

    def pepHitsByStart(self):
        for h in self.pdb.get_protein_hits_by_start(self.id):
            yield h

    def pepHitsByEnd(self):
        for h in self.pdb.get_protein_hits_by_end(self.id):
            yield h

    def packedHits(self,filter=None):
        therows = [ ]
        for ph in self.pepHitsByStart():
            if filter != None and ph.peptide not in filter:
                continue
            inserted = False
            for r in therows:
                if r[-1].end < ph.start():
                    r.append(ph)
                    inserted = True
                    break
            if not inserted:
                therows.append([])
                therows[-1].append(ph)
        return therows

    def maxnonolap(self,filter=None):
        lastnonolappos = -1
        nonolap = 0;
        for ph in self.pepHitsByEnd():
            if filter != None and ph.peptide not in filter:
                continue
            if ph.start >= lastnonolappos:
                lastnonolappos = ph.end
                nonolap += 1
        return nonolap

    def coverage(self,filter=None):
        assert self.length != None
        aacov = 0
        lastcovpos = -1
        for ph in self.pepHitsByStart():
            if filter != None and ph.peptide not in filter:
                continue
            if ph.start >= lastcovpos:
                aacov += len(ph.peptide)
                lastcovpos = ph.end
            elif ph.end > lastcovpos:
                aacov += ph.end-lastcovpos
                lastcovpos = ph.end
        return float(aacov)/self.length

    def hitset(self,filter=None):
        return set([ph.peptide for ph in self.pepHits() if filter != None and ph.peptide in filter])

    def count(self,filter=None):
        return len(self.hitset(filter))

class Hit:
    def __init__(self,*args):
        self.peptide,self.start,self.end,self.leftaa,self.rightaa,self.deltamass,self.substitutions = args

if __name__ == "__main__":

    import sys, os

    protdb = ProteinDatabase(flag='n')
    protdb.load_peptides(map(str.strip,sys.stdin.readlines()))
    print(protdb.peptide_count())
    del protdb
