
from sqlobject import *
from sqlobject.dberrors import *
from sqlobject.sqlbuilder import LEFTJOINOn, LEFTJOIN, Select, Insert
import sys, os, tempfile, os.path, shutil, re, time
from collections import defaultdict
import pickle, base64, hashlib
from itertools import islice
from operator import attrgetter,itemgetter

defaults = dict(memory=False,
                clear=False,
                scheme='sqlite',
                host=None,
                user=None,
                password=None,
                port=None,
                dbname=None,
                filename=None,
                fromfile=False,
                clone=None,
                processConnection=False,
                threadConnection=True,
                debug=False)

def getopts(opts,kw):
    d = dict()
    for k in list(defaults.keys()):
        d[k] = getattr(opts,k)
    d.update(kw)
    return d

def init(opts=None,**kw):
    if not opts:
        theopts = {}
        theopts.update(defaults)
    else:
        theopts = getopts(opts,kw)
    theopts.update(kw)
    filename = None
    if theopts.get('memory',True):
        connuri = 'sqlite:/:memory:'
    elif theopts.get('scheme') == 'sqlite':
        if theopts.get('filename'):
            if theopts.get('clone'):
                assert(theopts.get('filename') != theopts.get('clone'))
                if os.path.exists(theopts.get('filename')):
                    os.unlink(theopts.get('filename'))
                shutil.copy2(theopts.get('clone'),theopts.get('filename'))
            elif theopts.get('clear'):
                if os.path.exists(theopts.get('filename')):
                    os.unlink(theopts.get('filename'))
            filename = os.path.abspath(theopts.get('filename'))
            connuri = 'sqlite://%s'%(filename,)
        else:
            fp,fn = tempfile.mkstemp(prefix='psmdb',suffix='.db3')
            os.close(fp)
            filename = os.path.abspath(fn)
            if theopts.get('clone'):
                shutil.copy(theopts.get('clone'),filename)
            connuri = 'sqlite://%s'%(filename,)
    else:
        connuri = '%(scheme)s://%(user)s:%(password)s@%(host)s:%(port)d/%(dbname)s'%theopts
    extras = []
    if not theopts.get('autoCommit',True):
        extras.append(('autoCommit',""))
    if theopts.get('debug',False):
        extras.append(('debug',1))
    if not theopts.get('cache',True):
        extras.append(('cache',""))
    if len(extras) > 0:
        connuri += '?' + '&'.join(["%s=%s"%(kv[0],kv[1]) for kv in extras])
    dbconnection = connectionForURI(connuri)
    if theopts.get('fromfile') and theopts.get('memory') and theopts.get('filename') and os.path.exists(theopts.get('filename')):
        if theopts.get('clone') and os.path.exists(theopts.get('clone')):
            filename = os.path.abspath(theopts.get('clone'))
        else:
            filename = os.path.abspath(theopts.get('filename'))
        conn0 = dbconnection.module.connect(filename)
        import sqlitebck
        print("Reading sqlite database from disk...", file=sys.stderr)
        sqlitebck.copy(conn0,dbconnection._memoryConn)
        print("...done.", file=sys.stderr)
    if theopts.get('processConnection',False):
        sqlhub.processConnection = dbconnection
    if theopts.get('threadConnection',False):
        sqlhub.threadConnection = dbconnection

    if theopts.get('clear'):
        for cls in tables:
            cls.dropTable(ifExists=True,connection=dbconnection)
    for cls in tables:
        if not theopts.get('cache',True):
            cls.sqlmeta.cacheValues = False
        cls.createTable(ifNotExists=True,connection=dbconnection)
    for cls in initTables:
        cls.initid(dbconnection)
    return dbconnection,filename

def flush(connection,filename):
    filename = os.path.abspath(filename)
    conn1 = connection.module.connect(filename)
    import sqlitebck
    print("Writing sqlite database to disk...", file=sys.stderr)
    sqlitebck.copy(connection._memoryConn,conn1)
    print("...done.", file=sys.stderr)

def doInitTables(connection):
    for cls in initTables:
        cls.initid(connection)

def memoize(key=None,precache=None,maxprecache=10000000,idstr="",reportinterval=10000):
    thecache = {}
    thecache['_calls'] = 0
    thecache['_misses'] = 0
    def wrap(f):
        def decorated_function(*args,**kw):
            if len(thecache) == 2:
                if precache:
                    # if idstr:
                    #     print >>sys.stderr, "%sPrecaching..."%(idstr + ": ",)
                    thecache.update(precache(kw.get('connection'),maxprecache))
                    if idstr:
                        print("%sCached %d entries"%(idstr + ": ",len(thecache)-2), file=sys.stderr)
            if key != None:
                _key = key(args,kw)
            else:
                _key = args
            if idstr:
                thecache['_calls'] += 1
            if _key not in thecache:
                if idstr:
                    thecache['_misses'] += 1
                # if len(thecache) > maxsize:
                #     thecache.clear()
                thecache[_key] = f(*args,**kw)
            if idstr and thecache['_calls'] % reportinterval == 0:
                print("%sMiss ratio: %.2f%% (calls: %d, misses: %d, cache size: %d)"%(idstr + ": ",100.0*float(thecache['_misses'])/thecache['_calls'],thecache['_calls'],thecache['_misses'],len(thecache)-2), file=sys.stderr)
            return thecache[_key]
        return decorated_function
    return wrap

def addarg(func):
    def wrap(f):
        def decorated_function(*args,**kw):
            args1 = list(args) + [func(args,kw)]
            return f(*args1,**kw)
        return decorated_function
    return wrap

def decoymap(d):
    if d == None:
        return 'NULL'
    return (bool(d)*1)

class MetadataContainer(SQLObject):
    metadata = PickleCol(length=8*1024,default=None)
    def hasdata(self,key):
        if self.metadata == None:
            return False
        return key in self.metadata
    def getdata(self,key,default=None):
        if self.metadata == None:
            return default
        return self.metadata.get(key,default)
    def setdata(self,key,value):
        d = self.metadata
        if d == None:
            d = {}
        d[key] = value
        self.metadata = d
    def unsetdata(self,key):
        d = self.metadata
        del d[key]
        self.metadata = d
    def updatedata(self,kw):
        d = self.metadata
        if d == None:
            d = {}
        d.update(kw)
        self.metadata = d
    def appenddata(self,key,value):
        d = self.metadata
        if d == None:
            d = {}
        if key not in d:
            d[key] = []
        d[key].append(value)
        self.metadata = d

class Analysis(MetadataContainer):
    name = StringCol(alternateID=True)
    psms = SQLMultipleJoin('PeptideSpectrumMatch')
    index1 = DatabaseIndex('name')

    @staticmethod
    @memoize()
    def insert(name,connection=None):
        try:
            return Analysis(name=name,connection=connection)
        except DuplicateEntryError:
            return Analysis.byName(name,connection=connection)

    @staticmethod
    def find(name,connection=None):
        try:
            f = Analysis.byName(name,connection=connection)
        except SQLObjectNotFound:
            try:
                f = Analysis(name=name,connection=connection)
            except DuplicateEntryError:
                f = Analysis.byName(name,connection=connection)
        return f

class SpectrumFileGroup(MetadataContainer):
    name = StringCol(alternateID=True)
    type = StringCol(default=None)
    files = RelatedJoin('SpectrumFile')
    index1 = DatabaseIndex('name')

    @staticmethod
    @memoize()
    def insert(name,type=None,connection=None):
        try:
            return SpectrumFileGroup(name=name,type=type,connection=connection)
        except DuplicateEntryError:
            return SpectrumFileGroup.byName(name,connection=connection)

    @staticmethod
    def find(name,type=None,connection=None):
        try:
            f = SpectrumFileGroup.byName(name,connection=connection)
        except SQLObjectNotFound:
            try:
                f = SpectrumFileGroup(name=name,type=type,connection=connection)
            except DuplicateEntryError:
                f = SpectrumFileGroup.byName(name,connection=connection)
        return f

class Sample(MetadataContainer):
    name = StringCol(alternateID=True)
    type = StringCol(default=None)
    files = MultipleJoin('SampleSpectrumFileLabel')
    index1 = DatabaseIndex('name')

    @staticmethod
    def insert(name,type=None,connection=None):
        try:
            return Sample(name=name,type=type,connection=connection)
        except DuplicateEntryError:
            return Sample.byName(name,connection=connection)

    @staticmethod
    def find(name,type=None,connection=None):
        try:
            f = Sample.byName(name,connection=connection)
        except SQLObjectNotFound:
            try:
                f = Sample(name=name,type=type,connection=connection)
            except DuplicateEntryError:
                f = Sample.byName(name,connection=connection)
        return f

class Label(MetadataContainer):
    name = StringCol(alternateID=True)
    type = StringCol(default=None)
    index1 = DatabaseIndex('name')

    def sortkey(self):
        m = re.search(r'^(\w+)-(\d+)([A-Z]?)$',self.name)
        assert m
        index = ["","N","C"].index(m.group(3))
        return (int(m.group(2)),index)

    @staticmethod
    def insert(name,type=None,connection=None):
        try:
            return Label(name=name,type=type,connection=connection)
        except DuplicateEntryError:
            return Label.byName(name,connection=connection)

    @staticmethod
    def find(name,type=None,connection=None):
        try:
            f = Label.byName(name,connection=connection)
        except SQLObjectNotFound:
            try:
                f = Label(name=name,type=type,connection=connection)
            except DuplicateEntryError:
                f = Label.byName(name,connection=connection)
        return f

class SampleSpectrumFileLabel(MetadataContainer):
    sample = ForeignKey('Sample',cascade=True)
    label = ForeignKey('Label',cascade=True)
    spectrumfile = ForeignKey('SpectrumFile',cascade=True)
    index1 = DatabaseIndex('sample','label','spectrumfile',unique=True)

    @staticmethod
    @memoize()
    def insert(sample,spectrumfile,label=None,connection=None):
        try:
            return SampleSpectrumFileLabel(sample=sample,label=label,spectrumfile=spectrumfile,
                                           connection=connection)
        except DuplicateEntryError:
            return SampleSpectrumFileLabel.selectBy(sample=sample,label=label,spectrumfile=spectrumfile,
                                                    connection=connection)

class SpectrumFile(MetadataContainer):
    groups = RelatedJoin('SpectrumFileGroup')
    samples = MultipleJoin('SampleSpectrumFileLabel',joinColumn='spectrumfile_id')
    name = StringCol(alternateID=True)
    spectra = MultipleJoin('Spectrum',joinColumn='file_id')

    def sample(self,label):
        try:
            ssfl = SampleSpectrumFileLabel.select(AND(SampleSpectrumFileLabel.q.spectrumfileID==self.id,
                                                      SampleSpectrumFileLabel.j.label,
                                                      Label.q.name==label))[0]
            return ssfl.sample
        except IndexError:
            return None

    @staticmethod
    @memoize()
    def insert(name,connection=None):
        try:
            return SpectrumFile(name=name,connection=connection)
        except DuplicateEntryError:
            return SpectrumFile.byName(name,connection=connection)

    @staticmethod
    def find(name,connection=None):
        try:
            f = SpectrumFile.byName(name,connection=connection)
        except SQLObjectNotFound:
            try:
                f = SpectrumFile(name=name,connection=connection)
            except DuplicateEntryError:
                f = SpectrumFile.byName(name,connection=connection)
        return f

    @staticmethod
    def access(name,connection=None):
        try:
            return SpectrumFile.byName(name,connection=connection)
        except SQLObjectNotFound:
            pass
        return None

    def label_type(self):
        ltype = None
        for s in self.samples:
            if ltype == None:
                ltype = s.label.type
            assert(ltype == s.label.type)
        return ltype

    def sample_names(self):
        return [s.sample.name for s in sorted(self.samples,key=lambda s: s.label.name)]

    def analytical_sample(self,delim=":"):
        return delim.join(self.sample_names())

class Spectrum(MetadataContainer):
    nextid = 1
    file = ForeignKey('SpectrumFile',cascade=True)
    start_scan = IntCol()
    end_scan   = IntCol()
    precursorMz = FloatCol(default=None)
    psms = MultipleJoin('PeptideSpectrumMatch')
    index1 = DatabaseIndex('file','start_scan','end_scan',unique=True)
    index2 = DatabaseIndex('precursorMz')

    def distinct_peptides(self):
        return set(psm.peptideIon.peptide for psm in self.psms)

    def nonequiv_peptides(self):
        nonequiv = {}
        for pep in set(psm.peptideIon.peptide for psm in self.psms):
            cannonseq = pep.cannon_ILKQ()
            if cannonseq not in nonequiv:
                nonequiv[cannonseq] = pep
        return set(nonequiv.values())

    def distinct_peptideIons(self):
        return set(psm.peptideIon for psm in self.psms)

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from spectrum;
        """)[0]
        if not maxid:
            maxid = 0
        Spectrum.nextid = maxid+1

    @staticmethod
    @memoize()
    def insert(fileid,start_scan,end_scan,precursorMz=None,metadata=None,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = Spectrum.nextid
            connection.query("""
                insert into spectrum (id,file_id,start_scan,end_scan,precursor_mz,metadata) values (%d,%d,%d,%d,%s,%s);
            """%(id,fileid,start_scan,end_scan,'NULL' if precursorMz == None else "%.10f"%precursorMz,
                "\""+base64.encodebytes(pickle.dumps(metadata)).decode()+"\"" if metadata else 'NULL'))
            Spectrum.nextid += 1
            return id
        except DuplicateEntryError:
            return Spectrum.select(AND(Spectrum.q.fileID==fileid,
                                       Spectrum.q.start_scan==start_scan,
                                       Spectrum.q.end_scan==end_scan),
                                   connection=connection)[0].id

    @staticmethod
    def find(file,start_scan,end_scan,connection=None):
        try:
            s = Spectrum.select(AND(Spectrum.q.file==file,
                                    Spectrum.q.start_scan==start_scan,
                                    Spectrum.q.end_scan==end_scan),
                                connection=connection)[0]
        except IndexError:
            try:
                s = Spectrum(file=file,start_scan=start_scan,end_scan=end_scan,
                             connection=connection)
            except DuplicateEntryError:
                s = Spectrum.select(AND(Spectrum.q.file==file,
                                        Spectrum.q.start_scan==start_scan,
                                        Spectrum.q.end_scan==end_scan),
                                    connection=connection)[0]
        return s

    @staticmethod
    def access(file,start_scan,end_scan,connection=None):
        try:
            Spectrum.select(AND(Spectrum.q.file==file,
                                Spectrum.q.start_scan==start_scan,
                                Spectrum.q.end_scan==end_scan),
                            connection=connection)[0]
        except IndexError:
            pass
        return None

    @staticmethod
    def normalize(quiet=False,connection=None):
        orphanids = set(map(attrgetter('id'),
                             Spectrum.select(PeptideSpectrumMatch.q.id == None,
                                             join=LEFTJOINOn(Spectrum,PeptideSpectrumMatch,
                                                             PeptideSpectrumMatch.j.spectrum),
                                             connection=connection)))
        # For file ForeignKey constraint
        orphanids1 = set(map(attrgetter('id'),
                              Spectrum.select(SpectrumFile.q.id == None,
                                              join=LEFTJOINOn(Spectrum,SpectrumFile,
                                                              Spectrum.j.file),
                                              connection=connection)))
        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        Spectrum.deleteMany(IN(Spectrum.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d spectra"%(count), file=sys.stderr)
        return count

class PeptideSpectrumMatch(MetadataContainer):
    nextid = 1
    spectrum = ForeignKey('Spectrum',cascade=True)
    peptideIon = ForeignKey('PeptideIon',cascade=True)
    analysis = ForeignKey('Analysis',cascade=True)
    index1 = DatabaseIndex('spectrum','peptideIon','analysis',unique=True)
    index2 = DatabaseIndex('value')
    index3 = DatabaseIndex('analysis')
    index4 = DatabaseIndex('peptideIon')
    value =  FloatCol(default=None)

    def setdelta(self):
        pi = self.peptideIon
        tmw = (pi.mw+1.0078*pi.charge)
        self.setdata('delta',(self.spectrum.precursorMz*pi.charge - tmw))
        self.setdata('deltappm',1e6*(self.delta)/tmw)

    def c13_delta_mass(self,minc13=0,maxc13=2):
        if not self.hasdata('delta'):
            self.setdelta()
        d = self.getdata('delta',None)
        i = int(round(d))
        if i < minc13:
            i = maxc13+1
            d1 = (d - minc13*1.0025)
        elif i > maxc13:
            i = maxc13+1
            d1 = (d - maxc13*1.0025)
        else:
            d1 = (d - i*1.0025)
        return i,d1

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from peptide_spectrum_match;
        """)[0]
        if not maxid:
            maxid = 0
        PeptideSpectrumMatch.nextid = maxid+1

    @staticmethod
    @memoize()
    def insert(spectrumid,peptideionid,analysisid=None,value=None,metadata=None,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = PeptideSpectrumMatch.nextid
            connection.query("""
                insert into peptide_spectrum_match (id,spectrum_id,peptide_ion_id,analysis_id,value,metadata)
                values (%d,%d,%d,%s,%s,%s)
            """%(id,spectrumid,peptideionid,analysisid if analysisid else 'NULL',
                 'NULL' if value == None else "%.10f"%value,
                 "\""+base64.encodebytes(pickle.dumps(metadata)).decode()+"\"" if metadata else 'NULL'))
            PeptideSpectrumMatch.nextid += 1
            return id
            # return PeptideSpectrumMatch(spectrumID=spectrumid,
            #                             peptideIon=peptideion,
            #                             analysis=analysis,
            #                             value=value,
            #                             metadata=metadata,
            #                             connection=connection)
        except DuplicateEntryError:
            return PeptideSpectrumMatch.select(AND(PeptideSpectrumMatch.q.spectrumID==spectrumid,
                                                   PeptideSpectrumMatch.q.peptideIonID==peptideionid,
                                                   PeptideSpectrumMatch.q.analysisID==analysisid),
                                               connection=connection)[0].id

    @staticmethod
    def find(spectrum,peptideion,analysis=None,connection=None):
        try:
            psm = PeptideSpectrumMatch.select(AND(PeptideSpectrumMatch.q.spectrum==spectrum,
                                                  PeptideSpectrumMatch.q.peptideIon==peptideion,
                                                  PeptideSpectrumMatch.q.analysis==analysis),
                                              connection=connection)[0]
        except IndexError:
            try:
                psm = PeptideSpectrumMatch(spectrum=spectrum,
                                           peptideIon=peptideion,
                                           analysis=analysis,
                                           connection=connection)
            except DuplicateEntryError:
                psm = PeptideSpectrumMatch.select(AND(PeptideSpectrumMatch.q.spectrum==spectrum,
                                                      PeptideSpectrumMatch.q.peptideIon==peptideion,
                                                      PeptideSpectrumMatch.q.analysis==analysis),
                                                  connection=connection)[0]
        return psm

    @staticmethod
    def filter(minvalue=None,maxvalue=None,connection=None):
        if maxvalue != None:
            PeptideSpectrumMatch.deleteMany(PeptideSpectrumMatch.q.value > maxvalue,connection=connection)
        if minvalue != None:
            PeptideSpectrumMatch.deleteMany(PeptideSpectrumMatch.q.value < minvalue,connection=connection)

    @staticmethod
    def normalize(quiet=False,connection=None):
        # spectrum foreign key constraint
        orphanids = set(map(attrgetter('id'),
                             PeptideSpectrumMatch.select(Spectrum.q.id == None,
                                                         join=LEFTJOINOn(PeptideSpectrumMatch,Spectrum,
                                                                         PeptideSpectrumMatch.j.spectrum),
                                                         connection=connection)))

        # peptideIon foreign key constraint
        orphanids1 = set(map(attrgetter('id'),
                              PeptideSpectrumMatch.select(PeptideIon.q.id == None,
                                                          join=LEFTJOINOn(PeptideSpectrumMatch,PeptideIon,
                                                                          PeptideSpectrumMatch.j.peptideIon),
                                                          connection=connection)))

        # analysis foreign key constraint
        orphanids2 = set(map(attrgetter('id'),
                              PeptideSpectrumMatch.select(Analysis.q.id == None,
                                                          join=LEFTJOINOn(PeptideSpectrumMatch,Analysis,
                                                                          PeptideSpectrumMatch.j.analysis),
                                                          connection=connection)))

        orphanids |= orphanids1
        orphanids |= orphanids2

        if len(orphanids) == 0:
            return 0

        PeptideSpectrumMatch.deleteMany(IN(PeptideSpectrumMatch.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d peptide spectrum matches"%(count), file=sys.stderr)
        return count


class PeptideIon(MetadataContainer):
    nextid = 1
    mz = FloatCol(default=None)
    charge = IntCol()
    mw = FloatCol(default=None)
    peptide = ForeignKey('Peptide',cascade=True)
    modifiedSites = MultipleJoin('ModificationSite')
    modString = StringCol()
    psms = MultipleJoin('PeptideSpectrumMatch')
    index1 = DatabaseIndex('peptide','modString','charge',unique=True)

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from peptide_ion;
        """)[0]
        if not maxid:
            maxid = 0
        PeptideIon.nextid = maxid+1

    @staticmethod
    def _mods2str(mods):
        if len(mods) == 0:
            return ""
        return ','.join(["%s%d:%s"%(m[1].aminoAcid,m[0],m[1].massDeltaString) for m in sorted(mods,key=lambda m: (m[0],m[1].massDelta))])

    @staticmethod
    @addarg(lambda a,k: PeptideIon._mods2str(a[1]))
    @memoize(idstr="PeptideIon",
             key=lambda a,k: (a[0], a[3], a[2]),
             precache=lambda c,n: map(lambda r: ((r[0],r[1],r[2]),r[3]),
                                     c.queryAll(c.sqlrepr(Select([PeptideIon.q.peptideID,
                                                                  PeptideIon.q.modString,
                                                                  PeptideIon.q.charge,
                                                                  PeptideIon.q.id],orderBy='id',start=0,end=n)))))

    def insert(peptideid,mods,z,modString,mw=None,mz=None,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = PeptideIon.nextid
            connection.query("""
                insert into peptide_ion (id,peptide_id,mod_string,charge,mw,mz) values (%d,%d,"%s",%s,%s,%s)
            """%(id,peptideid,modString,'NULL' if z == None else "%d"%z,'NULL' if mw == None else "%.10f"%mw,'NULL' if mz == None else "%.10f"%mz))
            # ion = PeptideIon(peptideID=peptideid,modString=modString,charge=z,mw=mw,mz=mz,connection=connection)
            for m in mods:
                ModificationSite.insert(id,m[0],m[1].id,connection=connection)
            PeptideIon.nextid += 1
            return id
        except DuplicateEntryError:
            return PeptideIon.select(AND(PeptideIon.q.peptideID==peptideid,
                                         PeptideIon.q.modString==modString,
                                         PeptideIon.q.charge==z),
                                     connection=connection)[0].id

    @staticmethod
    def find(peptide,mods,z,connection=None):
        modString=PeptideIon._mods2str(mods)
        new = False
        try:
            ion = PeptideIon.select(AND(PeptideIon.q.peptide==peptide,
                                        PeptideIon.q.modString==modString,
                                        PeptideIon.q.charge==z),
                                    connection=connection)[0]
        except IndexError:
            try:
                ion = PeptideIon(peptide=peptide,modString=modString,charge=z,connection=connection)
                new = True
            except DuplicateEntryError:
                ion = PeptideIon.select(AND(PeptideIon.q.peptide==peptide,
                                            PeptideIon.q.modString==modString,
                                            PeptideIon.q.charge==z),
                                        connection=connection)[0]
        if new:
            for m in mods:
                try:
                    ms = ModificationSite.find(ion,m[0],m[1].id,connection=connection)
                except DuplicateEntryError:
                    pass
        return ion

    def analytical_samples(self):
        samples = set()
        for s in self.spectra():
            samples.add(s.file.analytical_sample())
        return samples

    def spectra(self):
        spec = set()
        for psm in self.psms:
            spec.add(psm.spectrum)
        return spec

    def modifiedPeptide(self):
        pl = list(self.peptide.sequence)
        for ms in self.modifiedSites:
            if 0 < ms.position <= len(pl):
                pl[ms.position-1] = pl[ms.position-1].lower()
        return ''.join(pl)

    def annotatedPeptide(self,modsym={},before=True):
        pl = list(self.peptide.sequence)
        for ms in self.modifiedSites:
            if 0 < ms.position <= len(pl):
                if ms.modification in modsym:
                    if before:
                        pl[ms.position-1] = modsym[ms.modification] + pl[ms.position-1]
                    else:
                        pl[ms.position-1] += modsym[ms.modification]
        return ''.join(pl)

    @staticmethod
    def normalize(quiet=False,connection=None):
        orphanids = set(map(attrgetter('id'),
                             PeptideIon.select(PeptideSpectrumMatch.q.id == None,
                                               join=LEFTJOINOn(PeptideIon,PeptideSpectrumMatch,
                                                               PeptideSpectrumMatch.j.peptideIon),
                                               connection=connection)))
        # peptide foreign key constraint
        orphanids1 = set(map(attrgetter('id'),
                              PeptideIon.select(Peptide.q.id == None,
                                                join=LEFTJOINOn(PeptideIon,Peptide,
                                                                PeptideIon.j.peptide),
                                                connection=connection)))

        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        PeptideIon.deleteMany(IN(PeptideIon.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d peptide ions"%(count), file=sys.stderr)
        return count

class ModificationSite(MetadataContainer):
    nextid = 1
    position = IntCol()
    peptideIon = ForeignKey('PeptideIon',cascade=True)
    modification = ForeignKey('Modification',cascade=True)
    index1 = DatabaseIndex('peptideIon','position','modification',unique=True)
    index2 = DatabaseIndex('modification')

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from modification_site;
        """)[0]
        if not maxid:
            maxid = 0
        ModificationSite.nextid = maxid+1

    @staticmethod
    @memoize(idstr="ModificationSite")
    def insert(peptideionid, position, modificationid, connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = ModificationSite.nextid
            connection.query("""
                insert into modification_site (id,position,peptide_ion_id,modification_id) values (%d,%d,%d,%d)
            """%(id,position,peptideionid,modificationid))
            ModificationSite.nextid += 1
            return id
            # return ModificationSite(peptideIon=peptideion,position=position,modification=modification,
            #                         connection=connection)
        except DuplicateEntryError:
            return ModificationSite.select(AND(ModificationSite.q.peptideIon == peptideion,
                                               ModificationSite.q.position == position,
                                               ModificationSite.q.modification == modification),
                                           connection=connection)[0].id

    @staticmethod
    def find(peptideion, position, modification,connection=None):
        try:
            ms = ModificationSite.select(AND(ModificationSite.q.peptideIon == peptideion,
                                             ModificationSite.q.position == position,
                                             ModificationSite.q.modification == modification),
                                         connection=connection)[0]
        except IndexError:
            try:
                ms = ModificationSite(peptideIon=peptideion,position=position,modification=modification,
                                      connection=connection)
            except DuplicateEntryError:
                ms = ModificationSite.select(AND(ModificationSite.q.peptideIon == peptideion,
                                                 ModificationSite.q.position == position,
                                                 ModificationSite.q.modification == modification),
                                             connection=connection)[0]
        return ms

    @staticmethod
    def normalize(quiet=False,connection=None):
        # peptideIon foreign key constraint
        orphanids = set(map(attrgetter('id'),
                             ModificationSite.select(PeptideIon.q.id == None,
                                                     join=LEFTJOINOn(ModificationSite,PeptideIon,
                                                                     ModificationSite.j.peptideIon),
                                                     connection=connection)))
        # modification foreign key constraint
        orphanids1 = set(map(attrgetter('id'),
                              ModificationSite.select(Modification.q.id == None,
                                                      join=LEFTJOINOn(ModificationSite,Modification,
                                                                      ModificationSite.j.modification),
                                                      connection=connection)))

        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        ModificationSite.deleteMany(IN(ModificationSite.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d modification sites"%(count), file=sys.stderr)
        return count



class Modification(MetadataContainer):
    aminoAcid = StringCol(length=1,varchar=False)
    massDelta = FloatCol()
    massDeltaString = StringCol()
    index1 = DatabaseIndex('aminoAcid','massDeltaString',unique=True)

    def peptideIons(self,connection=None):
        lastpi = None
        for ms in ModificationSite.select(ModificationSite.q.modification==self,
                                          orderBy=ModificationSite.q.peptideIon,
                                          connection=connection):
            if ms.peptideIonID != lastpi:
                yield ms.peptideIon
                lastpi = ms.peptideIonID

    @staticmethod
    def _deltostr(delta):
        return "%+.3f"%(delta,)

    @staticmethod
    @addarg(lambda a,k: Modification._deltostr(a[1]))
    @memoize(key=lambda a,k: (a[0],a[2]))
    def insert(aa,delta,deltaStr,connection=None):
        try:
            return Modification(aminoAcid=aa,massDeltaString=deltaStr,massDelta=delta,
                                connection=connection)
        except DuplicateEntryError:
            return Modification.select(AND(Modification.q.aminoAcid==aa,
                                           Modification.q.massDeltaString==deltaStr),
                                       connection=connection)[0]

    @staticmethod
    def access(aa,delta,connection=None):
        massDeltaString = Modification._deltostr(delta)
        try:
            m = Modification.select(AND(Modification.q.aminoAcid==aa,
                                        Modification.q.massDeltaString==massDeltaString),
                                    connection=connection)[0]
        except IndexError:
            m = None
        return m

    @staticmethod
    def find(aa,delta,connection=None):
        massDeltaString = Modification._deltostr(delta)
        try:
            m = Modification.select(AND(Modification.q.aminoAcid==aa,
                                        Modification.q.massDeltaString==massDeltaString),
                                    connection=connection)[0]
        except IndexError:
            try:
                m = Modification(aminoAcid=aa,massDeltaString=massDeltaString,massDelta=delta,
                                 connection=connection)
            except DuplicateEntryError:
                m = Modification.select(AND(Modification.q.aminoAcid==aa,
                                            Modification.q.massDeltaString==massDeltaString),
                                        connection=connection)[0]
        return m

class Peptide(MetadataContainer):
    nextid = 1
    sequence = StringCol(alternateID=True)
    decoy = BoolCol(default=None)
    ions = MultipleJoin('PeptideIon',joinColumn='peptide_id')
    alignments = MultipleJoin('Alignment',joinColumn='peptide_id')
    membership = MultipleJoin('Peptide2PeptideGroup',joinColumn='peptide_id')
    index1 = DatabaseIndex(dict(column=sequence,length=50))

    @staticmethod
    def isobaric(p1,p2,tolerance=1e-5):
        from load import mt
        return (abs(mt.peptide_mr(p1.sequence)-mt.peptide_mr(p2.sequence)) < tolerance)

    def cannon_IL(self):
        return self.sequence.replace('I','L')

    def cannon_KQ(self):
        s = self.sequence
        return s[:-1].replace('K','Q')+s[-1]

    def cannon_ILKQ(self):
        s = self.sequence.replace('I','L')
        return s[:-1].replace('K','Q')+s[-1]

    @staticmethod
    def equiv_IL(p1,p2):
        return (p1.cannon_IL() == p2.cannon_IL())

    @staticmethod
    def equiv_KQ(p1,p2):
        return (p1.cannon_KQ() == p2.cannon_KQ())

    @staticmethod
    def equiv_ILKQ(p1,p2):
        return (p1.cannon_ILKQ() == p2.cannon_ILKQ())

    def analytical_samples(self):
        samples = set()
        for f in self.spectrum_files():
            samples.add(f.analytical_sample())
        return samples

    def spectrum_files(self):
        files = set()
        for s in self.spectra():
            files.add(s.file)
        return files

    def spectra(self):
        spec = set()
        for ions in self.ions:
            for psm in ions.psms:
                spec.add(psm.spectrum)
        return spec

    def psms(self,spectrum=None):
        psms = set()
        for ions in self.ions:
            for psm in ions.psms:
                if not spectrum or psm.spectrum == spectrum:
                    psms.add(psm)
        return psms

    def proteins(self):
        return list(set([a.protein for a in self.alignments]))

    def genes(self):
        connection = self._connection
        q = connection.sqlrepr(Select(
            [ProteinGroup.q.id],
            AND(Alignment.q.peptideID==self.id,
                ProteinGroup.q.type=='Gene',
                Protein2ProteinGroup.j.group,
                Alignment.q.proteinID == Protein2ProteinGroup.q.proteinID)))
        l = set()
        for r in connection.queryAll(q):
            l.add(r[0])
        return [ProteinGroup.get(id,connection=connection) for id in l]

    def alignmentsto(self,protein,locus=None,loci=None):
        if locus != None:
            return [a for a in self.alignments if a.protein == protein and a.start <= locus < a.end]
        elif loci != None:
            return [a for a in self.alignments if a.protein == protein and a.start <= min(loci) and max(loci) < a.end]
        return [a for a in self.alignments if a.protein == protein]

    def groups(self):
        return list(m.group for m in self.membership)

    def modsites(self,mods,psmcondition=None):
        positions = set()
        for psm in self.psms():
            if psmcondition == None or psmcondition(psm):
                ion = psm.peptideIon
                for ms in ion.modifiedSites:
                    if ms.modification in mods:
                        positions.add(ms.position)
        return positions

    @staticmethod
    def gene2peptides0(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [ProteinGroup.q.id,
             func.GROUP_CONCAT(func.DISTINCT(Alignment.q.peptideID))],
             AND(Alignment.j.protein,
                 Protein2ProteinGroup.j.protein,
                 Protein2ProteinGroup.j.group,
                 ProteinGroup.q.type=='Gene'),
            groupBy=[ProteinGroup.q.id]))
        for r in connection.queryAll(q):
            yield r[0],set(map(int,r[1].split(',')))

    @staticmethod
    def gene2peptides(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [ProteinGroup.q.id,Alignment.q.peptideID],
             AND(ProteinGroup.q.type=='Gene',
                 Protein2ProteinGroup.j.group,
                 Alignment.q.proteinID == Protein2ProteinGroup.q.proteinID)))
        d = defaultdict(set)
        for r in connection.queryAll(q):
            d[r[0]].add(r[1])
        return d

    @staticmethod
    def peptide2genes(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Alignment.q.peptideID,ProteinGroup.q.id],
             AND(ProteinGroup.q.type=='Gene',
                 Protein2ProteinGroup.j.group,
                 Alignment.q.proteinID == Protein2ProteinGroup.q.proteinID)))
        d = defaultdict(set)
        for r in connection.queryAll(q):
            d[r[0]].add(r[1])
        return d

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from peptide;
        """)[0]
        if not maxid:
            maxid = 0
        Peptide.nextid = maxid+1

    @staticmethod
    @memoize(key=lambda a,k: a[0],
             precache=lambda c,n: c.queryAll(c.sqlrepr(Select([Peptide.q.sequence,Peptide.q.id],orderBy='id',start=0,end=n))),
             idstr="Peptide")
    def insert(peptide,decoy=None,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = Peptide.nextid
            connection.query("insert into peptide (id,sequence,decoy) values (%d,\"%s\",%s);"%(id,peptide,decoymap(decoy)))
            Peptide.nextid += 1
            return id
            # return Peptide(sequence=peptide,connection=connection)
        except DuplicateEntryError:
            return Peptide.select(Peptide.q.sequence==peptide,
                                  connection=connection)[0].id

    @staticmethod
    def find(peptide,decoy=None,connection=None):
        try:
            p = Peptide.select(Peptide.q.sequence==peptide,
                                       connection=connection)[0]
        except IndexError:
            try:
                p = Peptide(sequence=peptide,decoy=decoymap(decoy),connection=connection)
            except DuplicateEntryError:
                p = Peptide.select(Peptide.q.sequence==peptide,
                                   connection=connection)[0]
        return p

    @staticmethod
    def access(peptide,connection=None):
        try:
            return Peptide.select(Peptide.q.sequence==peptide,
                                          connection=connection)[0]
        except IndexError:
            pass
        return None

    @staticmethod
    def normalize(alignedOnly=True,quiet=False,connection=None):
        orphanids = set()
        if alignedOnly:
            orphanids = set(map(attrgetter('id'),
                                 Peptide.select(Alignment.q.id == None,
                                                join=LEFTJOINOn(Peptide,Alignment,Alignment.j.peptide),
                                                connection=connection)))

        orphanids1 = set(map(attrgetter('id'),
                              Peptide.select(PeptideIon.q.id == None,
                                             join=LEFTJOINOn(Peptide,PeptideIon,PeptideIon.j.peptide),
                                             connection=connection)))

        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        Peptide.deleteMany(IN(Peptide.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d peptides"%(count), file=sys.stderr)
        return count

class Peptide2PeptideGroup(SQLObject):
    nextid = 1
    peptide = ForeignKey('Peptide',cascade=True)
    group = ForeignKey('PeptideGroup',cascade=True)
    index1 = DatabaseIndex('peptide','group',unique=True)
    index1 = DatabaseIndex('group','peptide')

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from peptide2_peptide_group;
        """)[0]
        if not maxid:
            maxid = 0
        Peptide2PeptideGroup.nextid = maxid+1

    @staticmethod
    @memoize()
    def insert(peptideid,groupid,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = Peptide2PeptideGroup.nextid
            connection.query("""
                insert into peptide2_peptide_group (id,peptide_id,group_id) values (%d,%d,%d)
            """%(id,peptideid,groupid))
            Peptide2PeptideGroup.nextid += 1
            return id
            #return Peptide2PeptideGroup(peptideID=peptideid,groupID=groupid,connection=connection).id
        except DuplicateEntryError:
            return Peptide2PeptideGroup.selectBy(peptideID=peptideid,groupID=groupid,connection=connection)[0].id

    @staticmethod
    def clear(connection=None):
        Peptide2PeptideGroup.deleteMany(None,connection=connection)

    @staticmethod
    def normalize(quiet=False,connection=None):
        # peptide foreign key
        orphanids = set(map(attrgetter('id'),
                             Peptide2PeptideGroup.select(Peptide.q.id == None,
                                                         join=LEFTJOINOn(Peptide2PeptideGroup,Peptide,
                                                                         Peptide2PeptideGroup.j.peptide),
                                                         connection=connection)))

        # peptide group foreign key
        orphanids1 = set(map(attrgetter('id'),
                              Peptide2PeptideGroup.select(PeptideGroup.q.id == None,
                                                          join=LEFTJOINOn(Peptide2PeptideGroup,PeptideGroup,
                                                                          Peptide2PeptideGroup.j.group),
                                                          connection=connection)))

        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        Peptide2PeptideGroup.deleteMany(IN(Peptide2PeptideGroup.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d peptide group memberships"%(count), file=sys.stderr)
        return count

class PeptideGroup(MetadataContainer):
    name = StringCol(alternateID=True,length=1024)
    type = StringCol(default=None)
    count = IntCol(default=0)
    members = MultipleJoin('Peptide2PeptideGroup',joinColumn='group_id')
    index0 = DatabaseIndex('name')
    index1 = DatabaseIndex('count')
    index2 = DatabaseIndex('type')

    def peptides(self):
        return list(m.peptide for m in self.members)

    def spectra(self):
        if len(self.members) == 0:
            return []
        thespectra = set()
        pep0 = self.members[0].peptide
        pepids = set(m.peptide for m in self.members)
        for s in Spectrum.select(AND(PeptideSpectrumMatch.j.spectrum,
                                     PeptideSpectrumMatch.j.peptideIon,
                                     PeptideIon.j.peptideID==pep0.id),
                                 connection=self._connection):
            if s.distinct_peptides() == pepids:
                thespectra.add(s)
        return thespectra

    @staticmethod
    def distinct_proteins(peptides,connection):
        pepids = list(map(attrgetter('id'),peptides))
        return Protein.select(AND(IN(Peptide.q.id,pepids),
                                  Alignment.j.peptide,Alignment.j.protein),
                              distinct=True,
                              connection=connection)

    def proteins(self):
        return PeptideGroup.distinct_proteins((m.peptide for m in self.members),self._connection)

    @staticmethod
    def insertSpectraGroups(type,prefix,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [func.DISTINCT(func.GROUP_CONCAT(func.DISTINCT(PeptideIon.q.peptideID)))],
            AND(PeptideSpectrumMatch.j.spectrum,PeptideSpectrumMatch.j.peptideIon),
            groupBy=[Spectrum.q.id]))
        try:
            nextid = PeptideGroup.select(orderBy='-id',connection=connection)[0].id;
        except IndexError:
            nextid = 0
        for r in connection.queryAll(q):
            nextid += 1
            pepids = sorted(map(int,r[0].split(',')))
            name = ':'.join(map(str,pepids))
            name = ":".join([prefix,hashlib.sha1(name.encode('utf8')).hexdigest().lower(),
                                    hashlib.md5(name.encode('utf8')).hexdigest().lower(),
                                    str(len(pepids))])
            assert len(name) < 1024, "name is too long: %d chars, %d peptide ids"%(len(name), len(pepids))
            try:
                connection.query(connection.sqlrepr(Insert(PeptideGroup.sqlmeta.table,
                                                           values=dict(list(zip(['id','name','type','count'],
                                                                           (nextid,name,type,len(pepids))))))))
            except DuplicateEntryError:
                continue
            for pepid in pepids:
                connection.query(connection.sqlrepr(Insert(Peptide2PeptideGroup.sqlmeta.table,
                                                           values=dict(list(zip(['peptide_id','group_id'],
                                                                           (pepid,nextid)))))))

    @staticmethod
    def protein2groups(type=None,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Alignment.q.proteinID,
             func.GROUP_CONCAT(func.DISTINCT(Peptide2PeptideGroup.q.groupID))],
            AND(Peptide2PeptideGroup.q.peptideID == Alignment.q.peptideID),
            groupBy=[Alignment.q.proteinID]))
        for r in connection.queryAll(q):
            yield r[0],set(map(int,r[1].split(',')))

    @staticmethod
    def gene2groups(type=None,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [ProteinGroup.q.id,
             func.GROUP_CONCAT(func.DISTINCT(Peptide2PeptideGroup.q.groupID))],
            AND(Peptide2PeptideGroup.q.peptideID == Alignment.q.peptideID,
                Alignment.q.proteinID == Protein2ProteinGroup.q.proteinID,
                Protein2ProteinGroup.j.group,
                ProteinGroup.q.type=='Gene'),
            groupBy=[ProteinGroup.q.id]))
        for r in connection.queryAll(q):
            yield r[0],set(map(int,r[1].split(',')))

    @staticmethod
    def group2proteins(type=None,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Peptide2PeptideGroup.q.groupID,
             func.GROUP_CONCAT(func.DISTINCT(Alignment.q.proteinID))],
            AND(Peptide2PeptideGroup.q.peptideID == Alignment.q.peptideID),
            groupBy=[Peptide2PeptideGroup.q.groupID]))
        for r in connection.queryAll(q):
            yield r[0],set(map(int,r[1].split(',')))

    @staticmethod
    def group2genes(type=None,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Peptide2PeptideGroup.q.groupID,
             func.GROUP_CONCAT(func.DISTINCT(ProteinGroup.q.id))],
            AND(Peptide2PeptideGroup.q.peptideID == Alignment.q.peptideID,
                Alignment.q.proteinID == Protein2ProteinGroup.q.proteinID,
                Protein2ProteinGroup.j.group,
                ProteinGroup.q.type=='Gene'),
            groupBy=[Peptide2PeptideGroup.q.groupID]))
        for r in connection.queryAll(q):
            yield r[0],set(map(int,r[1].split(',')))

    @staticmethod
    def groupprotein2peptides(type=None,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Peptide2PeptideGroup.q.groupID,Alignment.q.proteinID,
             func.GROUP_CONCAT(func.DISTINCT(Peptide2PeptideGroup.q.peptideID))],
            AND(Peptide2PeptideGroup.q.peptideID == Alignment.q.peptideID),
            groupBy=[Peptide2PeptideGroup.q.groupID,Alignment.q.proteinID]))
        for r in connection.queryAll(q):
            yield (r[0],r[1]),set(map(int,r[2].split(',')))

    @staticmethod
    def groupgene2peptides(type=None,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Peptide2PeptideGroup.q.groupID,ProteinGroup.q.id,
             func.GROUP_CONCAT(func.DISTINCT(Peptide2PeptideGroup.q.peptideID))],
            AND(Peptide2PeptideGroup.q.peptideID == Alignment.q.peptideID,
                Alignment.q.proteinID == Protein2ProteinGroup.q.proteinID,
                Protein2ProteinGroup.j.group,
                ProteinGroup.q.type=='Gene'),
            groupBy=[Peptide2PeptideGroup.q.groupID,ProteinGroup.q.id]))
        for r in connection.queryAll(q):
            yield (r[0],r[1]),set(map(int,r[2].split(',')))

    @staticmethod
    def group2peptides(type=None,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Peptide2PeptideGroup.q.groupID,
             func.GROUP_CONCAT(func.DISTINCT(Peptide2PeptideGroup.q.peptideID))],
            groupBy=[Peptide2PeptideGroup.q.groupID]))
        for r in connection.queryAll(q):
            yield r[0],set(map(int,r[1].split(',')))

    @staticmethod
    def getSpecCounts(prefix,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [func.GROUP_CONCAT(func.DISTINCT(PeptideIon.q.peptideID)),Spectrum.q.id],
            AND(PeptideSpectrumMatch.j.spectrum,PeptideSpectrumMatch.j.peptideIon),
            groupBy=[Spectrum.q.id]))
        counts=defaultdict(int)
        for r in connection.queryAll(q):
            pepids = sorted(map(int,r[0].split(',')))
            name = ':'.join(map(str,pepids))
            name = ":".join([prefix,hashlib.sha1(name.encode('utf8')).hexdigest().lower(),
                                    hashlib.md5(name.encode('utf8')).hexdigest().lower(),
                                    str(len(pepids))])
            counts[name] += 1
        for pgname,cnt in counts.items():
            pg = PeptideGroup.access(pgname,connection=connection)
            if pg:
                yield (pg.id,cnt)

    @staticmethod
    def getMinFDR(prefix,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [func.GROUP_CONCAT(func.DISTINCT(PeptideIon.q.peptideID)),Spectrum.q.id,func.MAX(PeptideSpectrumMatch.q.value)],
            AND(PeptideSpectrumMatch.j.spectrum,PeptideSpectrumMatch.j.peptideIon),
            groupBy=[Spectrum.q.id]))
        minfdr=defaultdict(lambda: 1e+20)
        for r in connection.queryAll(q):
            pepids = sorted(map(int,r[0].split(',')))
            name = ':'.join(map(str,pepids))
            name = ":".join([prefix,hashlib.sha1(name.encode('utf8')).hexdigest().lower(),
                                    hashlib.md5(name.encode('utf8')).hexdigest().lower(),
                                    str(len(pepids))])
            minfdr[name] = min(minfdr[name],r[2])
        for pgname,fdr in minfdr.items():
            pg = PeptideGroup.access(pgname,connection=connection)
            if pg:
                yield (pg.id,fdr)

    @staticmethod
    def getSpecMetrics(type,connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        if type != None:
            sb = PeptideGroup.selectBy(type=type,connection=connection)
        else:
            sb = PeptideGroup.select(connection=connection)
        for pg in sb:
            spectra = pg.spectra()
            minmaxval = 1e+20
            for s in spectra:
                maxval = max(list(map(attrgetter('value'),s.psms)))
                minmaxval = min(minmaxval,maxval)
            yield pg.id,len(spectra),minmaxval

    def addPeptide(self,pep):
        try:
            Peptide2PeptideGroup(peptide=pep,group=self,connection=pep._connection)
        except DuplicateEntryError:
            pass
        self.count = len(self.members)

    def removePeptide(self,pep):
        Peptide2PeptideGroup.deleteBy(peptide=pep,group=self,connection=pep._connection)
        self.count = len(self.members)

    @staticmethod
    @memoize()
    def insert(name,type=None,connection=None):
        try:
            return PeptideGroup(name=name,type=type,connection=connection).id
        except DuplicateEntryError:
            return PeptideGroup.byName(name,connection=connection).id

    @staticmethod
    def find(name,type=None,connection=None):
        try:
            g = PeptideGroup.byName(name,connection=connection)
        except SQLObjectNotFound:
            try:
                g = PeptideGroup(name=name,type=type,connection=connection)
            except DuplicateEntryError:
                g = PeptideGroup.byName(name,connection=connection)
        return g

    @staticmethod
    def access(name,connection=None):
        try:
            return PeptideGroup.byName(name,connection=connection)
        except SQLObjectNotFound:
            pass
        return None

    @staticmethod
    def clear(connection=None):
        PeptideGroup.deleteMany(None,connection=connection)

    @staticmethod
    def normalize(quiet=False,connection=None):
        # members related join
        orphanids = set(map(attrgetter('id'),
                             PeptideGroup.select(Peptide2PeptideGroup.q.id == None,
                                                 join=LEFTJOINOn(PeptideGroup,Peptide2PeptideGroup,
                                                                 Peptide2PeptideGroup.j.group),
                                                 connection=connection)))

        if len(orphanids) == 0:
            return 0

        PeptideGroup.deleteMany(IN(PeptideGroup.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d peptide groups"%(count), file=sys.stderr)
        return count



class Alignment(MetadataContainer):
    nextid = 1
    peptide = ForeignKey('Peptide',cascade=True)
    protein = ForeignKey('Protein',cascade=True)
    start = IntCol(default=-1) # NULL's are considered to mismatch!
    end = IntCol(default=-1)
    substitutionSites = MultipleJoin('SubstitutionSite')
    index1 = DatabaseIndex('peptide')
    index2 = DatabaseIndex('protein','peptide','start','end',unique=True)

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from alignment;
        """)[0]
        if not maxid:
            maxid = 0
        Alignment.nextid = maxid+1

    @staticmethod
    @memoize(idstr="Alignment",maxprecache=10000000,
             precache=lambda c,n: map(lambda r: ((r[0],r[1],r[2],r[3]),r[4]),
                                     c.queryAll(c.sqlrepr(Select([Alignment.q.peptideID,
                                                                  Alignment.q.proteinID,
                                                                  Alignment.q.start,
                                                                  Alignment.q.end,
                                                                  Alignment.q.id],orderBy='id',start=0,end=n)))))
    def insert(peptideid,proteinid,start=-1,end=-1,metadata=None,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = Alignment.nextid
            connection.query("""
                insert into alignment (id,peptide_id,protein_id,start,end,metadata) values (%d,%d,%d,%d,%d,%s)
            """%(id,peptideid,proteinid,start,end,
                 "\""+base64.encodebytes(pickle.dumps(metadata)).decode()+"\"" if metadata else 'NULL'))
            Alignment.nextid += 1
            return id
            # return Alignment(peptideID=peptideid,protein=protein,start=start,end=end,
            #                  metadata=metadata,connection=connection)
        except DuplicateEntryError:
            return Alignment.select(AND(Alignment.q.peptideID == peptideid,
                                        Alignment.q.proteinID == proteinid,
                                        Alignment.q.start == start,
                                        Alignment.q.end == end),
                                    connection=connection)[0].id

    @staticmethod
    def find(peptide,protein,start=None,end=None,connection=None):
        try:
            if start:
                a = Alignment.select(AND(Alignment.q.peptide == peptide,
                                         Alignment.q.protein == protein,
                                         Alignment.q.start == start,
                                         Alignment.q.end == end),
                                     connection=connection)[0]
            else:
                a = Alignment.select(AND(Alignment.q.peptide == peptide,
                                         Alignment.q.protein == protein),
                                     connection=connection)[0]
        except IndexError:
            try:
                a = Alignment(peptide=peptide,protein=protein,start=start,end=end,
                              connection=connection)
            except DuplicateEntryError:
                if start:
                    a = Alignment.select(AND(Alignment.q.peptide == peptide,
                                             Alignment.q.protein == protein,
                                             Alignment.q.start == start,
                                             Alignment.q.end == end),
                                         connection=connection)[0]
                else:
                    a = Alignment.select(AND(Alignment.q.peptide == peptide,
                                             Alignment.q.protein == protein),
                                         connection=connection)[0]
        return a

    @staticmethod
    def normalize(quiet=False,connection=None):
        orphanids = set(map(attrgetter('id'),
                             Alignment.select(Peptide.q.id == None,
                                              join=LEFTJOINOn(Alignment,Peptide,Alignment.j.peptide),
                                              connection=connection)))
        # print >>sys.stderr, len(orphanids)
        orphanids1 = set(map(attrgetter('id'),
                              Alignment.select(Protein.q.id == None,
                                               join=LEFTJOINOn(Alignment,Protein,Alignment.j.protein),
                                               connection=connection)))
        # print >>sys.stderr, len(orphanids1)
        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        Alignment.deleteMany(IN(Alignment.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d alignments"%(count), file=sys.stderr)
        return count

class SubstitutionSite(MetadataContainer):
    position = IntCol()
    alignment = ForeignKey('Alignment',cascade=True)
    substitution = ForeignKey('Substitution',cascade=True)
    index1 = DatabaseIndex('alignment','position','substitution',unique=True)
    index2 = DatabaseIndex('substitution')

    def peptide(self):
        return self.alignment.peptide

    def altpeptide(self):
        pepseq = self.alignment.peptide.sequence
        pos = self.position - 1
        delpepseq = self.substitution.peptideAminoAcids
        insproseq = self.substitution.proteinAminoAcids
        newseq = pepseq[:pos]+insproseq+pepseq[pos+len(delpepseq):]
        return newseq

    @staticmethod
    @memoize()
    def insert(alignmentid, position, substitutionid, connection=None):
        try:
            return SubstitutionSite(alignmentID=alignmentid,position=position,substitutionID=substitutionid,
                                    connection=connection)
        except DuplicateEntryError:
            return SubstitutionSite.select(AND(SubstitutionSite.q.alignmentid == alignmentid,
                                               SubstitutionSite.q.position == position,
                                               SubstitutionSite.q.substitutionid == substitutionid),
                                           connection=connection)[0]

    @staticmethod
    def find(alignment, position, substitution, connection=None):
        try:
            ms = SubstitutionSite.select(AND(SubstitutionSite.q.alignment == alignment,
                                             SubstitutionSite.q.position == position,
                                             SubstitutionSite.q.substitution == substitution),
                                         connection=connection)[0]
        except IndexError:
            try:
                ms = SubstitutionSite(alignment=alignment,position=position,substitution=substitution,
                                      connection=connection)
            except DuplicateEntryError:
                ms = SubstitutionSite.select(AND(SubstitutionSite.q.alignment == alignment,
                                                 SubstitutionSite.q.position == position,
                                                 SubstitutionSite.q.substitution == substitution),
                                             connection=connection)[0]
        return ms

    @staticmethod
    def normalize(quiet=False,connection=None):
        # alignment foreign key constraint
        orphanids = set(map(attrgetter('id'),
                             SubstitutionSite.select(Alignment.q.id == None,
                                                     join=LEFTJOINOn(SubstitutionSite,Alignment,
                                                                     SubstitutionSite.j.alignment),
                                                     connection=connection)))
        # modification foreign key constraint
        orphanids1 = set(map(attrgetter('id'),
                              SubstitutionSite.select(Substitution.q.id == None,
                                                      join=LEFTJOINOn(SubstitutionSite,Substitution,
                                                                      SubstitutionSite.j.substitution),
                                                      connection=connection)))

        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        SubstitutionSite.deleteMany(IN(SubstitutionSite.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d substitution sites"%(count), file=sys.stderr)
        return count


class Substitution(MetadataContainer):
    peptideAminoAcids = StringCol()
    proteinAminoAcids = StringCol()
    massDelta = FloatCol()
    index1 = DatabaseIndex('peptideAminoAcids','proteinAminoAcids',unique=True)
    sites = MultipleJoin('SubstitutionSite')

    def alignments(self,connection=None):
        lastal = None
        for ss in SubstitutionSite.select(SubstitutionSite.q.substitution==self,
                                           orderBy=SubstitutionSite.q.alignment,
                                           connection=connection):
            if ss.alignmentID != lastal:
                yield ss.alignment
                lastal = ss.alignmentID

    @staticmethod
    @memoize(key=lambda a,k: (a[0], a[1]))
    def insert(pepaas,proaas,massdelta,connection=None):
        try:
            return Substitution(peptideAminoAcids=pepaas,proteinAminoAcids=proaas,
                                massDelta=massdelta,connection=connection)
        except DuplicateEntryError:
            return Substitution.select(AND(Substitution.q.peptideAminoAcids==pepaas,
                                           Substitution.q.proteinAminoAcids==proaas),
                                       connection=connection)[0]

    @staticmethod
    def find(pepaas,proaas,connection=None):
        try:
            s = Substitution.select(AND(Substitution.q.peptideAminoAcids==pepaas,
                                        Substitution.q.proteinAminoAcids==proaas),
                                    connection=connection)[0]
        except IndexError:
            try:
                s = Substitution(peptideAminoAcids=pepaas,proteinAminoAcids=proaas,
                                 connection=connection)
            except DuplicateEntryError:
                s = Substitution.select(AND(Substitution.q.peptideAminoAcids==pepaas,
                                            Substitution.q.proteinAminoAcids==proaas),
                                        connection=connection)[0]
        return s

    @staticmethod
    def normalize(quiet=False,connection=None):
        # substitution site foreign key constraint
        orphanids = set(map(attrgetter('id'),
                             Substitution.select(SubstitutionSite.q.id == None,
                                                 join=LEFTJOINOn(Substitution,SubstitutionSite,
                                                                 SubstitutionSite.j.substitution),
                                                 connection=connection)))

        if len(orphanids) == 0:
            return 0

        Substitution.deleteMany(IN(Substitution.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d substitutions"%(count), file=sys.stderr)
        return count


class Protein(MetadataContainer):
    nextid = 1
    accession = StringCol(alternateID=True)
    decoy = BoolCol(default=False)
    index1 = DatabaseIndex(dict(column=accession,length=50))
    alignments = MultipleJoin('Alignment',joinColumn='protein_id')
    membership = MultipleJoin('Protein2ProteinGroup',joinColumn='protein_id')

    def short_accession(self):
        return self.getdata('short_accession')

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from protein;
        """)[0]
        if not maxid:
            maxid = 0
        Protein.nextid = maxid+1

    def peptides(self,locus=None,loci=None):
        if locus != None:
            return [a.peptide for a in [a for a in self.alignments if a.start <= locus < a.end]]
        if loci != None:
            return [a.peptide for a in [a for a in self.alignments if a.start <= min(loci) and max(loci) < a.end]]
        return [a.peptide for a in self.alignments]

    def alignmentsto(self,peptide):
        return [a for a in self.alignments if a.peptide == peptide]

    def unshared_peptides(self):
        unsh = set()
        for pep in self.peptides():
            if len(pep.proteins()) == 1:
                unsh.add(pep)
        return unsh

    def peptidegroups(self,type):
        pgs = []
        for pep in self.peptides():
            pgs.extend(pep.groups())
        return pgs

    def groups(self):
        return list(m.group for m in self.membership)

    @staticmethod
    @memoize(idstr="Protein",
             key=lambda a,k: a[0],
             precache=lambda c,n: c.queryAll(c.sqlrepr(Select([Protein.q.accession,
                                                               Protein.q.id],orderBy='id',start=0,end=n))))
    def insert(protein,decoy=None,metadata=None,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = Protein.nextid
            connection.query("""
                insert into protein (id,accession,decoy,metadata) values (%s,"%s",%s,%s)
            """%(id,protein,decoymap(decoy),
                "\""+base64.encodebytes(pickle.dumps(metadata)).decode()+"\"" if metadata else 'NULL'))
            Protein.nextid += 1
            return id
            # return Protein(accession=protein,metadata=metadata,connection=connection)
        except DuplicateEntryError:
            return Protein.select(Protein.q.accession==protein,
                                  connection=connection)[0].id

    @staticmethod
    def find(protein,decoy=None,connection=None):
        try:
            p = Protein.select(Protein.q.accession==protein,
                               connection=connection)[0]
        except IndexError:
            try:
                p = Protein(accession=protein, decoy=decoy, connection=connection)
            except DuplicateEntryError:
                p = Protein.select(Protein.q.accession==protein,
                                   connection=connection)[0]
        return p

    @staticmethod
    def access(protein,connection=None):
        try:
            return Protein.select(Protein.q.accession==protein,
                                  connection=connection)[0]
        except IndexError:
            pass
        return None

    @staticmethod
    def clear(connection=None):
        Protein.deleteMany(None,connection=connection)

    @staticmethod
    def normalize(quiet=False,connection=None):
        orphanids = set(map(attrgetter('id'),
                             Protein.select(Alignment.q.id == None,
                                            join=LEFTJOINOn(Protein,Alignment,Alignment.j.protein),
                                            connection=connection)))
        if len(orphanids) == 0:
            return 0

        Protein.deleteMany(IN(Protein.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d proteins"%(count), file=sys.stderr)
        return count

class Protein2ProteinGroup(SQLObject):
    nextid = 1
    protein = ForeignKey('Protein',cascade=True)
    group = ForeignKey('ProteinGroup',cascade=True)
    index1 = DatabaseIndex('protein','group',unique=True)

    @staticmethod
    def initid(connection):
        maxid = connection.queryOne("""
            select max(id) from protein2_protein_group;
        """)[0]
        if not maxid:
            maxid = 0
        Protein2ProteinGroup.nextid = maxid+1

    @staticmethod
    @memoize()
    def insert(proteinid,groupid,connection=None):
        try:
            if not connection:
                connection = sqlhub.getConnection()
            id = Protein2ProteinGroup.nextid
            connection.query("""
                insert into protein2_protein_group (id,protein_id,group_id) values (%d,%d,%d)
            """%(id,proteinid,groupid))
            Protein2ProteinGroup.nextid += 1
            return id
            #return Protein2ProteinGroup(proteinID=proteinid,groupID=groupid,connection=connection).id
        except DuplicateEntryError:
            return Protein2ProteinGroup.selectBy(proteinID=proteinid,groupID=groupid,connection=connection)[0].id

    @staticmethod
    def normalize(quiet=False,connection=None):
        # protein foreign key
        orphanids = set(map(attrgetter('id'),
                             Protein2ProteinGroup.select(Protein.q.id == None,
                                                         join=LEFTJOINOn(Protein2ProteinGroup,Protein,
                                                                         Protein2ProteinGroup.j.protein),
                                                         connection=connection)))

        # protein group foreign key
        orphanids1 = set(map(attrgetter('id'),
                              Protein2ProteinGroup.select(ProteinGroup.q.id == None,
                                                          join=LEFTJOINOn(Protein2ProteinGroup,ProteinGroup,
                                                                          Protein2ProteinGroup.j.group),
                                                          connection=connection)))

        orphanids |= orphanids1

        if len(orphanids) == 0:
            return 0

        Protein2ProteinGroup.deleteMany(IN(Protein2ProteinGroup.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d protein group memberships"%(count), file=sys.stderr)
        return count

class ProteinGroup(MetadataContainer):
    name = StringCol(alternateID=True)
    type = StringCol(default=None)
    members = MultipleJoin('Protein2ProteinGroup',joinColumn='group_id')

    def proteins(self):
        return list(m.protein for m in self.members)

    @staticmethod
    def distinct_peptides(proteins,connection=None):
        prids = list(map(attrgetter('id'),proteins))
        return Peptide.select(AND(IN(Protein.q.id,prids),
                                  Alignment.j.peptide,Alignment.j.protein),
                              distinct=True,
                              connection=connection)

    @staticmethod
    def unshared_peptides(proteins,connection=None):
        prids = set(map(attrgetter('id'),proteins))
        unshpeps = []
        for pep in ProteinGroup.distinct_peptides(proteins,connection):
            if len(set(map(attrgetter('id'),pep.proteins()))&prids) == 1:
                unshpeps.append(pep)
        return unshpeps

    def peptides(self):
        return ProteinGroup.distinct_peptides((m.protein for m in self.members),self._connection)

    def component(self,type=None):
        comp = set()
        comppeps = set()
        compprots = set()
        prgs = set([self])
        while len(prgs) > 0:
            prg = prgs.pop()
            comp.add(prg)
            peps = prg.peptides()
            for pep in peps:
                comppeps.add(pep)
                for pr in pep.proteins():
                    compprots.add(pr)
                    for prg in pr.groups():
                        if type and prg.type != type:
                            continue
                        if prg not in comp and prg not in prgs:
                            prgs.add(prg)
        return comp,compprots,comppeps

    def addProtein(self,pr):
        try:
            Protein2ProteinGroup(protein=pr,group=self,connection=pr._connection)
        except DuplicateEntryError:
            pass

    def removeProtein(self,pr):
        Protein2ProteinGroup.deleteBy(protein=pr,group=self,connection=pr._connection)

    @staticmethod
    @memoize()
    def insert(name,type=None,connection=None):
        try:
            return ProteinGroup(name=name,type=type,connection=connection).id
        except DuplicateEntryError:
            return ProteinGroup.byName(name,connection=connection).id

    @staticmethod
    def find(name,type=None,connection=None):
        try:
            g = ProteinGroup.byName(name,connection=connection)
        except SQLObjectNotFound:
            try:
                g = ProteinGroup(name=name,type=type,connection=connection)
            except DuplicateEntryError:
                g = ProteinGroup.byName(name,connection=connection)
        return g

    @staticmethod
    def access(name,connection=None):
        try:
            return ProteinGroup.byName(name,connection=connection)
        except SQLObjectNotFound:
            pass
        return None

    @staticmethod
    def normalize(quiet=False,connection=None):
        # members related join
        orphanids = set(map(attrgetter('id'),
                             ProteinGroup.select(Protein2ProteinGroup.q.id == None,
                                                 join=LEFTJOINOn(ProteinGroup,Protein2ProteinGroup,
                                                                 Protein2ProteinGroup.j.group),
                                                 connection=connection)))

        if len(orphanids) == 0:
            return 0

        ProteinGroup.deleteMany(IN(ProteinGroup.q.id,orphanids),connection=connection)

        count = len(orphanids)
        if not quiet:
            print("Removed %d protein groups"%(count), file=sys.stderr)
        return count

class FileProteinPSMCounts(SQLObject):
    file = ForeignKey('SpectrumFile')
    protein = ForeignKey('Protein')
    count = IntCol(default=0)
    index1 = DatabaseIndex('file','protein',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            FileProteinPSMCounts.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def counts(connection=None):
        return defaultdict(int,map(lambda r: ((r.fileID,r.proteinID),r.count),FileProteinPSMCounts.select(connection=connection)))

    @staticmethod
    def clear(connection=None):
        FileProteinPSMCounts.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select(
            [Spectrum.q.fileID,Alignment.q.proteinID,func.COUNT(func.DISTINCT(PeptideSpectrumMatch.q.spectrumID))],
            AND(Alignment.j.peptide,PeptideIon.j.peptide,PeptideSpectrumMatch.j.peptideIon,
                PeptideSpectrumMatch.j.spectrum),
            groupBy=[Spectrum.q.fileID,Alignment.q.proteinID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(FileProteinPSMCounts.sqlmeta.table,
                                                       values=dict(list(zip(['file_id','protein_id','count'],r))))))


class PeptideProteinCounts(SQLObject):
    peptide = ForeignKey('Peptide')
    count = IntCol(default=0)
    index1 = DatabaseIndex('peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        try:
            PeptideProteinCounts.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def clear(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        PeptideProteinCounts.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        PeptideProteinCounts.clear(connection)
        q = connection.sqlrepr(Select(
            [Alignment.q.peptideID,func.COUNT(func.DISTINCT(Alignment.q.proteinID))],
            groupBy=[Alignment.q.peptideID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(PeptideProteinCounts.sqlmeta.table,
                                                       values=dict(list(zip(['peptide_id','count'],r))))))

    @staticmethod
    def unshared_peptides(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select([PeptideProteinCounts.q.peptideID],PeptideProteinCounts.q.count == 1))
        return set(map(itemgetter(0),connection.queryAll(q)))

class PeptideGeneCounts(SQLObject):
    peptide = ForeignKey('Peptide')
    count = IntCol(default=0)
    index1 = DatabaseIndex('peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        try:
            PeptideGeneCounts.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def clear(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        PeptideGeneCounts.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        PeptideGeneCounts.clear(connection)
        q = connection.sqlrepr(Select(
            [Alignment.q.peptideID,func.COUNT(func.DISTINCT(Protein2ProteinGroup.q.groupID))],
            AND(Alignment.q.proteinID == Protein2ProteinGroup.q.proteinID,
                Protein2ProteinGroup.j.group,
                ProteinGroup.q.type == 'Gene'),
            groupBy=[Alignment.q.peptideID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(PeptideGeneCounts.sqlmeta.table,
                                                       values=dict(list(zip(['peptide_id','count'],r))))))

    @staticmethod
    def unshared_peptides(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select([PeptideGeneCounts.q.peptideID],PeptideGeneCounts.q.count == 1))
        return set(map(itemgetter(0),connection.queryAll(q)))

class PeptidePSMCounts(SQLObject):
    peptide = ForeignKey('Peptide')
    count = IntCol(default=0)
    index1 = DatabaseIndex('peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            PeptidePSMCounts.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def counts(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select([PeptidePSMCounts.q.peptideID,PeptidePSMCounts.q.count]))
        # start = time.time()
        d = defaultdict(int,connection.queryAll(q))
        # print >>sys.stderr, "Elapsed: ", time.time()-start
        return d

    @staticmethod
    def counts1(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        # start = time.time()
        q = connection.sqlrepr(Select([PeptideIon.q.peptideID,PeptideSpectrumMatch.q.spectrumID],
                                      PeptideSpectrumMatch.j.peptideIon))
        d = defaultdict(set)
        for r in connection.queryAll(q):
            d[r[0]].add(r[1])
        r = defaultdict(int)
        for k,v in d.items():
            r[k] = len(v)
        # print >>sys.stderr, "Elapsed: ", time.time()-start
        return r

    @staticmethod
    def clear(connection=None):
        PeptidePSMCounts.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        PeptidePSMCounts.clear(connection)
        q = connection.sqlrepr(Select(
            [PeptideIon.q.peptideID,func.COUNT(func.DISTINCT(PeptideSpectrumMatch.q.spectrumID))],
            PeptideSpectrumMatch.j.peptideIon,
            groupBy=[PeptideIon.q.peptideID]))
        first = True;
        start = time.time()
        for r in connection.queryAll(q):
            if first:
                # print >>sys.stderr, "Elapsed: ", time.time()-start
                start = time.time()
                first = False
            connection.query(connection.sqlrepr(Insert(PeptidePSMCounts.sqlmeta.table,
                                                       values=dict(list(zip(['peptide_id','count'],r))))))
        # print >>sys.stderr, "Elapsed: ", time.time()-start

    @staticmethod
    def compute1(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select([PeptideIon.q.peptideID,PeptideSpectrumMatch.q.spectrumID],
                                       PeptideSpectrumMatch.j.peptideIon))
        d=defaultdict(set)
        start = time.time()
        for r in connection.queryAll(q):
            d[r[0]].add(r[1])
        start = time.time()
        for k,v in d.items():
            connection.query(connection.sqlrepr(Insert(PeptidePSMCounts.sqlmeta.table,
                                                       values=dict(list(zip(['peptide_id','count'],(k,len(v))))))))

class PeptidePSMMinDelta:
    @staticmethod
    def minvalues1(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        # start = time.time()
        q = connection.sqlrepr(Select([PeptideIon.q.peptideID,PeptideSpectrumMatch.q.delta],
                                       PeptideSpectrumMatch.j.peptideIon))
        d = defaultdict(set)
        for r in connection.queryAll(q):
            d[r[0]].add(r[1])
        r = defaultdict(lambda: 1e+20)
        for k,v in d.items():
            r[k] = min(v)
        # print >>sys.stderr, "Elapsed: ", time.time()-start
        return r


class PeptidePSMMinValue(SQLObject):
    peptide = ForeignKey('Peptide')
    minvalue = FloatCol(default=1e+20)
    index1 = DatabaseIndex('peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            PeptidePSMMinValue.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def minvalues(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        q = connection.sqlrepr(Select([PeptidePSMMinValue.q.peptideID,PeptidePSMMinValue.q.minvalue]))
        return defaultdict(lambda: 1e+20,connection.queryAll(q))

    @staticmethod
    def minvalues1(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        # start = time.time()
        q = connection.sqlrepr(Select([PeptideIon.q.peptideID,PeptideSpectrumMatch.q.value],
                                       PeptideSpectrumMatch.j.peptideIon))
        d = defaultdict(set)
        for r in connection.queryAll(q):
            d[r[0]].add(r[1])
        r = defaultdict(lambda: 1e+20)
        for k,v in d.items():
            r[k] = min(v)
        # print >>sys.stderr, "Elapsed: ", time.time()-start
        return r

    @staticmethod
    def clear(connection=None):
        PeptidePSMMinValue.deleteMany(None,connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        PeptidePSMMinValue.clear(connection)
        q = connection.sqlrepr(Select(
            [PeptideIon.q.peptideID,func.MIN(PeptideSpectrumMatch.q.value)],
            PeptideSpectrumMatch.j.peptideIon,
            groupBy=[PeptideIon.q.peptideID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(PeptidePSMMinValue.sqlmeta.table,
                                                       values=dict(list(zip(['peptide_id','minvalue'],r))))))

class FilePeptidePSMMinValue(SQLObject):
    file = ForeignKey('SpectrumFile')
    peptide = ForeignKey('Peptide')
    minvalue = FloatCol(default=1e+20)
    index1 = DatabaseIndex('file','peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            FilePeptidePSMMinValue.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def minvalues(connection=None):
        return defaultdict(lambda: 1e+20,map(lambda r: ((r.fileID,r.peptideID),r.minvalue),FilePeptidePSMMinValue.select(connection=connection)))

    @staticmethod
    def clear(connection=None):
        FilePeptidePSMMinValue.deleteMany(None,connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        FilePeptidePSMMinValue.clear(connection)
        q = connection.sqlrepr(Select(
            [Spectrum.q.fileID,PeptideIon.q.peptideID,func.MIN(PeptideSpectrumMatch.q.value)],
            AND(PeptideSpectrumMatch.j.peptideIon,
                PeptideSpectrumMatch.j.spectrum),
            groupBy=[Spectrum.q.fileID,PeptideIon.q.peptideID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(FilePeptidePSMMinValue.sqlmeta.table,
                                                       values=dict(list(zip(['file_id','peptide_id','minvalue'],r))))))

class FilePeptidePSMCounts(SQLObject):
    file = ForeignKey('SpectrumFile')
    peptide = ForeignKey('Peptide')
    count = IntCol(default=0)
    index1 = DatabaseIndex('file','peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            FilePeptidePSMCounts.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def counts(connection=None):
        return defaultdict(int,map(lambda r: ((r.fileID,r.peptideID),r.count),FilePeptidePSMCounts.select(connection=connection)))

    @staticmethod
    def clear(connection=None):
        FilePeptidePSMCounts.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        FilePeptidePSMCounts.clear(connection)
        q = connection.sqlrepr(Select(
            [Spectrum.q.fileID,PeptideIon.q.peptideID,func.COUNT(func.DISTINCT(PeptideSpectrumMatch.q.spectrumID))],
            AND(PeptideSpectrumMatch.j.peptideIon,
                PeptideSpectrumMatch.j.spectrum),
            groupBy=[Spectrum.q.fileID,PeptideIon.q.peptideID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(FilePeptidePSMCounts.sqlmeta.table,
                                                       values=dict(list(zip(['file_id','peptide_id','count'],r))))))

class FilePeptideIonCounts(SQLObject):
    file = ForeignKey('SpectrumFile')
    peptide = ForeignKey('Peptide')
    count = IntCol(default=0)
    index1 = DatabaseIndex('file','peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            FilePeptideIonCounts.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def counts(connection=None):
        return defaultdict(int,map(lambda r: ((r.fileID,r.peptideID),r.count),FilePeptideIonCounts.select(connection=connection)))

    @staticmethod
    def clear(connection=None):
        FilePeptideIonCounts.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        FilePeptideIonCounts.clear(connection)
        q = connection.sqlrepr(Select(
            [Spectrum.q.fileID,PeptideIon.q.peptideID,func.COUNT(func.DISTINCT(PeptideIon.q.id))],
            AND(PeptideSpectrumMatch.j.peptideIon,
                PeptideSpectrumMatch.j.spectrum),
            groupBy=[Spectrum.q.fileID,PeptideIon.q.peptideID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(FilePeptideIonCounts.sqlmeta.table,
                                                       values=dict(list(zip(['file_id','peptide_id','count'],r))))))

class FileIonPSMCounts(SQLObject):
    file = ForeignKey('SpectrumFile')
    ion = ForeignKey('Protein')
    count = IntCol(default=0)
    index1 = DatabaseIndex('file','ion',unique=True)

    @staticmethod
    def clear(connection=None):
        FileIonPSMCounts.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        FileIonPSMCounts.clear(connection)
        q = connection.sqlrepr(Select(
            [Spectrum.q.fileID,PeptideSpectrumMatch.q.peptideIonID,
             func.COUNT(func.DISTINCT(PeptideSpectrumMatch.q.spectrumID))],
            PeptideSpectrumMatch.j.spectrum,
            groupBy=[Spectrum.q.fileID,PeptideSpectrumMatch.q.peptideIonID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(FileIonPSMCounts.sqlmeta.table,
                                                       values=dict(list(zip(['file_id','ion_id','count'],r))))))

class FilePeptideSpectra(SQLObject):
    file = ForeignKey('SpectrumFile')
    peptide = ForeignKey('Peptide')
    spectra = StringCol(length=32*1024,default=0)
    index1 = DatabaseIndex('file','peptide',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            FilePeptideSpectra.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def clear(connection=None):
        FilePeptideSpectra.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        FilePeptideSpectra.clear(connection)
        q = connection.sqlrepr(Select(
            [Spectrum.q.fileID,PeptideIon.q.peptideID,func.GROUP_CONCAT(func.DISTINCT(Spectrum.q.id))],
            AND(PeptideSpectrumMatch.j.peptideIon,
                PeptideSpectrumMatch.j.spectrum),
            groupBy=[Spectrum.q.fileID,PeptideIon.q.peptideID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(FilePeptideSpectra.sqlmeta.table,
                                                       values=dict(list(zip(['file_id','peptide_id','spectra'],r))))))
class FilePeptideIonSpectra(SQLObject):
    file = ForeignKey('SpectrumFile')
    peptideion = ForeignKey('PeptideIon')
    spectra = StringCol(length=32*1024,default=0)
    index1 = DatabaseIndex('file','peptideion',unique=True)

    @staticmethod
    def empty(connection=None):
        try:
            FilePeptideIonSpectra.select(connection=connection)[0]
        except IndexError:
            return True
        return False

    @staticmethod
    def clear(connection=None):
        FilePeptideIonSpectra.deleteMany(None,connection=connection)

    @staticmethod
    def compute(connection=None):
        if not connection:
            connection = sqlhub.getConnection()
        FilePeptideIonSpectra.clear(connection)
        q = connection.sqlrepr(Select(
            [Spectrum.q.fileID,PeptideSpectrumMatch.j.peptideIonID,
             func.GROUP_CONCAT(func.DISTINCT(Spectrum.q.id))],
            AND(PeptideSpectrumMatch.j.peptideIon,
                PeptideSpectrumMatch.j.spectrum),
            groupBy=[Spectrum.q.fileID,PeptideSpectrumMatch.j.peptideIonID]))
        for r in connection.queryAll(q):
            connection.query(connection.sqlrepr(Insert(FilePeptideIonSpectra.sqlmeta.table,
                                                       values=dict(list(zip(['file_id','peptideion_id','spectra'],r))))))

tables = list(map(eval,list(map(str.strip,"""
Analysis
SpectrumFileGroup
Sample
Label
SampleSpectrumFileLabel
SpectrumFile
Spectrum
PeptideSpectrumMatch
PeptideIon
ModificationSite
Modification
Peptide
Peptide2PeptideGroup
PeptideGroup
Alignment
SubstitutionSite
Substitution
Protein
Protein2ProteinGroup
ProteinGroup
PeptidePSMCounts
PeptideProteinCounts
PeptideGeneCounts
PeptidePSMMinValue
FileProteinPSMCounts
FilePeptidePSMCounts
FilePeptidePSMMinValue
FilePeptideIonCounts
FileIonPSMCounts
FilePeptideSpectra
FilePeptideIonSpectra
""".split()))))

initTables = list(map(eval,list(map(str.strip,"""
Spectrum
PeptideSpectrumMatch
PeptideIon
Peptide
Peptide2PeptideGroup
Protein
Protein2ProteinGroup
Alignment
ModificationSite
""".split()))))

def dump_counts(fh=sys.stderr,connection=None):
    for t in tables:
        c = t.select(connection=connection).count()
        if c > 0:
            print("%25s: %d"%(t.__name__,c), file=fh)

def normalize(connection,alignedPeptidesOnly=True,quiet=False,vacuum=False):
    # dump_counts(connection=connection)
    while True:
        count = 0
        count += Spectrum.normalize(quiet=quiet,connection=connection)
        count += PeptideSpectrumMatch.normalize(quiet=quiet,connection=connection)
        count += PeptideIon.normalize(quiet=quiet,connection=connection)
        count += ModificationSite.normalize(quiet=quiet,connection=connection)
        count += Peptide.normalize(alignedOnly=alignedPeptidesOnly,quiet=quiet,connection=connection)
        count += Peptide2PeptideGroup.normalize(quiet=quiet,connection=connection)
        count += PeptideGroup.normalize(quiet=quiet,connection=connection)
        count += Alignment.normalize(quiet=quiet,connection=connection)
        count += SubstitutionSite.normalize(quiet=quiet,connection=connection)
        count += Substitution.normalize(quiet=quiet,connection=connection)
        count += Protein.normalize(quiet=quiet,connection=connection)
        count += Protein2ProteinGroup.normalize(quiet=quiet,connection=connection)
        count += ProteinGroup.normalize(quiet=quiet,connection=connection)
        if count == 0:
            break
    PeptideProteinCounts.clear(connection=connection)
    PeptideGeneCounts.clear(connection=connection)
    PeptidePSMCounts.clear(connection=connection)
    PeptidePSMMinValue.clear(connection=connection)
    FileProteinPSMCounts.clear(connection=connection)
    FilePeptidePSMCounts.clear(connection=connection)
    FilePeptidePSMMinValue.clear(connection=connection)
    FilePeptideIonCounts.clear(connection=connection)
    FilePeptideSpectra.clear(connection=connection)
    FilePeptideIonSpectra.clear(connection=connection)
    FileIonPSMCounts.clear(connection=connection)
    Peptide2PeptideGroup.clear(connection=connection)
    PeptideGroup.clear(connection=connection)
    # dump_counts(connection=connection)
    if vacuum:
        connection.query('VACUUM;')
