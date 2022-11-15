
from . import model
from .load import SearchResultLoad
from collections import defaultdict
import os, os.path, sys, re
from operator import attrgetter, itemgetter
import itertools

from peptidescan.PeptideRemapper import PeptideRemapper

class PSMDb:

    def __init__(self,filename=None,quiet=False,**kw):

        kw['processConnection'] = False
        kw['threadConnection'] = False
        kw['filename'] = filename
        # kw['debug'] = True
        self.quiet = quiet
        self.tofile = False
        if filename and os.path.exists(filename) and kw.get('memory'):
            kw['fromfile'] = True
        if filename and kw.get('memory') and not kw.get('readonly',True):
            self.tofile = True
            try:
                del kw['readonly']
            except KeyError:
                pass

        self.conn,self.filename = model.init(**kw)

        self.todelete = []

        if not filename:
            self.todelete.append(self.filename)

        self.model = model

    def __del__(self):
        if self.tofile:
            if os.path.exists(self.filename):
                os.unlink(self.filename)
            model.flush(self.conn,self.filename)
        self.conn.close()
        for f in self.todelete:
            try:
                os.unlink(f)
            except:
                continue

    def clone(self,**kw):
        assert hasattr(self,'filename')
        kw['clone'] = self.filename
        return PSMDb(**kw)

    def begin(self,nopragma=False):
        trans = self.conn.transaction()
        pragmas = """
        PRAGMA journal_mode=OFF;
        PRAGMA automatic_indexing=OFF;
        PRAGMA cache_size=30000;
        PRAGMA temp_store=MEMORY;
        PRAGMA page_size=32768;
        PRAGMA synchronous=OFF;
        PRAGMA foreign_key=ON;
        """.split('\n')
        if nopragma:
            pragmas = ""
        for pr in pragmas:
            if pr:
                trans.query(pr)
        return trans

    def load(self,filename,psmfilter=None,noprotein=False,accrule=None,reload=False,analysisname=None,analysisroot=None,normalize=True):
        if not psmfilter:
            from .filter import NoFilter
            psmfilter = NoFilter()
        trans = self.begin()
        if SearchResultLoad(filename,psmfilter,noprotein,accrule,reload,analysisname,analysisroot,connection=trans) and normalize:
            model.normalize(trans,alignedPeptidesOnly=(not noprotein),quiet=self.quiet,vacuum=False)
        trans.commit(close=True)

    def normalize(self,noprotein=False,vacuum=False):
        trans = self.begin(nopragma=vacuum)
        model.normalize(trans,alignedPeptidesOnly=(not noprotein),quiet=self.quiet,vacuum=vacuum)
        trans.commit(close=True)

    def analyses(self):
        return model.Analysis.select(connection=self.conn)

    def spectrumfiles(self):
        return model.SpectrumFile.select(connection=self.conn)

    def spectrumfileid(self,name):
        return model.SpectrumFile.access(name,connection=self.conn)

    def spectrumfileid_name(self,id):
        return model.SpectrumFile.get(id,connection=self.conn).name

##     def newSpectrumFileGroup(self,name,spectrumfileids,type=None):
##         sfgrp = model.SpectrumFileGroup.find(name,type=type,connection=self.conn)
##         for sfid in spectrumfileids:
##             sfgrp.addSpectrumFile(model.SpectrumFile.get(sfid,connection=self.conn))
##         return sfgrp

##     def spectrumfilegroupid_name(self,id):
##      return model.SpectrumFileGroup.get(id,connection=self.conn).name

##     def spectrumfilegroups(self):
##         return model.SpectrumFileGroup.select(connection=self.conn)

##     def spectrumfilegroup(self,id):
##         return model.SpectrumFileGroup.get(id,connection=self.conn)

    def newSample(self,name,spectrumfiles,labels=None,type=None):
        s = model.Sample.find(name,type=type,connection=self.conn)
        if labels == None:
            for sf in spectrumfileids:
                model.SampleSpectrumFileLabel.insert(s,sf,connection=self.conn)
        else:
            assert len(spectrumfiles) == len(labels)
            for sf,labstr in zip(spectrumfiles,labels):
                lab = model.Label.find(labstr,labstr.split('-')[0],connection=self.conn)
                model.SampleSpectrumFileLabel.insert(s,sf,lab,connection=self.conn)
        return s

    def deleteSamples(self,type=None):
        q = model.NOT(False)
        if type:
            q = (model.Sample.q.type==type)
        for sfgrp in model.Sample.select(q,connection=self.conn):
            sfgrp.destroySelf()

    def sample_name(self,id):
        model.Sample.get(id,connection=self.conn).name

    def samples(self):
        return model.Sample.select(connection=self.conn)

    def sample(self,id):
        return model.Sample.get(id,connection=self.conn)

    def analytical_samples(self):
        analsamps = set()
        for sf in self.spectrumfiles():
            analsamps.add(sf.analytical_sample())
        return analsamps

    def deleteNoSampleSpectrumFiles(self):
        toremove=set()
        for sf in self.spectrumfiles():
            if len(sf.sample_names()) == 0:
                toremove.add(sf.id)
        if len(toremove) > 0:
            trans = self.begin()
            model.SpectrumFile.deleteMany(model.IN(model.SpectrumFile.q.id,toremove),connection=trans)
            model.normalize(connection=trans,quiet=self.quiet,alignedPeptidesOnly=False)
            trans.commit(close=True)

    def labels(self):
        return model.Label.select(connection=self.conn)

##     def deleteSpectrumFileGroups(self,type=None):
##         q = model.NOT(False)
##         if type:
##             q = (model.SpectrumFileGroup.q.type==type)
##         for sfgrp in model.SpectrumFileGroup.select(q,connection=self.conn):
##             sfgrp.destroySelf()

    def proteinids(self):
        return TableIDIter(model.Protein,self.conn)

    def proteins(self,nodecoys=False,connection=None):
        if not connection:
            connection = self.conn
        if nodecoys:
            return model.Protein.select(model.OR(model.Protein.q.decoy==None,model.Protein.q.decoy==False), connection=connection)
        return model.Protein.select(connection=connection)

    def protein(self,id=None,accession=None):
        assert id or accession, "At least one of protein id or accession must be provided!"
        if not id:
            id = self.proteinid(accession)
        if not id:
            return None
        return model.Protein.get(id,connection=self.conn)

    def proteinid(self,accession):
        try:
            return model.Protein.access(accession,connection=self.conn).id
        except AttributeError:
            pass
        return None

    def proteinids_byshortaccession(self,short_accession):
        ids = set()
        for pr in model.Protein.select(model.Protein.q.accession.startswith(short_accession),
                                       connection=self.conn):
            if pr.accession == short_accession or pr.short_accession() == short_accession:
                ids.add(pr.id)
        return ids

    def protein_count(self):
        return model.Protein.select(connection=self.conn).count()

    def proteinid_defline(self,id):
        md = model.Protein.get(id,connection=self.conn).metadata
        if md:
            return md.get('defline','')
        return ''

    def proteinid_length(self,id):
        md = model.Protein.get(id,connection=self.conn).metadata
        if md:
            return md.get('length',0)
        return 0

    def proteinid_accession(self,id):
        return model.Protein.get(id,connection=self.conn).accession

    def proteinid_description(self,id):
        md = model.Protein.get(id,connection=self.conn).metadata
        if md:
            return md.get('description','')
        return ''

    def peptideids(self):
        return TableIDIter(model.Peptide,self.conn)

    def proteinid_filter(self,keep=None,discard=None):
        trans = self.begin()
        for pr in model.Protein.select(connection=trans):
            if discard != None and pr.id in discard:
                pr.destroySelf()
                continue
            if keep != None and pr.id not in keep:
                pr.destroySelf()
                continue
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def geneid_filter(self,keep=None,discard=None):
        trans = self.begin()
        for pg in model.ProteinGroup.select(model.ProteinGroup.q.type=='Gene',connection=trans):
            if discard != None and pg.id in discard:
                pg.destroySelf()
                continue
            if keep != None and pg.id not in keep:
                pg.destroySelf()
                continue
        for pr in model.Protein.select(connection=trans):
            genes = 0
            for pg in pr.groups():
                if pg.type == 'Gene':
                    genes = 1
                    break
            if genes == 0:
                pr.destroySelf()
                continue
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def peptides(self,nodecoys=False,chunksize=1000,connection=None):
        if not connection:
            connection=self.conn
        if nodecoys:
            decoyconst = model.OR(model.Peptide.q.decoy==False,model.Peptide.q.decoy==None)
        else:
            decoyconst = (1 == 1)
        for i in range(self.peptide_minid(),self.peptide_maxid()+1,chunksize):
            for pep in model.Peptide.select(model.AND(decoyconst,
                                                      model.Peptide.q.id>=i,
                                                      model.Peptide.q.id<i+chunksize),
                                            connection=connection):
                yield pep
        return

    def peptide_maxid(self):
        try:
            return model.Peptide.select(orderBy='-id',connection=self.conn)[0].id
        except IndexError:
            return 0

    def peptide_minid(self):
        try:
            return model.Peptide.select(orderBy='id',connection=self.conn)[0].id
        except IndexError:
            return 0

    def distinct_peptides(self,proteins):
        return model.ProteinGroup.distinct_peptides(proteins)

    def unshared_peptideids(self,proteins=None):
        if proteins == None:
            if model.PeptideProteinCounts.empty(connection=self.conn):
                trans = self.begin()
                model.PeptideProteinCounts.compute(connection=trans)
                trans.commit(close=True)
            return model.PeptideProteinCounts.unshared_peptides(connection=self.conn)
        return [pep.id for pep in model.ProteinGroup.unshared_peptides(proteins,connection=self.conn)]

    def unshared_peptideids_bygene(self):
        if model.PeptideGeneCounts.empty(connection=self.conn):
            trans = self.begin()
            model.PeptideGeneCounts.compute(connection=trans)
            trans.commit(close=True)
        return model.PeptideGeneCounts.unshared_peptides(connection=self.conn)

    def peptide(self,id=None,sequence=None,connection=None):
        assert id or sequence, "At least one of peptide id or sequence must be provided!"
        if not connection:
            connection=self.conn
        if not id:
            id = self.peptideid(sequence,connection=connection)
            if not id:
                return None
            else:
                return id
        return model.Peptide.get(id,connection=connection)

    def peptide_count(self):
        return model.Peptide.select(connection=self.conn).count()

    def peptide_filter(self,condition,normalize=True):
        removed = 0
        trans = self.begin()
        for pep in self.peptides(connection=trans):
            if not condition(pep):
                pep.destroySelf()
                removed += 1
        if removed > 0 and normalize:
            model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def peptideid_sequence(self,id):
        return model.Peptide.get(id,connection=self.conn).sequence

    def alignmentids2pos(self):
        for al in model.Alignment.select(connection=self.conn):
            if al.start >= 0:
                yield (al.proteinID,al.peptideID),(al.start,al.end)

    def alignmentids2subst(self):
        for pr in self.proteins():
            for al in pr.alignments:
                if len(al.substitutionSites) > 0:
                    substsite = al.substitutionSites[0]
                    subst = substsite.substitution
                    yield (al.peptideID,pr.id),\
                          (subst.peptideAminoAcids,substsite.position,subst.proteinAminoAcids,subst.massDelta)

    def ions(self,chunksize=1000,connection=None):
        if not connection:
            connection=self.conn
        i = 0
        while True:
            any = False
            for ion in model.PeptideIon.select(model.AND(model.PeptideIon.q.id>=i,
                                                         model.PeptideIon.q.id<i+chunksize),
                                               connection=connection):
                any = True
                yield ion
            if not any:
                break
            i += chunksize
        return

    def ionids(self):
        return TableIDIter(model.PeptideIon)

    def ion(self,id):
        try:
            return model.PeptideIon.get(id,connection=self.conn)
        except model.SQLObjectNotFound:
            return None

    def ion_count(self):
        return model.PeptideIon.select(connection=self.conn).count()

    def ion_filter(self,condition):
        trans = self.begin()
        for ion in self.ions(connection=trans):
            if not condition(ion):
                ion.destroySelf()
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def modifications(self):
        return model.Modification.select(connection=self.conn)

    def modification(self,aa,delta):
        return model.Modification.access(aa,delta,connection=self.conn)

    def modString_mods(self,modString):
        retval = []
        for ms in modString.split(','):
            m = re.search(r'^([A-Z[])(\d+):(.*)$',ms)
            mod = model.Modification.find(m.group(1),float(m.group(3)),connection=self.conn)
            retval.append((int(m.group(2)),mod))
        return retval

    def psmids(self):
        return TableIDIter(model.PeptideSpectrumMatch)

    def psms_byvalue(self):
        for psm in  model.PeptideSpectrumMatch.select(orderBy='value',connection=self.conn):
            yield psm

    def psms(self,chunksize=1000,connection=None):
        if not connection:
            connection = self.conn
        for i in range(self.psm_minid(),self.psm_maxid()+1,chunksize):
            for psm in model.PeptideSpectrumMatch.select(model.AND(model.PeptideSpectrumMatch.q.id>=i,
                                                                   model.PeptideSpectrumMatch.q.id<i+chunksize),
                                                         connection=connection):
                yield psm
        return

    def psm_filter(self,minvalue=None,maxvalue=None,condition=None):
        trans = self.begin()
        if minvalue != None or maxvalue != None:
            model.PeptideSpectrumMatch.filter(minvalue,maxvalue,connection=trans)
        if condition:
            for psm in model.PeptideSpectrumMatch.select(connection=trans):
                if not condition(psm):
                    psm.destroySelf()
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def peptideid(self,sequence,connection=None):
        if not connection:
            connection = self.conn
        return model.Peptide.access(sequence,connection=connection)

    def psm_maxid(self):
        try:
            return model.PeptideSpectrumMatch.select(orderBy='-id',connection=self.conn)[0].id
        except IndexError:
            return 0

    def psm_minid(self):
        try:
            return model.PeptideSpectrumMatch.select(orderBy='id',connection=self.conn)[0].id
        except IndexError:
            return 0

    def psm_count(self):
        return model.PeptideSpectrumMatch.select(connection=self.conn).count()

    def psm(self,id):
        try:
            return model.PeptideSpectrumMatch.get(id,connection=self.conn)
        except model.SQLObjectNotFound:
            return None

    def spectrum_count(self):
        return model.Spectrum.select(connection=self.conn).count()

    def spectra_minid(self):
        try:
            return model.Spectrum.select(orderBy='id',connection=self.conn)[0].id
        except IndexError:
            return 0

    def spectra_maxid(self):
        try:
            return model.Spectrum.select(orderBy='-id',connection=self.conn)[0].id
        except IndexError:
            return 0

    def spectrum(self,id):
        try:
            return model.Spectrum.get(id,connection=self.conn)
        except model.SQLObjectNotFound:
            return None

    def spectra(self,chunksize=1000,connection=None):
        if not connection:
            connection=self.conn
        for i in range(self.spectra_minid(),self.spectra_maxid()+1,chunksize):
            for s in model.Spectrum.select(model.AND(model.Spectrum.q.id>=i,
                                                     model.Spectrum.q.id<i+chunksize),
                                          connection=self.conn):
                yield s
        return

    def deletePeptideGroups(self,type=None,connection=None):
        if not connection:
            connection = self.conn
        if type:
            model.PeptideGroup.deleteBy(type=type,connection=connection)
        else:
            model.PeptideGroup.deleteMany(None,connection=connection)
        model.Peptide2PeptideGroup.normalize(quiet=True,connection=connection)


    def peptideGroupsBySpectra(self):
        if self.peptidegroup_count('SpectraEquivClass') > 0:
            return
        trans = self.begin()
        model.PeptideGroup.insertSpectraGroups('SpectraEquivClass','SPE',connection=trans)
        trans.commit(close=True)

    def grouped_peptideids(self):
        self.peptideGroupsBySpectra()
        trans = self.begin()
        peptideids = set()
        for pg in model.PeptideGroup.select(model.AND(model.PeptideGroup.q.type=='SpectraEquivClass',
                                                      model.PeptideGroup.q.count > 1),
                                            connection=trans):
            for pep in pg.peptides():
                singleton = False
                for pg1 in pep.groups():
                    if pg1.count == 1:
                        singleton = True
                if not singleton:
                    peptideids.add(pep.id)
        trans.commit(close=True)
        return peptideids

    def removeMultiPeptideSpectra(self):
        trans = self.begin()

        badspecids = set()
        for pg in model.PeptideGroup.select(model.AND(model.PeptideGroup.q.type=='SpectraEquivClass',
                                                      model.PeptideGroup.q.count > 1),
                                      connection=trans):
            badspecids.update(s.id for s in pg.spectra())
        model.Spectrum.deleteMany(model.IN(model.Spectrum.q.id,badspecids),connection=trans)
        model.normalize(quiet=self.quiet,connection=trans)
        trans.commit(close=True)

    def removeAmbiguousMultiPeptideSpectra(self):
        trans = self.begin()

        badspecids = set()
        for pg in model.PeptideGroup.select(model.AND(model.PeptideGroup.q.type=='SpectraEquivClass',
                                                      model.PeptideGroup.q.count > 1),
                                      connection=trans):
            maxgrps = max([len(p.groups()) for p in pg.peptides()])
            if maxgrps > 1:
                badspecids.update(s.id for s in pg.spectra())
        model.Spectrum.deleteMany(model.IN(model.Spectrum.q.id,badspecids),connection=trans)
        model.normalize(quiet=self.quiet,connection=trans)
        trans.commit(close=True)

    def peptidegroup(self,id):
        return model.PeptideGroup.get(id,connection=self.conn)

    def peptidegroupids(self,type=None):
        if type != None:
            return TableIDIter(model.PeptideGroup,self.conn,
                               condition=(model.PeptideGroup.q.type == type))
        else:
            return TableIDIter(model.PeptideGroup,self.conn)

    def peptidegroup_count(self,type=None):
        q = None
        if type:
            q = (model.PeptideGroup.q.type==type)
        return model.PeptideGroup.select(q,connection=self.conn).count()

    def peptidegroups(self,type=None):
        q = None
        if type:
            q = (model.PeptideGroup.q.type==type)
        return model.PeptideGroup.select(q,connection=self.conn)

    def proteinid2peptidegroupids(self,type=None):
        return dict(model.PeptideGroup.protein2groups(type,connection=self.conn))

    def geneid2peptidegroupids(self,type=None):
        return dict(model.PeptideGroup.gene2groups(type,connection=self.conn))

    def geneid2peptideids(self):
        return model.Peptide.gene2peptides(connection=self.conn)

    def peptidegroupid2proteinids(self,type=None):
        return dict(model.PeptideGroup.group2proteins(type,connection=self.conn))

    def peptidegroupid2geneids(self,type=None):
        return dict(model.PeptideGroup.group2genes(type,connection=self.conn))

    def peptidegroupid2peptideids(self,type=None):
        return dict(model.PeptideGroup.group2peptides(type,connection=self.conn))

    def peptidegroupidproteinid2peptideids(self,type=None):
        return dict(model.PeptideGroup.groupprotein2peptides(type,connection=self.conn))

    def peptidegroupidgeneid2peptideids(self,type=None):
        return dict(model.PeptideGroup.groupgene2peptides(type,connection=self.conn))

    def peptidegroupid_peptides(self,id):
        return self.peptidegroup(id).peptides()

    def peptidegroupid_psmcounts(self,type=None):
        return dict(model.PeptideGroup.getSpecCounts('SPE',connection=self.conn))

    def peptidegroupid_psmminvalue(self,type=None):
        return dict(model.PeptideGroup.getMinFDR('SPE',connection=self.conn))

    def proteinid2peptideids(self):
        return RelatedJoinIDMapping(model.Alignment.q.proteinID,model.Alignment.q.peptideID,connection=self.conn)

##     def geneid2peptideids(self):
##         return RelatedJoinIDMapping1(model.ProteinGroup.q.id,
##                                      model.Alignment.q.peptideID,
##                                      model.AND(
##                                      model.Alignment.j.protein,
##                                      model.Protein2ProteinGroup.j.protein,
##                                         model.Protein2ProteinGroup.j.group,
##                                         model.ProteinGroup.q.type == 'Gene'),
##                                      connection=self.conn)

    def proteinids2peptideids(self,proteinids=None,peptideids=None):
        assert proteinids or peptideids;
        q = True
        if proteinids:
            q = model.AND(q,model.IN(model.Alignment.q.proteinID,proteinids))
        if peptideids:
            q = model.AND(q,model.IN(model.Alignment.q.peptideID,peptideids))
        als = defaultdict(set)
        for al in model.Alignment.select(q,connection=self.conn):
            als[al.proteinID].add(al.peptideID)
        return als

    def alignment_positions(self,proteinids=None,peptideids=None):
        assert proteinids or peptideids;
        q = model.AND(model.Alignment.q.start >= 0,
                      model.Alignment.q.end >= 0)
        if proteinids:
            q = model.AND(q,model.IN(model.Alignment.q.proteinID,proteinids))
        if peptideids:
            q = model.AND(q,model.IN(model.Alignment.q.peptideID,peptideids))
        als = {}
        for al in model.Alignment.select(q,connection=self.conn):
            als[(al.proteinID,al.peptideID)] = (al.start,al.end)
        return als

    def proteinid_coverage(self,proteinid,peptideids=None):
        plen = self.proteinid_length(proteinid)
        if not plen:
            raise RuntimeError("No protein length, so no coverage")
        alignments = self.alignment_positions(proteinids=[proteinid],peptideids=peptideids)
        if len(alignments) == 0:
            raise RuntimeError("No alignments, so no coverage")
        aacov = 0
        lastcovpos = -1
        for ps,pe in sorted(alignments.values()):
            if ps >= lastcovpos:
                aacov += (pe-ps)
                lastcovpos = pe
            elif pe > lastcovpos:
                aacov += (pe-lastcovpos)
                lastcovpos = pe
        return 100.0*(float(aacov)/plen)

    def peptideid2proteinids(self):
        return RelatedJoinIDMapping(model.Alignment.q.peptideID,model.Alignment.q.proteinID,connection=self.conn)

    def peptideid2geneids(self):
        return model.Peptide.peptide2genes(connection=self.conn)

    def protein_filter(self,keep=None,discard=None,condition=None,normalize=True):
        trans = self.begin()
        if not keep and not discard and not condition:
            return
        removed=0
        for pr in model.Protein.select(connection=trans):
            if keep and pr.accession not in keep:
                pr.destroySelf()
                removed += 1
            if discard and pr.accession in discard:
                pr.destroySelf()
                removed += 1
            if condition and condition(pr) == False:
                pr.destroySelf()
                removed += 1
        if removed > 0 and normalize:
            model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def protein_map(self,protmap):
        trans = self.begin()
        for prg in model.ProteinGroup.select(connection=trans):
            prg.destroySelf()
        for pr in model.Protein.select(connection=trans):
            if pr.accession not in protmap and pr.short_accession() not in protmap:
                pr.destroySelf()
                continue
            peps = pr.peptides()
            oldacc = pr.accession
            if pr.accession in protmap:
                oldacc1 = pr.accession
                newprs = protmap[oldacc1]
            else:
                oldacc1 = pr.short_accession()
                newprs = protmap[oldacc1]
            pr.destroySelf()
            for newpr in newprs:
                print(newpr)
                acc,gene = list(map(newpr.get,('accession','gene')))
                assert(acc)
                dbpr = model.Protein.find(acc,connection=trans)
                for pep in peps:
                    al = model.Alignment.find(pep,dbpr,connection=trans)
                prmd = {}
                for k,v in list(newpr.items()):
                    if k in ('accession',):
                        continue
                    if v in (oldacc,oldacc1):
                        continue
                    if v:
                        prmd[k] = v
                oldaccs = set([_f for _f in dbpr.getdata('mapped','').split(';') if _f])
                oldaccs.add(oldacc)
                prmd['mapped'] = ';'.join(sorted(oldaccs))
                dbpr.updatedata(prmd)
                if gene:
                    dbgene = model.ProteinGroup.find(gene,type='Gene',connection=trans)
                    if dbpr not in dbgene.proteins:
                        dbgene.addProtein(dbpr)
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def protein_orthologs(self,protmap):
        trans = self.begin()
        for pr in model.Protein.select(connection=trans):
            pr2 = None
            if pr.accession in protmap:
                pracc2 = protmap[pr.accession]
                pr2 = model.Protein.access(pracc2,connection=trans)
            elif pr.short_accession() in protmap:
                pracc2 = protmap[pr.short_accession()]
                pr2 = [pr1 for pr1 in model.Protein.select(model.AND(model.Protein.q.id != pr.id,
                                                            model.Protein.q.accession.startswith(pracc2)),
                                                       connection=trans) if pr1.short_accession() == pracc2 and \
                                         pr1.getdata('organism') != pr.getdata('organism')]
                if len(pr2) == 1:
                    pr2 = pr2[0]
                else:
                    pr2 = None
            if pr2 == None:
                continue

            if pr.peptides() != pr2.peptides():
                continue
            good = True
            for pep in pr.peptides():
                if pep.alignmentsto(pr)[0].metadata != pep.alignmentsto(pr2)[0].metadata:
                    good = False
                    break
            if not good:
                continue

            peps  = set(pr.peptides())
            peps &= set(pr2.peptides())

            acc = '|'.join([pr2.accession,pr.accession])
            dbpr = model.Protein.find(acc,connection=trans)
            if pr.decoy == pr2.decoy:
                dbpr.decoy = pr.decoy

            for pep in peps:
                dbal = model.Alignment.find(pep,dbpr,connection=trans)
                dbal.updatedata(pep.alignmentsto(pr)[0].metadata)
            for pg in pr.groups():
                model.Protein2ProteinGroup.insert(dbpr.id,pg.id,connection=trans)
            for pg in pr2.groups():
                model.Protein2ProteinGroup.insert(dbpr.id,pg.id,connection=trans)
            prmd = {'ortholog-pair': True,
                    'item1-acc': pr.accession,
                    'item2-acc': pr2.accession,
                    'item1-decoy': pr.decoy,
                    'item2-decoy': pr2.decoy,
                    'item1-metadata': pr.metadata,
                    'item2-metadata': pr2.metadata,
                    'item1-alignments':  [(al.peptideID,al.start,al.end,al.metadata) for al in pr.alignments],
                    'item2-alignments':  [(al.peptideID,al.start,al.end,al.metadata) for al in pr2.alignments],
                    }
            for k,v in list(pr2.metadata.items()):
                if pr2.getdata(k) != pr.getdata(k):
                    if ';' not in pr2.getdata(k) + pr.getdata(k):
                        prmd[k] = '; '.join([pr2.getdata(k),pr.getdata(k)])
                    else:
                        prmd[k] = ';;'.join([pr2.getdata(k),pr.getdata(k)])
                else:
                    prmd[k] = pr2.getdata(k)
            dbpr.updatedata(prmd)
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def undo_protein_orthologs(self):
        trans = self.begin()
        for pr in model.Protein.select(connection=trans):
            if not pr.getdata('ortholog-pair',False):
                continue
            dbpr1 = model.Protein.find(pr.getdata('item1-acc'),connection=trans)
            dbpr2 = model.Protein.find(pr.getdata('item2-acc'),connection=trans)
            dbpr1.decoy = pr.getdata('item1-decoy')
            dbpr2.decoy = pr.getdata('item2-decoy')
            dbpr1.updatedata(pr.getdata('item1-metadata'))
            dbpr2.updatedata(pr.getdata('item2-metadata'))
            for pepid,st,ed,md in pr.getdata('item1-alignments'):
                model.Alignment.insert(peptideid=pepid,proteinid=dbpr1.id,start=st,end=ed,metadata=md,connection=trans)
            for pepid,st,ed,md in pr.getdata('item2-alignments'):
                model.Alignment.insert(peptideid=pepid,proteinid=dbpr2.id,start=st,end=ed,metadata=md,connection=trans)
            pr.destroySelf()
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def gene_orthologs(self,protmap):
        trans = self.begin()
        for gn in model.ProteinGroup.select(model.ProteinGroup.q.type=='Gene',connection=trans):
            gn2 = None
            if gn.name in protmap:
                gnacc2 = protmap[gn.name]
                gn2 = model.ProteinGroup.access(gnacc2,connection=trans)
            if gn2 == None:
                continue
            assert(gn2.type == 'Gene')

            if set(gn.peptides()) != set(gn2.peptides()):
                continue

            peps  = set(gn.peptides())
            peps &= set(gn2.peptides())

            name = '|'.join([gn2.name,gn.name])
            dbgn = model.ProteinGroup.find(name,type='Gene',connection=trans)
            for pr in gn.proteins():
                model.Protein2ProteinGroup.insert(pr.id,dbgn.id,connection=trans)
            for pr in gn2.proteins():
                model.Protein2ProteinGroup.insert(pr.id,dbgn.id,connection=trans)
            gnmd = {'ortholog-pair': True,
                    'item1-name': gn.name,
                    'item2-name': gn2.name,
                    'item1-metadata': gn.metadata,
                    'item2-metadata': gn2.metadata,
                    'item1-proteinids': [pr.id for pr in gn.proteins()],
                    'item2-proteinids': [pr.id for pr in gn2.proteins()],
                    }
            for k,v in list(gn2.metadata.items()):
                if gn2.getdata(k) != gn.getdata(k):
                    if isinstance(gn2.getdata(k),list) or isinstance(gn.getdata(k),list):
                        gnmd[k] = gn2.getdata(k,[]) + [""] + gn.getdata(k,[])
                    else:
                        if ';' not in gn2.getdata(k,"") + gn.getdata(k,""):
                            gnmd[k] = '; '.join([gn2.getdata(k,""),gn.getdata(k,"")])
                        else:
                            gnmd[k] = ';;'.join([gn2.getdata(k,""),gn.getdata(k,"")])
                else:
                    gnmd[k] = gn2.getdata(k)
            dbgn.updatedata(gnmd)
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def undo_gene_orthologs(self):
        trans = self.begin()
        for gn in model.ProteinGroup.select(model.ProteinGroup.q.type=='Gene',connection=trans):
            if not gn.getdata('ortholog-pair',False):
                continue
            dbgn1 = model.ProteinGroup.find(gn.getdata('item1-name'),type='Gene',connection=trans)
            dbgn2 = model.ProteinGroup.find(gn.getdata('item2-name'),type='Gene',connection=trans)
            dbgn1.updatedata(gn.getdata('item1-metadata'))
            dbgn2.updatedata(gn.getdata('item2-metadata'))
            for prid in gn.getdata('item1-proteinids'):
                model.Protein2ProteinGroup.insert(prid,dbgn1.id,connection=trans)
            for prid in gn.getdata('item2-proteinids'):
                model.Protein2ProteinGroup.insert(prid,dbgn2.id,connection=trans)
            gn.destroySelf()
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def gene_map(self,genemap,remove=False):
        trans = self.begin()
        for prg in model.ProteinGroup.select(model.ProteinGroup.q.type=='Gene',connection=trans):
            prg.destroySelf()
        for pr in model.Protein.select(connection=trans):
            if remove and pr.accession not in genemap and pr.short_accession() not in genemap:
                pr.destroySelf()
                continue
            if pr.accession in genemap:
                genes = genemap[pr.accession]
            else:
                genes = genemap[pr.short_accession()]
            for gene in genes:
                dbgene = model.ProteinGroup.find(gene['accession'],type='Gene',connection=trans)
                gmd = {}
                for k,v in list(gene.items()):
                    if v in (None,"",dbgene.name,pr.accession,pr.short_accession()):
                        continue
                    gmd[k] = v
                dbgene.updatedata(gmd)
                dbgene.addProtein(pr)
        model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def peptide_remap(self,seqdb,pracc,aasubst=0,dnamut=0,validdeltas=None,validsubst=None,update=False,decoy=None,noprotein=False):
        trans = self.begin()
        if not update:
            model.Protein.clear(connection=trans)
            model.normalize(alignedPeptidesOnly=False,connection=trans,quiet=self.quiet)
        peptides = map(lambda pep: pep.sequence,model.Peptide.select(connection=trans))
        remapper = PeptideRemapper(peptides,seqdb,pracc,blocksize=10000,aasubst=aasubst,dnamut=dnamut)
        for pep in model.Peptide.select(connection=trans):
            for mapped in remapper.proteins(pep.sequence):
                acc,desc,laa,start,end,raa,defline,subst,deltamass,prlen = itemgetter(0,1,2,3,4,5,6,8,7,9)(mapped)
                # print pep.sequence,acc,desc,laa,start,end,raa,defline,subst,deltamass,prlen
                subst = [_f for _f in subst.split(',') if _f]
                if len(subst) > 1:
                    continue
                if len(subst) == 1:
                    good = False
                    m = re.search("^([A-Z]+)(\d+)\-\>([A-Z])+\((\d+)\)$",subst[0])
                    assert m
                    if validdeltas != None:
                        for d in validdeltas:
                            if abs(deltamass - d) <= 0.01:
                                good = True
                                break
                    if validsubst != None:
                        for f,t in validsubst:
                            if m.group(1) == f and m.group(3) == t:
                                good = True
                                break
                    if validsubst == None and validdeltas == None:
                        good = True
                    if not good:
                        continue
                shacc = None
                url = None
                gene = None
                entry = None
                org = None
                if defline:
                    shacc = pracc.shortacc(acc)
                    url = pracc.url(defline)
                    gene = pracc.gene(defline)
                    entry = pracc.entry(defline)
                    org = pracc.org(defline)
                prmd = {}
                prmd['description'] = desc
                prmd['defline'] = defline
                prmd['length'] = prlen
                if shacc:
                    prmd['short_accession'] = shacc
                if url:
                    prmd['url'] = url
                if org:
                    prmd['organism'] = org
                if gene:
                    prmd['gene'] = gene
                if entry:
                    prmd['entry'] = entry
                dbpr = model.Protein.access(acc,connection=trans)
                # print dbpr
                if not dbpr:
                    dbprid = model.Protein.insert(acc,decoy=decoy,metadata=prmd,connection=trans)
                else:
                    dbprid = dbpr.id
                    if prlen > dbpr.getdata('length',0):
                        dbpr.setdata('length',prlen)
                    if not dbpr.hasdata('description') or \
                       pracc.prefer((0,desc),(0,dbpr.getdata('description'))) < 0 or \
                       pracc.prefer((0,defline),(0,dbpr.getdata('defline'))) < 0:
                        dbpr.updatedata(prmd)
                # start = -1; end = -1;
                almd = {}
                if laa:
                    almd['left-amino-acid'] = laa[-1]
                    almd['left-context'] = laa
                if raa:
                    almd['right-amino-acid'] = raa[0]
                    almd['right-context'] = raa
                alid = model.Alignment.insert(pep.id,dbprid,start,end,metadata=almd,connection=trans)
                if decoy != None:
                    pep.decoy = decoy
                if len(subst) > 0:
                    # Can only be one...
                    m = re.search("^([A-Z]+)(\d+)\-\>([A-Z])+\((\d+)\)$",subst[0])
                    assert m
                    pos = int(m.group(2))
                    k = int(m.group(4))
                    sid = model.Substitution.insert(m.group(1),m.group(3),float(deltamass),connection=trans)
                    model.SubstitutionSite.insert(alid,pos,sid,connection=trans)
                if gene:
                    pgid = model.ProteinGroup.insert(gene,type='Gene',connection=trans)
                    model.Protein2ProteinGroup.insert(dbprid,pgid,connection=trans)
                if entry:
                    pgid = model.ProteinGroup.insert(entry,type='Entry',connection=trans)
                    model.Protein2ProteinGroup.insert(dbprid,pgid,connection=trans)
        model.normalize(connection=trans,alignedPeptidesOnly=(not noprotein),quiet=self.quiet)
        trans.commit(close=True)

    def deleteProteinGroups(self,type):
        for prgrp in model.ProteinGroup.selectBy(type=type,connection=self.conn):
            prgrp.destroySelf()

    def protein_group(self,id):
        return model.ProteinGroup.get(id,connection=self.conn)

    def getProteinGroup(self,name,type=None):
        prg = model.ProteinGroup.access(name,connection=self.conn)
        if prg and (not type or prg.type == type):
            return prg
        return None

    def newProteinGroup(self,name,proteinids,type=None):
        prgrp = model.ProteinGroup.find(name,type,connection=self.conn)
        prgrp.type = type
        for pr in prgrp.proteins():
            if pr.id not in proteinids:
                prgrp.removeProtein(pr,connection=self.conn)
        for prid in proteinids:
            prgrp.addProtein(model.Protein.get(prid,connection=self.conn))
        return prgrp

    def proteingroups(self,type=None,connection=None):
        if not connection:
            connection = self.conn
        q = None
        if type:
            q = (model.ProteinGroup.q.type==type)
        return model.ProteinGroup.select(q,connection=connection)

    def genes(self,nodecoys=False,connection=None):
        genes = self.proteingroups(type='Gene',connection=connection)
        if nodecoys:
            return itertools.iterfilter(lambda g: not g.getdata('decoy',False), genes)
        return genes

    def geneids(self):
        return list(map(attrgetter('id'),self.proteingroups(type='Gene')))

    def gene_filter(self,keep=None,discard=None,condition=None,normalize=True):
        self.proteingroup_filter(keep=keep,discard=discard,condition=condition,type='Gene',normalize=normalize)

    def proteingroup_count(self,type=None):
        q = None
        if type:
            q = (model.ProteinGroup.q.type==type)
        return model.ProteinGroup.select(q,connection=self.conn).count()

    def gene_count(self):
        return self.proteingroup_count(type='Gene')

    def geneid_name(self,id):
        return self.protein_group(id).name

    def geneid_description(self,id):
        return self.protein_group(id).getdata('description',"")

    def proteingroup_filter(self,keep=None,discard=None,condition=None,type=None,normalize=True):
        trans = self.begin()
        if not keep and not discard and not condition and not type:
            return
        removed=0
        for prg in self.proteingroups(type=type,connection=trans):
            if keep != None:
                if prg.name not in keep:
                    prg.destroySelf()
                    removed += 1
                continue
            if discard != None:
                if prg.name in discard:
                    prg.destroySelf()
                    removed += 1
                continue
            if condition != None:
                if not condition(prg):
                    prg.destroySelf()
                    removed += 1
                continue
            if type:
                if prg.type == type:
                    prg.destroySelf()
                    removed += 1
                continue
        if removed > 0 and normalize:
            model.normalize(connection=trans,quiet=self.quiet)
        trans.commit(close=True)

    def update(self,db2,connection=None):
        if not connection:
            connection=self.conn
        model.doInitTables(connection)
        for a2 in db2.analyses():
            a1 = model.Analysis.insert(a2.name,connection)
            if a1.hasdata('Filter'):
                assert a1.metadata == a2.metadata
            else:
                a1.updatedata(a2.metadata)
        for s2 in db2.samples():
            s1 = model.Sample.find(s2.name,type=s2.type,connection=connection)
        for l2 in db2.labels():
            l1 = model.Label.find(l2.name,type=l2.type,connection=connection)
        for sf2 in db2.spectrumfiles():
            sf1 = model.SpectrumFile.insert(sf2.name,connection=connection)
            for s2 in sf2.samples:
                l1 = model.Label.selectBy(name=s2.label.name,connection=connection)[0]
                s1 = model.Sample.selectBy(name=s2.sample.name,connection=connection)[0]
                model.SampleSpectrumFileLabel.insert(s1,sf1,l1,connection=connection)
        for p2 in db2.peptides():
            p = model.Peptide.insert(p2.sequence,p2.decoy,connection=connection)
        for m2 in db2.modifications():
            m1 = model.Modification.insert(m2.aminoAcid,m2.massDelta,connection=connection)
        for i,psm2 in enumerate(db2.psms()):
            s2 = psm2.spectrum
            fid = model.SpectrumFile.find(s2.file.name,connection=connection).id
            sid = model.Spectrum.insert(fid,s2.start_scan,s2.end_scan,precursorMz=s2.precursorMz,metadata=s2.metadata,connection=connection)
            # print i,psm2,s2,sid,model.Spectrum.get(sid,connection=connection)
            aid = model.Analysis.insert(psm2.analysis.name,connection=connection).id
            pi2 = psm2.peptideIon
            p1 = model.Peptide.access(pi2.peptide.sequence,connection=connection)
            mods = []
            for ms2 in pi2.modifiedSites:
                mods.append((ms2.position,self.modification(ms2.modification.aminoAcid,ms2.modification.massDelta)))
            piid = model.PeptideIon.insert(p1.id,mods,pi2.charge,mw=pi2.mw,mz=pi2.mz,connection=connection)
            model.PeptideSpectrumMatch.insert(sid,piid,aid,psm2.value,metadata=psm2.metadata,connection=connection)
        for pg2 in db2.proteingroups():
            pg = model.ProteinGroup.insert(pg2.name,type=pg2.type,connection=connection)
        for p2 in db2.proteins():
            pid = model.Protein.insert(p2.accession,p2.decoy,p2.metadata,connection=connection)
            for al2 in p2.alignments:
                pepid = model.Peptide.access(al2.peptide.sequence,connection=connection).id
                al = model.Alignment.insert(pepid,pid,al2.start,al2.end,metadata=al2.metadata,connection=connection)
            for pg2 in p2.groups():
                gid = model.ProteinGroup.insert(pg2.name,type=pg2.type,connection=connection)
                p2pgid = model.Protein2ProteinGroup.insert(pid,gid,connection=connection)

        model.normalize(quiet=False,connection=connection)

    def dump(self,wh=sys.stderr):
        model.dump_counts(wh,connection=self.conn)

    def peptideid_psmcounts0(self):
        if model.PeptidePSMCounts.empty(connection=self.conn):
            trans = self.begin()
            model.PeptidePSMCounts.compute(connection=trans)
            trans.commit(close=True)
        return model.PeptidePSMCounts.counts(connection=self.conn)

    def peptideid_psmcounts(self):
        return model.PeptidePSMCounts.counts1(connection=self.conn)

    def peptideid_psmminvalue0(self):
        if model.PeptidePSMMinValue.empty(connection=self.conn):
            trans = self.begin()
            model.PeptidePSMMinValue.compute(connection=trans)
            trans.commit(close=True)
        return model.PeptidePSMMinValue.minvalues(connection=self.conn)

    def peptideid_psmminvalue(self):
        return model.PeptidePSMMinValue.minvalues1(connection=self.conn)

    def peptideid_psmmindelta(self):
        return model.PeptidePSMMinDelta.minvalues1(connection=self.conn)

    def fileid_peptideid_rows(self):
        if not model.FilePeptidePSMCounts.empty(connection=self.conn):
            return model.FilePeptidePSMCounts.select(connection=self.conn).count()
        elif not model.FilePeptidePSMMinValue.empty(connection=self.conn):
            return model.FilePeptidePSMMinValue.select(connection=self.conn).count()
        self.fileid_peptideid_psmcounts()
        return model.FilePeptidePSMCounts.select(connection=self.conn).count()

    def fileid_peptideid_psmcounts(self):
        if model.FilePeptidePSMCounts.empty(connection=self.conn):
            trans = self.begin()
            model.FilePeptidePSMCounts.compute(connection=trans)
            trans.commit(close=True)
        return FileValueMapping(model.FilePeptidePSMCounts.q.peptide,
                                model.FilePeptidePSMCounts.q.count,
                                connection=self.conn)

    def fileid_peptideid_psmminvalue(self):
        if model.FilePeptidePSMMinValue.empty(connection=self.conn):
            trans = self.begin()
            model.FilePeptidePSMMinValue.compute(connection=trans)
            trans.commit(close=True)
        return FileValueMapping(model.FilePeptidePSMMinValue.q.peptide,
                                model.FilePeptidePSMMinValue.q.minvalue,
                                connection=self.conn)

    def fileid_peptideid_ioncounts(self):
        if model.FilePeptideIonCounts.empty(connection=self.conn):
            trans = self.begin()
            model.FilePeptideIonCounts.compute(connection=trans)
            trans.commit(close=True)
        return FileValueMapping(model.FilePeptideIonCounts.q.peptide,
                                model.FilePeptideIonCounts.q.count,
                                connection=self.conn)

    def fileid_peptideid_specfunc(self,fn):
        if model.FilePeptideSpectra.empty(connection=self.conn):
            trans = self.begin()
            model.FilePeptideSpectra.compute(connection=trans)
            trans.commit(close=True)
        return FileSpectraFunctionMapping(model.FilePeptideSpectra.q.peptide,fn,
                                          connection=self.conn)

    def fileid_peptideionid_specfunc(self,fn):
        if model.FilePeptideIonSpectra.empty(connection=self.conn):
            trans = self.begin()
            model.FilePeptideIonSpectra.compute(connection=trans)
            trans.commit(close=True)
        return FileSpectraFunctionMapping(model.FilePeptideIonSpectra.q.peptideion,fn,
                                          connection=self.conn)

    def fileid_proteinid_psmcounts(self):
        if model.FileProteinPSMCounts.empty(connection=self.conn):
            trans = self.begin()
            model.FileProteinPSMCounts.compute(connection=trans)
            trans.commit(close=True)
        return FileValueMapping(model.FileProteinPSMCounts.q.protein,
                                model.FileProteinPSMCounts.q.count,
                                connection=self.conn)

from sqlobject.sqlbuilder import Select, func, SQLTrueClause
import collections


class TableIDIter(collections.Set):
    def __init__(self,cls,connection,condition=None):
        self._cls = cls
        self._conn = connection
        self._condition = condition if condition else SQLTrueClause
    def __contains__(self,x):
        for r in self._conn.queryAll(self._conn.sqlrepr(Select(self._cls.q.id,where=model.AND(self._condition,self._cls.q.id==x),limit=1))):
            return True
        return False
    def __iter__(self):
        return map(itemgetter(0),self._conn.queryAll(self._conn.sqlrepr(Select(self._cls.q.id,where=self._condition))))
    def __len__(self):
        return self._cls.select(self._condition,connection=self._conn).count()

class RelatedJoinIDMapping(collections.Mapping):
    def __init__(self,frm,to,connection):
        self._cls = frm.soClass
        self._conn = connection
        self._from = frm
        self._to = to
    def __contains__(self,x):
        for r in self._conn.queryAll(self._conn.sqlrepr(Select(self._from,where=(self._from==x),limit=1))):
            return True
        return False
    def __iter__(self):
        return map(itemgetter(0),self._conn.queryAll(self._conn.sqlrepr(Select(func.DISTINCT(self._from)))))
    def __len__(self):
        return self._conn.queryOne(self._conn.sqlrepr(Select(func.COUNT(func.DISTINCT(self._from)),limit=1)))[0]
    def __getitem__(self,x):
        return set(map(itemgetter(0),self._conn.queryAll(self._conn.sqlrepr(Select(self._to,where=(self._from==x))))))

class RelatedJoinIDMapping1(collections.Mapping):
    def __init__(self,frm,to,clause,connection):
        self._conn = connection
        self._from = frm
        self._to = to
        self._clause = clause

    def __contains__(self,x):
        for r in self._conn.queryAll(self._conn.sqlrepr(Select(self._from,where=(self._from==x),limit=1))):
            return True
        return False
    def __iter__(self):
        return map(itemgetter(0),self._conn.queryAll(self._conn.sqlrepr(Select(func.DISTINCT(self._from)))))
    def __len__(self):
        return self._conn.queryOne(self._conn.sqlrepr(Select(func.COUNT(func.DISTINCT(self._from)),limit=1)))[0]
    def __getitem__(self,x):
        return set(map(itemgetter(0),self._conn.queryAll(self._conn.sqlrepr(Select(self._to,where=model.AND(self._from==x,self._clause))))))

class FileValueMapping(collections.Mapping):
    def __init__(self,by,value,connection):
        self._cls = by.soClass
        self._col = by
        self._file  = self._cls.q.file
        self._value = value
        self._conn = connection

    def __contains__(self,x):
        for r in self._conn.queryAll(self._conn.sqlrepr(Select(self._cls.q.id,where=model.AND(self._file==x[0],
                                                                                              self._col==x[1]),limit=1))):
            return True
        return False
    def iteritems(self):
        return map(lambda t: ((t[0],t[1]),t[2]), self._conn.queryAll(self._conn.sqlrepr(Select((self._file,self._col,self._value)))))
    def __iter__(self):
        return iter(self._conn.queryAll(self._conn.sqlrepr(Select((self._file,self._col)))))
    def __len__(self):
        return self._cls.select(connection=self._conn).count()
    def __getitem__(self,x):
        for r in self._conn.queryAll(self._conn.sqlrepr(Select(self._value,where=model.AND(self._file==x[0],self._col==x[1])))):
            return r[0]
        return 0

class FileSpectraFunctionMapping(collections.Mapping):
    def __init__(self,by,fn,connection):
        self._cls = by.soClass
        self._col = by
        self._file  = self._cls.q.file
        self._spectra = self._cls.q.spectra
        self._fn = fn
        self._conn = connection

    def __contains__(self,x):
        for r in self._conn.queryAll(self._conn.sqlrepr(Select(self._cls.q.id,where=model.AND(self._file==x[0],
                                                                                              self._col==x[1]),limit=1))):
            return True
        return False
    def iteritems(self):
        return map(lambda t: ((t[0],t[1]),self._fn(t[1],[model.Spectrum.get(int(id),connection=self._conn) for id in t[2].split(',')])), self._conn.queryAll(self._conn.sqlrepr(Select((self._file,self._col,self._spectra)))))
    def __iter__(self):
        return iter(self._conn.queryAll(self._conn.sqlrepr(Select((self._file,self._col)))))
    def __len__(self):
        return self._cls.select(connection=self._conn).count()
    def __getitem__(self,x):
        for r in self._conn.queryAll(self._conn.sqlrepr(Select(self._spectra,where=model.AND(self._file==x[0],self._col==x[1])))):
            return self._fn(x[1],[model.Spectrum.get(int(id),connection=self._conn) for id in r[0].split(',')])
        return None
