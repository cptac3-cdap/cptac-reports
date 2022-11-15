#!/bin/env python

from .model import *
import gzip

from .PepXML import PepXMLReader
from .MZIdentML import MZIdentMLReader
from .idpdb import IDPDB
from peptidescan.PeptideRemapper import FirstWord, TestAccessionRules
from fragger.MassTable import MonoisotopicMassTable
from .TandemResults import TandemResultsReader
from .DIATSV import DIATSVReader

mt = MonoisotopicMassTable()

def parserSwitch(f,alignments=True):
    for extn in ('pepxml','pep.xml'):
        if f.lower().endswith(extn):
            return PepXMLReader(open(f,'r'))
        elif f.lower().endswith(extn+'.gz'):
            return PepXMLReader(gzip.open(f,'r'))
    for extn in ('mzid','mzIdentML'):
        if f.lower().endswith(extn):
            return MZIdentMLReader(open(f,'r'),alignments=alignments,filename=f)
        elif f.lower().endswith(extn+'.gz'):
            return MZIdentMLReader(gzip.open(f,'r'),alignments=alignments,filename=f)
    for extn in ('idpdb',):
        if f.lower().endswith(extn):
            return IDPDB(f)
    if '.elib.' in f:
        return DIATSVReader(f,alignments=alignments)
    if '.msgfplus.' in f:
        if f.lower().endswith('.gz'):
            return MZIdentMLReader(gzip.open(f,'r'),alignments=alignments,filename=f)
        else:
            return MZIdentMLReader(open(f,'r'),alignments=alignments,filename=f)
    try:
        return TandemResultsReader(f)
    except RuntimeError:
        pass
    return None

def SearchResultLoad(f,psmfilter,noprotein=False,accrule=None,reload=False,analysisname=None,analysisroot=None,connection=None):
    autoacc = True
    if accrule != None:
        autoacc = False

    if analysisname:
        aname = analysisname
        reload = True
    elif analysisroot:
        assert os.path.isdir(analysisroot)
        assert f.startswith(analysisroot+os.sep)
        aname = f[len(analysisroot)+1:]
    else:
        aname = os.path.split(f)[1]

    if not reload and Analysis.selectBy(name=aname,connection=connection).count() > 0:
        print("Analysis %s already loaded"%(aname,), file=sys.stderr)
        return False
    aid = Analysis.insert(aname,connection=connection).id
    a = Analysis.get(aid,connection=connection)
    a.updatedata({'Filter':psmfilter.__class__.__name__, 'Threshold':psmfilter.loose})

    for res in parserSwitch(f,alignments=(not noprotein)):
        for psm in res.getPSMs():
            if psmfilter and not psmfilter.retain(psm):
                continue
            specfileid = SpectrumFile.insert(res.get('basename'),connection=connection).id
            start_scan = int(str(res.get('start_scan')))
            end_scan = int(str(res.get('end_scan')))
            smd = {}
            z = psm.get('assumed_charge')
            pmr = psm.get('precursor_neutral_mass')
            pmz = None
            if z != None and pmr != None:
                z = int(z)
                pmr = float(pmr)
                pmz = (pmr+1.0078*z)/z
            rt = res.get('retention_time')
            if rt != None:
                smd['retention_time'] = rt
            precint = res.get('precursor_intensity')
            if precint != None:
                smd['precursor_intensity'] = precint
            precarea = res.get("precursor_area")
            if precarea != None:
                smd['precursor_area'] = precarea
            for k in (113,114,115,116,117,118,119,121):
                for ty in (4,8,""):
                    k1 = "iTRAQ%s-%s_reporter_ion_intensity"%(ty,k)
                    k2 = "iTRAQ%s-%s_reporter_ion_massdelta"%(ty,k)
                    if res.has(k1):
                        smd[k1] = res.get(k1)
                    if res.has(k2):
                        smd[k2] = res.get(k2)
            for k in (126,127,128,129,130,131,"126N","126C","127N","127C","128N","128C","129N","129C","130N","130C","131N","131C","132N","132C","133N","133C","134N","134C","135N",):
                for ty in (2,6,10,11,16,18,"XX",""):
                    k1 = "TMT%s-%s_reporter_ion_intensity"%(ty,k)
                    k2 = "TMT%s-%s_reporter_ion_massdelta"%(ty,k)
                    if res.has(k1):
                        smd[k1] = res.get(k1)
                    if res.has(k2):
                        smd[k2] = res.get(k2)
            if len(smd) == 0:
                smd = None
            spectrumid = Spectrum.insert(specfileid,start_scan,end_scan,precursorMz=pmz,metadata=smd,connection=connection)
            peptideid=Peptide.insert(psm.get('peptide'),psm.get('decoy'),connection=connection)
            cmr = mt.peptide_mr(psm.get('peptide'))
            mods = []
            if psm.get('mods') not in ("-",None):
                for m in psm.get('mods').split(','):
                    aapos,delta = m.split(':')
                    aa = aapos[0]
                    pos = int(aapos[1:])
                    delta = float(delta)
                    cmr += delta
                    mods.append((pos,Modification.insert(aa,delta,connection=connection)))
            mw = cmr
            mz = None
            if z:
                mz = (cmr+1.0078*z)/z
            peptideionid=PeptideIon.insert(peptideid,mods,z,mw=mw,mz=mz,connection=connection)
            if psmfilter:
                value = psmfilter.psmqvalue(psm)
            psmmd={}
            if psmfilter:
                for k,v in list(psmfilter.psmcriteria(psm).items()):
                    psmmd[k] = v
            if pmr:
                psmmd['massdelta'] = pmr-cmr
            if psm.get('abundance'):
                psmmd['abundance'] = psm.get('abundance')
            if psm.get('numFragments'):
                psmmd['numFragments'] = psm.get('numFragments')
            for k in psm.get('metrics',[]):
                if k not in psmmd:
                    psmmd[k] = psm.get(k)
            psmdbid = PeptideSpectrumMatch.insert(spectrumid,peptideionid,aid,value,metadata=psmmd,connection=connection)
            if noprotein:
                continue
            if not psm.get('protein'):
                print(psm)
            praccs = psm.get('protein').split(',')
            if psm.has('protein_descr'):
                prdescs = psm.get('protein_descr').split('\t')
                prdescs += [None]*len(praccs)
            else:
                prdescs = [None]*len(praccs)
            if psm.has('protein_org'):
                prorgs = psm.get('protein_org').split(',')
                prorgs += [None]*len(praccs)
            else:
                prorgs = [None]*len(praccs)
            if psm.has('protein_source'):
                prsrcs = psm.get('protein_source').split(',')
                prsrcs += [None]*len(praccs)
            else:
                prsrcs = [None]*len(praccs)
            if psm.has('peptide_offset'):
                offsets = psm.get('peptide_offset').split(',')
                offsets += [None]*len(praccs)
            else:
                offsets = [None]*len(praccs)
            if psm.has('num_tol_term'):
                ntts = str(psm.get('num_tol_term')).split(',')
                ntts += [None]*len(praccs)
            else:
                ntts = [None]*len(praccs)
            if psm.has('peptide_prev_aa'):
                laas = psm.get('peptide_prev_aa').split(',')
                laas += [None]*len(praccs)
            else:
                laas = [None]*len(praccs)
            if psm.has('peptide_next_aa'):
                raas = psm.get('peptide_next_aa').split(',')
                raas += [None]*len(praccs)
            else:
                raas = [None]*len(praccs)
            for pracc,prdesc,prorg,prsrc,offset,ntt,laa,raa in zip(praccs,prdescs,prorgs,prsrcs,offsets,ntts,laas,raas):
                prmd = None
                if prdesc:
                    if not accrule and autoacc:
                        accrule = TestAccessionRules(prdesc)
                        if isinstance(accrule,FirstWord):
                            accrule = None
                        autoacc = False
                    if accrule:
                        prdefline = prdesc
                        ruleacc,prdesc = accrule(prdefline)
                        prmd = {}
                        if prdesc != None:
                            prmd['defline'] = prdefline
                            prmd['description'] = prdesc
                            pracc = ruleacc
                            shacc = None
                            url = None
                            gene = None
                            entry = None
                            org = None
                            if prdefline:
                                shacc = accrule.shortacc(ruleacc)
                                url = accrule.url(prdefline)
                                gene = accrule.gene(prdefline)
                                entry = accrule.entry(prdefline)
                                org = accrule.org(prdefline)
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
                        else:
                            prmd['description'] = prdefline
                    else:
                        prmd = {}
                        prmd['description'] = prdesc
                    if prorg and 'organism' not in prmd:
                        prmd['organism'] = prorg
                    if prsrc and 'source' not in prmd:
                        prmd['source'] = prsrc
                    # print prdefline
                    # print pracc,prmd
                else:
                    m = re.search(r'^(.*)[.-]\d$',pracc)
                    if m:
                        prmd = {}
                        prmd['short_accession'] = m.group(1)

                prdbid = Protein.insert(pracc,decoy=psm.get('decoy'),metadata=prmd,connection=connection)
                start=-1; end=-1
                if offset:
                    start = int(offset); end = int(offset)+len(psm.get('peptide'))
                alignmd = {}
                if ntt != None:
                    alignmd['specific-termini'] = ntt
                if laa != None:
                    alignmd['left-amino-acid'] = laa
                if raa != None:
                    alignmd['right-amino-acid'] = raa
                alignid = Alignment.insert(peptideid,prdbid,start,end,metadata=alignmd,connection=connection)

    return True

if __name__ == '__main__':

    import sys

    conn,todele=init(filename='psmdb.db3',clear=True)

    from .filter import SpectralFDR
    psmfilter = SpectralFDR(0.1)

    for f in sys.argv[1:]:
        SearchResultLoad(f, psmfilter, connection=conn)

    dump_counts()
