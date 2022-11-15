#!/usr/bin/env python27

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
 
from dataset import CSVFileTable, XLSFileTable, XLSFileDataset, ZIPFileDataset, ParsimonyCSV
import sys, os, os.path, copy, re, math, gzip, bz2
from math import log
from collections import defaultdict
from operator import itemgetter, attrgetter
from functools import partial
from psmdb.PepXML import PepXMLReader, ParsimonyPepXML
from peptidescan.PeptideRemapper import *
from peptidescan.OutOfCoreTable import *
from itertools import combinations
import numpy as np
import scipy.optimize

import atexit
def onexit():
    print("Type <Enter> to Exit...", end=' ', file=sys.stderr)
    sys.stderr.flush()
    input()

from optparse_gui import OptionParser, OptionGroup, GUI, UserCancelledError, ProgressText
VERSION = '0.9999'

if GUI() and len(sys.argv) == 1:
    from optparse_gui import OptionParserGUI
    parser = OptionParserGUI(version=VERSION)
    atexit.register(onexit)
    error_kwargs = {'exit': False}
else:
    parser = OptionParser(version=VERSION)
    error_kwargs = {}

from psmdb.PSMDb import PSMDb

parser.add_option("-d","--database",type="file",dest="database",remember=True,
                  help="PSM database file",default="",notNone=True,
                  name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
# parser.add_option("-c","--common",type="string",dest="common",remember=True,
#                   help="Ratios with respect to common sample or reporter-ion.",
#                   default="iTRAQ-117",notNone=True,
#                   name="Common Sample/Channel")
parser.add_option("--labels",type="choice",dest="labels",choices=["TMT10","TMT10+2","TMT11","TMT11+2","TMT16","TMT16+2","TMT18","TMT6","iTRAQ","iTRAQ4","iTRAQ8"],default="TMT10",
                  help="Labeling reagent type. One of TMT6, TMT10, TMT10+2, TMT11, TMT11+2, TMT16, TMT16+2, TMT18, iTRAQ, iTRAQ4, iTRAQ8. Default: TMT10.")
parser.add_option("--nocorrections",action="store_true",dest="nocorrections",default=False,
                  help="Do not apply isotope corrections even if available. Default: False.")
parser.add_option("--correction",type="choice",dest="correction",default="NNLS+Check",
                  help="Strategy for computing corrected isotopic-label intensity values. One of None, LS, LS+Check, NNLS, NNLS+Check. Default: NNLS+Check.",
                  choices=["None","LS","LS+Check","NNLS","NNLS+Check"])
parser.add_option("--count",type="int",dest="count",remember=True,
                  help="Observations in at least this many samples.", default=0,
                  name="Observations")
parser.add_option("--ratio",type="float",dest="ratio",remember=True,
                  help="Observations in at least this proportion of the samples.", default=0.0,
                  name="Ratio")
parser.add_option("--report",type="choice",dest="report",default='Peptide',
                  help="Report to produce. One of Peptide, ModifiedPeptide, PeptideForm, PeptideIon, Phosphopeptide, PhosphopeptideSite, Phosphosite, PhosphositeCombination, Glycopeptide, GlycositeCombination, Acetylpeptide, AcetylsiteCombination, Ubiquitylpeptide, UbiquitylsiteCombination. Default: Peptide.",
                  choices=['Peptide','ModifiedPeptide','PeptideIon', 'PeptideForm','Phosphopeptide','PhosphopeptideSite','Phosphosite','PhosphositeCombination', 'Glycopeptide', 'GlycositeCombination','Acetylpeptide','AcetylsiteCombination','Ubiquitylpeptide', 'UbiquitylsiteCombination'],
                  remember=True,name="Report")
parser.add_option("--specific",action="store_true",dest="specific",remember=True,
                  help="Specific (tryptic) peptides only. Default: False.", name="Specfific?", default=False)
parser.add_option("--peptide",type="string",dest="peptide",remember=True,
                  help="Focus on peptide", default="",
                  name="Peptide")
parser.add_option("--locus",type="string",dest="locus",remember=True,
                  help="Focus on locus", default="",
                  name="Locus")
parser.add_option("--protein",type="string",dest="protein",remember=True,
                  help="Focus on protein accession", default="",
                  name="Protein")
parser.add_option("--peptideregex",type="string",dest="pepseqregex",remember=True,
                  help="Focus on peptides with matching sequences", default="",
                  name="Peptide Sequence Regular Expression")
parser.add_option("--peptidentermaas",type="string",dest="pepseqntermaas",remember=True,
                  help="Focus on peptides with specific initial AAs", default="",
                  name="Peptide Sequence N-term AAs")
## parser.add_option("--modifiedpeptide",action="store_true",
##                   dest="modifiedpeptide",default=False,remember=True,
##                   name="Modified Peptide",help="Do not collapse modifications. Default: False.")
## parser.add_option("--peptide",action="store_true",
##                   dest="peptide",default=False,remember=True,
##                   name="Peptide",help="Collapse all modifications. Default: False.")
## parser.add_option("--phosphopeptide",action="store_true",
##                   dest="phosphopeptide",default=False,remember=True,
##                   name="Phosphopeptide",help="Collapse non-phosphopeptide modifications. Default: False.")
## parser.add_option("--phosphopeptidesite",action="store_true",
##                   dest="phosphopeptidesite",default=False,remember=True,
##                   name="Phosphopeptide Site",help="Collapse non-phosphopeptide modifications and aggregate at each phosphorylated site. Default: False.")
parser.add_option("--localized",action="store_true",
                  dest="fullylocalized",default=False,remember=True,
                  name="Localized",help="Fully localized PSMs only. Default: False.")
parser.add_option("--noties",action="store_true",
                  dest="noties",default=False,remember=True,
                  name="No PSM ties (often modification rearrangements",help="No PSM ties. Default: False.")
parser.add_option("--extra",action="store_true",
                  dest="extraheaders",default=False,remember=True,
                  name="Extra headers",help="Extra headers. Default: False.")
parser.add_option("--allratios",action="store_true",
                  dest="allratios",default=False,remember=True,
                  name="All Ratios", help="Output all ratio columns, even if not observed. Default: False.")
parser.add_option("--append",action="store_true",
                  dest="append",default=False,remember=True,
                  name="Append to output file",help="Append to output file. Default: False.")
parser.add_option("-q","--quiet",action="store_true",
                  dest="quiet",default=False,remember=True,
                  name="Quiet",help="Quiet. Default: False.")
parser.add_option("-o","--output",type="savefile",
                  dest="output",default="",remember=True,
                  help="Output file.",
                  name="Output",
                  filetypes=[("Counts Table","*.csv;*.tsv;*.xlsx;*.xls"),
                             ("CSV File","*.csv"),
                             ("TSV File","*.tsv"),
                             ("Excel File","*.xlsx;*.xls")])

opts = None
while True:
    if 'exit' in error_kwargs:
        try:
            opts,args = parser.parse_args(opts=opts)
        except UserCancelledError:
            os._exit(0);
    else:
        opts,args = parser.parse_args()

    break

if opts.correction == "None":
    opts.nocorrections = True

itraq4labels=[_f for _f in """
    iTRAQ-114
    iTRAQ-115
    iTRAQ-116
    iTRAQ-117
""".split() if _f]

itraq8labels=[_f for _f in """
    iTRAQ-113
    iTRAQ-114
    iTRAQ-115
    iTRAQ-116
    iTRAQ-117
    iTRAQ-118
    iTRAQ-119
    iTRAQ-121
""".split() if _f]

tmt10labels=[_f for _f in """
    TMT10-126
    TMT10-127N
    TMT10-127C
    TMT10-128N
    TMT10-128C
    TMT10-129N
    TMT10-129C
    TMT10-130N
    TMT10-130C
    TMT10-131
""".split() if _f]

tmt10extralabels=[_f for _f in """
    TMTXX-131C
    TMTXX-132N
""".split() if _f]

tmt6labels=[_f for _f in """
    TMT6-126
    TMT6-127
    TMT6-128
    TMT6-129
    TMT6-130
    TMT6-131
""".split() if _f]

tmt11labels=[_f for _f in """
    TMT11-126C
    TMT11-127N
    TMT11-127C
    TMT11-128N
    TMT11-128C
    TMT11-129N
    TMT11-129C
    TMT11-130N
    TMT11-130C
    TMT11-131N
    TMT11-131C
""".split() if _f]

tmt11extralabels=[_f for _f in """
    TMTXX-132N
    TMTXX-132C
""".split() if _f]

tmt16labels=[_f for _f in """
    TMT16-126C
    TMT16-127N
    TMT16-127C
    TMT16-128N
    TMT16-128C
    TMT16-129N
    TMT16-129C
    TMT16-130N
    TMT16-130C
    TMT16-131N
    TMT16-131C
    TMT16-132N
    TMT16-132C
    TMT16-133N
    TMT16-133C
    TMT16-134N
""".split() if _f]

tmt16extralabels=[_f for _f in """
    TMTXX-134C
    TMTXX-135N
""".split() if _f]

tmt18labels=[_f for _f in """
    TMT18-126C
    TMT18-127N
    TMT18-127C
    TMT18-128N
    TMT18-128C
    TMT18-129N
    TMT18-129C
    TMT18-130N
    TMT18-130C
    TMT18-131N
    TMT18-131C
    TMT18-132N
    TMT18-132C
    TMT18-133N
    TMT18-133C
    TMT18-134N
    TMT18-134C
    TMT18-135N
""".split() if _f]

extralabels = []
if opts.labels == 'TMT10':
    labelnames = tmt10labels
elif opts.labels == 'TMT10+2':
    labelnames = tmt10labels
    extralabels = tmt10extralabels
    assert(not opts.nocorrections)
    assert(opts.correction in ("NNLS","NNLS+Check"))
elif opts.labels == 'TMT6':
    labelnames = tmt6labels
elif opts.labels == 'TMT11':
    labelnames = tmt11labels
elif opts.labels == 'TMT11+2':
    labelnames = tmt11labels
    extralabels = tmt11extralabels
    assert(not opts.nocorrections)
    assert(opts.correction in ("NNLS","NNLS+Check"))
elif opts.labels == 'TMT16':
    labelnames = tmt16labels
elif opts.labels == 'TMT16+2':
    labelnames = tmt16labels
    extralabels = tmt16extralabels
    assert(not opts.nocorrections)
    assert(opts.correction in ("NNLS","NNLS+Check"))
elif opts.labels == 'TMT18':
    labelnames = tmt18labels
elif opts.labels in ('iTRAQ','iTRAQ4'):
    labelnames = itraq4labels
elif opts.labels == 'iTRAQ8':
    labelnames = itraq8labels
else:
    raise RuntimeError("Bad label type %s"%opts.labels)

progress = ProgressText(quiet=(opts.quiet or (not opts.output)))

psmdb = PSMDb(filename=opts.database,quiet=opts.quiet)

progress.stage("Extracting experimental design information")

fileidlabel2sampid = dict()
fileidlabel2ratio = dict()
ratios = set()
rationame = dict()
isratio = dict()
ratioorder = dict()
denominators = defaultdict(list)
sfid2corrmat = dict()
for sf in psmdb.spectrumfiles():
    asname = sf.getdata('analyticalsample')
    if not asname:
        asname = ":".join(sf.sample_names())
    corrections = np.zeros(shape=(len(labelnames),9))
    correctioncount = 0
    for ssfl in sf.samples:
        sid = ssfl.sample.name
        label = ssfl.label.name
        fileidlabel2sampid[(sf.id,label)] = sid
    for ssfl in sf.samples:
        sid = ssfl.sample.name
        label = ssfl.label.name

        for denom,ordinal in ssfl.getdata('ratio_denominator',[]):
            denominators[(sf.id,label)].append(denom)
            dsid = fileidlabel2sampid.get((sf.id,denom),None)
            assert dsid != None or denom == "1.0"
            if dsid not in ("POOL",None):
                ratio = "%s/%s"%(sid,dsid)
            else:
                ratio = "%s"%(sid,)
            fileidlabel2ratio[(sf.id,label,denom)] = ratio
            rationame[ratio] = ratio
            isratio[ratio] = (denom != "1.0")
            ratios.add(ratio)
            ratioorder[ratio] = (asname,ordinal)

        correction = ssfl.getdata('isotope_corrections')
        if correction != None and not opts.nocorrections:
            assert opts.labels in ('TMT10','TMT10+2','TMT11','TMT11+2','TMT16','TMT16+2','TMT18')
            correction = list(map(float,correction))
            correctioncount += 1
            correction.insert(4,100.0)
            rowi = labelnames.index(label) #labelnames must be in correct order!
            corrections[rowi,:] = correction

    if correctioncount == len(labelnames):
        corrections /= corrections.sum(axis=1,keepdims=True)
        correctionmat = np.zeros(shape=(len(labelnames),len(labelnames)+len(extralabels)))
        for i in range(len(labelnames)):
            for k in range(0,9):
                j = i+(-4+k)
                if j >= 0 and j < correctionmat.shape[1]:
                    correctionmat[i,j] = corrections[i,k]
        sfid2corrmat[sf.id] = correctionmat.transpose()
    else:
        assert correctioncount == 0

nratios = len(ratios)

progress.done()

progress.stage("Determining protein and peptide lists for output")

if opts.peptide:
    peps = [ psmdb.peptideid(opts.peptide) ]
    pepids = [p.id for p in peps]
    prids = psmdb.proteinids()
elif opts.locus:
    pracc,sites = opts.locus.split(':',1)
    pr = psmdb.protein(accession=pracc)
    sites=list(map(int,sites.split(':')))
    peps = pr.peptides(loci=sites)
    pepids = [p.id for p in peps]
    prids = set([pr.id])
elif opts.protein:
    pr = psmdb.protein(accession=opts.protein)
    peps = pr.peptides()
    pepids = [p.id for p in peps]
    prids = set([pr.id])
elif opts.pepseqregex or opts.pepseqntermaas:
    if opts.pepseqntermaas:
        opts.pepseqregex = r"^[%s]"%(opts.pepseqntermaas,)
    pepseqregex = re.compile(opts.pepseqregex)
    pepids = set()
    for pep in psmdb.peptides():
        if pepseqregex.search(pep.sequence):
            pepids.add(pep.id)
    prids = psmdb.proteinids()
else:
    pepids = psmdb.peptideids()
    prids = psmdb.proteinids()

progress.done()
progress.message("Proteins: %d, peptides: %d."%(len(prids),len(pepids)))
if len(pepids) == 0:
    sys.exit(0)

def reporter_ion_extract(pepid,spectra,labelnames,corrections):
    bymz = defaultdict(lambda:tuple([0]+[None]*len(labelnames)))
    for s in spectra:
        if opts.noties and len(s.psms) > 1:
            continue
        correction = corrections.get(s.fileID)
        for psm in s.psms:
            ion = psm.peptideIon
            if ion.peptideID == pepid:
                key1 = ('mz',"%.2f"%(ion.mz,))
                key2 = ('ionid',ion.id)
                reporters = []
                for lab in labelnames + extralabels:
                    reporters.append(s.getdata('%s_reporter_ion_intensity'%lab,0))
                sumrep = sum(reporters)
                if sumrep <= 0:
                    continue
                if correction is not None:
                    reporters = np.array(reporters)
                    if opts.correction in ("LS","LS+Check"):
                        reporters1 = np.linalg.solve(correction,reporters)
                    elif opts.correction in ("NNLS","NNLS+Check"):
                        reporters1,rnorm1 = scipy.optimize.nnls(correction,reporters)
                    else:
                        raise RuntimeError('Bad correction? '+opts.correction)
                    if opts.correction in ("LS+Check","NNLS+Check"):
                        reporters2 = np.array(reporters1)
                        reporters2[reporters2<(max(reporters2)/1024)] = 0
                        reporters2[reporters[:len(labelnames)]<=0] = 0
                        expreps2 = np.dot(correction, reporters2)
                        rnorm2 = np.linalg.norm(expreps2 - reporters,2)
                        relrnorm2 = rnorm2/np.sum(reporters)
                        allclose2 = np.all((expreps2 - reporters)/np.max(reporters) < 0.05) and (relrnorm2 < 0.10)
                        if not allclose2:
                            continue
                        reporters = list(reporters2)
                    else:
                        reporters = list(reporters1)
                bymz[key1] = max(bymz[key1],tuple([sumrep]+reporters+[key1]))
                if not opts.fullylocalized or psm.getdata('CPTAC-CDAP:FullyLocalized') != 'N':
                    bymz[key2] = max(bymz[key2],tuple([sumrep]+reporters+[key2]))

    return list(bymz.values())

phosphomods = set([psmdb.modification('S',79.966),
                   psmdb.modification('T',79.966),
                   psmdb.modification('Y',79.966)])

glycomods = set([psmdb.modification('N',0.984)])

acetylmods = set([psmdb.modification('K',42.011)])

progress.stage("Retrieve reporter ions by file and peptide ion")
from functools import partial
ionvalues = psmdb.fileid_peptideid_specfunc(partial(reporter_ion_extract,
                                                    labelnames=labelnames,
                                                    corrections=sfid2corrmat))
progress.done()

def isdecoy(it):
    if hasattr(it,'decoy'):
        return (it.decoy == True)
    return (it.getdata('decoy',False) == True)


modsites = defaultdict(set)
if opts.report == 'PhosphopeptideSite':
    progress.stage("Retrieve potential modification sites for each peptide",len(pepids))
    for pepid in pepids:
        pep = psmdb.peptide(pepid);
        if opts.fullylocalized:
            modsites[pepid] = pep.modsites(phosphomods,(lambda psm: psm.getdata('CPTAC-CDAP:FullyLocalized') != 'N'))
        else:
            modsites[pepid] = pep.modsites(phosphomods)
        progress.update()
    progress.done()
elif opts.report == 'Phosphosite':
    progress.stage("Retrieve potential modification sites for each protein",2*len(pepids))
    prmodsites = defaultdict(set)
    for pepid in pepids:
        pep = psmdb.peptide(pepid);
        if opts.fullylocalized:
            mss = pep.modsites(phosphomods,(lambda psm: psm.getdata('CPTAC-CDAP:FullyLocalized') != 'N'))
        else:
            mss = pep.modsites(phosphomods)
        for pr in pep.proteins():
            als = pep.alignmentsto(pr)
            if len(als) > 1 or als[0].start in (-1,None):
                continue
            al = als[0]
            for ms in mss:
                prmodsites[pr.id].add(al.start+ms-1)
        progress.update()
    for pepid in pepids:
        pep = psmdb.peptide(pepid);
        for pr in pep.proteins():
            als = pep.alignmentsto(pr)
            if len(als) > 1 or als[0].start in (-1,None):
                continue
            al = als[0]
            for ms in prmodsites[pr.id]:
                site =  ms-al.start+1
                if site >= 1 and site <= len(pep.sequence):
                    modsites[pepid].add(ms-al.start+1)
        progress.update()
    progress.done()

ratios1 = defaultdict(list)
seenlogratios = set()
logratioitems = set()
logratioionids = defaultdict(set)
progress.stage("Average peptide ion iTRAQ ratios.",len(ionvalues))
modmap = dict(k='K',m='M',n='N',q='Q',c='C')
# keep Asn, deamidation
modmap1 = dict(k='K',m='M',q='Q',c='C')
#
modmap2 = {psmdb.modification('K',42.011): 'a'}
modmap3 = {psmdb.modification('K',114.043): 'u'}
for i,pepid in enumerate(pepids):
    pep = psmdb.peptide(pepid); pepseq = pep.sequence;
    if opts.specific and pep.getdata('tryptic-termini') != 'specific':
        continue
    if any(map(isdecoy,pep.proteins())) or any(map(isdecoy,pep.genes())):
        continue
    for sf in psmdb.spectrumfiles():
        if ionvalues[(sf.id,pepid)] == None:
            continue
        for ionintensities in ionvalues[(sf.id,pepid)]:
            type = ionintensities[-1][0]
            items = []
            ion = None
            if type == 'mz':
                # items.append(pepseq + ':'+ionintensities[-1][1])
                continue
            elif type == 'ionid':
                ion = psmdb.ion(ionintensities[-1][1])
                # items.append(':'.join([ion.modifiedPeptide(),"%d+"%ion.charge]))
                # items.append(':'.join([ion.modifiedPeptide()]))
                # if 'y' in ion.modifiedPeptide():
                #     continue
                # items.append(':'.join([ion.modifiedPeptide().replace('m','M').replace('k','K')]))
                # items.append(':'.join([pep.sequence]))
                # items.append(':'.join([ion.modifiedPeptide().replace('m','M')[-12:]]))
                # print ion.modString,ion.charge, ionintensities
                if opts.report == 'Peptide':
                    items.append((pepseq,pepid))
                elif opts.report == 'PhosphopeptideSite':
                    modpep = ion.modifiedPeptide()
                    modpep = ''.join([modmap.get(aa,aa) for aa in modpep])
                    # any = False
                    for i in modsites[pepid]:
                        psl = list(pepseq)
                        if modpep[i-1] in 'sty':
                            psl.insert(i,'[+80.0]')
                        else:
                            psl.insert(i,'[+0.00]')
                        items.append((''.join(psl),pepid,i))
                        # any = True
                    # if not any:
                    #     items.append((pepseq,pepid,None))
                elif opts.report == 'Phosphosite':
                    modpep = ion.modifiedPeptide()
                    modpep = ''.join([modmap.get(aa,aa) for aa in modpep])
                    for pr in pep.proteins():
                        if pr.id not in prids:
                            continue
                        if opts.specific:
                            anyspec = False
                            for al in pep.alignmentsto(pr):
                                if al.getdata('tryptic-termini',None) == 'specific':
                                    anyspec = True
                                    break
                            if not anyspec:
                                continue
                        als = pep.alignmentsto(pr)
                        if len(als) == 1 and als[0].start not in (-1,None):
                            st = als[0].start
                            for i in modsites[pepid]:
                                if modpep[i-1] in 'sty':
                                    it = ':'.join([pr.accession,"%s%d"%(pepseq[i-1],st+i-1),"+80.0"])
                                else:
                                    it = ':'.join([pr.accession,"%s%d"%(pepseq[i-1],st+i-1),"+0.00"])
                                items.append((it,pr.id,st+i-1))
                elif opts.report == 'PhosphositeCombination':
                    modpep = ion.modifiedPeptide()
                    modpep = ''.join([modmap.get(aa,aa) for aa in modpep])
                    for pr in pep.proteins():
                        if pr.id not in prids:
                            continue
                        for al in pep.alignmentsto(pr):
                            if opts.specific and al.getdata('tryptic-termini',None) != 'specific':
                                continue
                            if al.start in (-1,None):
                                continue
                            phossites = []
                            phossitepos = []
                            st = al.start
                            for i,aa in enumerate(modpep):
                                if aa in 'sty':
                                    phossites.append("%s%d"%(modpep[i],st+i))
                                    phossitepos.append(st+i)
                            if len(phossites) == 0:
                                continue
                            it = ':'.join([pr.accession,"".join(phossites)])
                            items.append((it,pr.id,",".join(map(str,phossitepos))))
                elif opts.report == 'AcetylsiteCombination':
                    modpep = ion.annotatedPeptide(modmap2,before=True)
                    modpep = modpep.replace('aK','k')
                    for pr in pep.proteins():
                        if pr.id not in prids:
                            continue
                        for al in pep.alignmentsto(pr):
                            if opts.specific and al.getdata('tryptic-termini',None) != 'specific':
                                continue
                            if al.start in (-1,None):
                                continue
                            acetsites = []
                            acetsitepos = []
                            st = al.start
                            for i,aa in enumerate(modpep):
                                if aa in 'k':
                                    acetsites.append("%s%d"%(modpep[i],st+i))
                                    acetsitepos.append(st+i)
                            if len(acetsites) == 0:
                                continue
                            it = ':'.join([pr.accession,"".join(acetsites)])
                            items.append((it,pr.id,",".join(map(str,acetsitepos))))
                elif opts.report == 'UbiquitylsiteCombination':
                    modpep = ion.annotatedPeptide(modmap3,before=True)
                    modpep = modpep.replace('uK','k')
                    for pr in pep.proteins():
                        if pr.id not in prids:
                            continue
                        for al in pep.alignmentsto(pr):
                            if opts.specific and al.getdata('tryptic-termini',None) != 'specific':
                                continue
                            if al.start in (-1,None):
                                continue
                            ubisites = []
                            ubisitepos = []
                            st = al.start
                            for i,aa in enumerate(modpep):
                                if aa in 'k':
                                    ubisites.append("%s%d"%(modpep[i],st+i))
                                    ubisitepos.append(st+i)
                            if len(ubisites) == 0:
                                continue
                            it = ':'.join([pr.accession,"".join(ubisites)])
                            items.append((it,pr.id,",".join(map(str,ubisitepos))))
                elif opts.report == 'Phosphopeptide':
                    modpep = ion.modifiedPeptide()
                    modpep = ''.join([modmap.get(aa,aa) for aa in modpep])
                    if re.search(r'[sty]',modpep):
                        items.append((modpep,pepid))
                elif opts.report == 'Acetylpeptide':
                    modpep = ion.annotatedPeptide(modmap2,before=True)
                    modpep = modpep.replace('aK','k')
                    if re.search(r'[k]',modpep):
                        items.append((modpep,pepid))
                elif opts.report == 'Ubiquitylpeptide':
                    modpep = ion.annotatedPeptide(modmap3,before=True)
                    modpep = modpep.replace('uK','k')
                    if re.search(r'[k]',modpep):
                        items.append((modpep,pepid))
                elif opts.report == 'Glycopeptide':
                    modpep = ion.modifiedPeptide()
                    for al in ion.peptide.alignments:
                        if opts.specific and al.getdata('tryptic-termini',None) != 'specific':
                            continue
                        if al.hasdata('right-context'):
                            modpep1 = modpep + al.getdata('right-context')
                        elif al.hasdata('right-amino-acid'):
                            modpep1 = modpep + al.getdata('right-amino-acid',)
                        else:
                            modpep1 = modpep
                        modpep1 += "XXXX"
                        modpep1 = [modmap1.get(aa,aa) for aa in modpep1]
                        for i in range(len(modpep1)):
                            if modpep1[i] == 'n':
                                if (modpep1[i+2] not in 'ST') or \
                                   (modpep1[i+1] in 'P') or \
                                   (modpep1[i+3] in 'P'):
                                    modpep1[i] = 'N'
                        modpep1 = ''.join(modpep1)
                        if 'n' in modpep1:
                            items.append((modpep1[:len(modpep)],pepid))
                            break
                elif opts.report == 'GlycositeCombination':
                    modpep = ion.modifiedPeptide()
                    for al in ion.peptide.alignments:
                        if opts.specific and al.getdata('tryptic-termini',None) != 'specific':
                            continue
                        if al.start in (-1,None):
                            continue
                        pr = al.protein
                        if pr.id not in prids:
                            continue
                        if al.hasdata('right-context'):
                            modpep1 = modpep + al.getdata('right-context')
                        elif al.hasdata('right-amino-acid'):
                            modpep1 = modpep + al.getdata('right-amino-acid')
                        else:
                            modpep1 = modpep
                        # so we can check at the end of the peptide without error
                        modpep1 += "XXXX"
                        modpep1 = [modmap1.get(aa,aa) for aa in modpep1]
                        for i in range(len(modpep)):
                            if modpep1[i] == 'n':
                                if (modpep1[i+2] not in 'ST') or \
                                   (modpep1[i+1] in 'P') or \
                                   (modpep1[i+3] in 'P'):
                                    modpep1[i] = 'N'
                        modpep1 = ''.join(modpep1)
                        if not 'n' in modpep1:
                            continue
                        glycosites = []
                        glycositepos = []
                        st = al.start
                        for i,aa in enumerate(modpep1[:len(modpep)]):
                            if aa == 'n':
                                glycosites.append("%s%d"%(modpep1[i],st+i))
                                glycositepos.append(st+i)
                        if len(glycosites) == 0:
                            continue
                        it = ':'.join([pr.accession,"".join(glycosites)])
                        items.append((it,pr.id,",".join(map(str,glycositepos))))
                elif opts.report == 'ModifiedPeptide':
                    items.append((ion.modifiedPeptide(),pepid))
                elif opts.report == 'PeptideIon':
                    items.append((":".join([ion.modifiedPeptide(),"%d+"%ion.charge]),pepid))
                elif opts.report == 'PeptideForm':
                    items.append((':'.join([pep.sequence,ion.modString]),pepid))
            weight = ionintensities[0]
            ionint = dict(list(zip(labelnames,ionintensities[1:])))
            fid = sf.id
            for lab in labelnames:
                numer = ionint[lab]
                if numer <= 0:
                    continue
                for denomlab in denominators[(fid,lab)]:
                    if denomlab != "1.0":
                        denom = ionint[denomlab]
                        if denom <= 0:
                            continue
                        logratio = log(numer/denom,2.0)
                    else:
                        logratio = numer

                    ratio = fileidlabel2ratio[(fid,lab,denomlab)]
                    seenlogratios.add(ratio)
                    for it in items:
                        logratioitems.add(it)
                        logratioionids[it].add(ion.id)
                        ratios1[(it,ratio)].append((weight,logratio))

        progress.update()
progress.done()

def cleanup(w,x):
        # print >>sys.stderr,w,x
    n0 = len(x)
    u = np.average(x)
    s = np.std(x)
    for i in range(len(x)-1,-1,-1):
        if x[i] > u+2*s or x[i] < u-2*s:
            del x[i]
            del w[i]
    n1 = len(x)
    return np.average(x),n0,n1

def cleanup1(w,x):
    # print >>sys.stderr,w,x
    n0 = len(x)
    return np.sum(x),n0,n0

progress.stage("Cleanup observations per peptide ion.",len(ratios1))
for k,v in list(ratios1.items()):
    w = list(map(itemgetter(0),v))
    x = list(map(itemgetter(1),v))
    if isratio[k[1]]:
        ratios1[k] = cleanup(w,x)
    else:
        ratios1[k] = cleanup1(w,x)
    progress.update()
progress.done()

if opts.allratios:
    progress.stage("Sorting logratios.")
    ratios = sorted(ratios,key=ratioorder.get)
else:
    progress.stage("Sorting observed logratios.")
    ratios = sorted(seenlogratios,key=ratioorder.get)
# print set(sampids1)-sampids
progress.done()

# headers = ['Item','Samples','OrigObs/Sample','TrimObs/Sample']

pep = psmdb.peptide(next(iter(pepids)))
pr = pep.proteins()[0]

if opts.report == 'Peptide':
    keyheader = 'Peptide'
elif opts.report == 'ModifiedPeptide':
    keyheader = 'ModifiedPeptide'
elif opts.report == 'PeptideIon':
    keyheader = 'PeptideIon'
elif opts.report == 'Phosphopeptide':
    keyheader = 'Phosphopeptide'
elif opts.report == 'Glycopeptide':
    keyheader = 'Glycopeptide'
elif opts.report == 'Acetylpeptide':
    keyheader = 'Acetylpeptide'
elif opts.report == 'Ubiquitylpeptide':
    keyheader = 'Ubuquitylpeptide'
elif opts.report == 'PhosphopeptideSite':
    keyheader = 'PhosphopeptideSite'
elif opts.report == 'Phosphosite':
    keyheader = 'Phosphosite'
elif opts.report == 'PhosphositeCombination':
    keyheader = 'Phosphosite'
elif opts.report == 'GlycositeCombination':
    keyheader = 'Glycosite'
elif opts.report == 'AcetylsiteCombination':
    keyheader = 'Acetylsite'
elif opts.report == 'UbiquitylsiteCombination':
    keyheader = 'Ubiquitylsite'
elif opts.report == 'PeptideForm':
    keyheader = 'PeptideForm'

def ratio2header(s):
    if '/' in rationame[s]:
        return "Log %s"%(rationame[s],)
    if not isratio[s]:
        return (rationame[s]+" Summed Intensity")
    return (rationame[s]+" Log Ratio")

headers = [keyheader]
for s in ratios:
    headers.extend([ratio2header(s)])
if opts.extraheaders:
    headers.extend(['Samples','OrigObs/Sample','TrimObs/Sample'])

def sitesortkey(s):
    acc,sites = s.split(':')
    sites = list(map(int,[_f for _f in re.split(r'[^\d]+',sites) if _f]))
    return acc,sites[0],sites[-1]

sortkey = lambda s: s
if opts.report in ('PhosphositeCombination','GlycositeCombination','AcetylsiteCombination','UbiquitylsiteCombination'):
    sortkey = sitesortkey

prisgene=False
if pr.hasdata('geneid') or not pr.hasdata('source'):
    headers.extend(['Gene'])
    prisgene=True
elif opts.report in ('Phosphosite','PhosphositeCombination','GlycositeCombination','AcetylsiteCombination','UbiquitylsiteCombination') :
    headers.extend(['Peptide','Gene'])
else:
    headers.extend(['Protein','Gene'])
if pr.hasdata('organism'):
    headers.append('Organism')

def rows(ratios,logratioitems):
    progress.stage("Construct output table rows",len(logratioitems))
    for key in sorted(logratioitems,key=lambda t: sortkey(t[0])):
        row = {keyheader: key[0]}
        nsamps = 0
        nobsorig = 0.
        nobstrim = 0.

        for s in ratios:
            if (key,s) in ratios1:
                row[ratio2header(s)] = ratios1[(key,s)][0]
                nobsorig += ratios1[(key,s)][1]
                nobstrim += ratios1[(key,s)][2]
                nsamps += 1
            else:
                row[ratio2header(s)] = ""
        if nsamps == 0:
            continue
        row['Samples'] = nsamps
        row['OrigObs/Sample'] = nobsorig/nsamps
        row['TrimObs/Sample'] = nobstrim/nsamps
        genes = []
        if opts.report not in ('Phosphosite','PhosphositeCombination','GlycositeCombination','AcetylsiteCombination','UbiquitylsiteCombination'):
            proteins = []
            orgs = []
            pep = psmdb.peptide(key[1])
            for pr in pep.proteins():
                if opts.report == 'PhosphopeptideSite' and key[2] != None:
                    als = pep.alignmentsto(pr)
                    if len(als) == 1 and als[0].start not in (-1,None):
                        proteins.append(':'.join([pr.accession,"%d"%(key[2]+als[0].start-1)]))
                    else:
                        proteins.append(pr.accession)
                else:
                    proteins.append(pr.accession)
                if pr.getdata('organism','') not in orgs:
                    orgs.append(pr.getdata('organism',''))
                for prg in pr.groups():
                    if prg.type != 'Gene':
                        continue
                    if prg.name not in genes:
                        genes.append(prg.name)
        elif opts.report == 'Phosphosite':
            pr = psmdb.protein(key[1])
            peptides = []
            for pep in pr.peptides(locus=key[2]):
                als = pep.alignmentsto(pr)
                if len(als) == 1 and als[0].start not in (-1,None):
                    peppos = key[2]-als[0].start+1
                    peptides.append((als[0].start,als[0].end,pep.sequence,peppos))
            peptides = ["%s:%s"%(it[2],it[3]) for it in sorted(peptides)]
            org = pr.getdata('organism','')
            for prg in pr.groups():
                if prg.type != 'Gene':
                    continue
                if prg.name not in genes:
                    genes.append(prg.name)
        elif opts.report in ('PhosphositeCombination','GlycositeCombination','AcetylsiteCombination','UbiquitylsiteCombination'):
            pr = psmdb.protein(key[1])
            prosites = list(map(int,key[2].split(',')))
            peptides = set()
            for ionid in logratioionids[key]:
                ion = psmdb.ion(ionid)
                modpep = ion.modifiedPeptide()
                if opts.report == 'PhosphositeCombination':
                    modpep = ''.join([modmap.get(aa,aa) for aa in modpep])
                elif opts.report == 'AcetylsiteCombination':
                    modpep = ion.annotatedPeptide(modmap2,before=True)
                    modpep = modpep.replace('aK','k')
                elif opts.report == 'UbiquitylsiteCombination':
                    modpep = ion.annotatedPeptide(modmap3,before=True)
                    modpep = modpep.replace('uK','k')
                for al in ion.peptide.alignmentsto(pr,loci=prosites):
                    if opts.specific:
                        assert (al.getdata('tryptic-termini',None) == 'specific')
                    if opts.report == 'GlycositeCombination':
                        if al.hasdata('right-context'):
                            modpep1 = modpep + al.getdata('right-context')
                        elif al.hasdata('right-amino-acid'):
                            modpep1 = modpep + al.getdata('right-amino-acid')
                        else:
                            modpep1 = modpep
                        modpep1 = [modmap1.get(aa,aa) for aa in modpep1]
                        for i in range(len(modpep1)-2):
                            if modpep1[i] == 'n' and modpep1[i+2] not in 'ST':
                                modpep1[i] = 'N'
                        modpep = ''.join(modpep1)[:len(modpep)]
                    peptides.add((al.start,al.end,modpep))
                    break
            peptides = [it[2] for it in sorted(peptides)]
            org = pr.getdata('organism','')
            for prg in pr.groups():
                if prg.type != 'Gene':
                    continue
                if prg.name not in genes:
                    genes.append(prg.name)
        if prisgene:
            row['Gene'] = ';'.join(sorted(proteins))
        elif opts.report in ('Phosphosite','PhosphositeCombination','GlycositeCombination','AcetylsiteCombination','UbiquitylsiteCombination'):
            row['Peptide'] = ';'.join(peptides)
            row['Gene'] = ';'.join(genes)
            row['Organism'] = org
        else:
            row['Protein'] = ';'.join(proteins)
            row['Gene'] = ';'.join(genes)
            row['Organism'] = ';'.join(orgs)
        if nsamps >= opts.ratio*len(ratios) and nsamps >= opts.count:
            yield row
        progress.update()
    progress.done()

if opts.output:

    from dataset import CSVFileTable, XLSXFileTable, XLSFileTable, TSVFileTable
    if opts.output.endswith('.csv'):
        t = CSVFileTable(filename=opts.output,headers=headers,append=opts.append)
        t.from_rows(rows(ratios,logratioitems))
    elif opts.output.endswith('.tsv'):
        t = TSVFileTable(filename=opts.output,headers=headers,append=opts.append)
        t.from_rows(rows(ratios,logratioitems))
    elif opts.output.endswith('.xls'):
        t = XLSFileTable(filename=opts.output,headers=headers,sheet='Log Ratios')
        t.from_rows(rows(retios,logratioitems))
    elif opts.output.endswith('.xlsx'):
        t = XLSXFileTable(filename=opts.output,headers=headers,sheet='Log Ratios')
        t.from_rows(rows(ratios,logratioitems))
