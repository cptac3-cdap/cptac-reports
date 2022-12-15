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
parser.add_option("--bygene",action="store_true",dest="bygene",default=False,remember=True,
                  name="By gene",help="Gene-based summary. Default: False.")
parser.add_option("--labels",type="choice",dest="labels",choices=["TMT10","TMT10+2","TMT11","TMT11+2","TMT16","TMT16+2","TMT18","TMT6","iTRAQ","iTRAQ4","iTRAQ8"],default="TMT10",
                  help="Labeling reagent type. One of TMT6, TMT10, TMT10+2, TMT11, TMT11+2, TMT16, TMT16+2, TMT18, iTRAQ, iTRAQ4, iTRAQ8. Default: TMT10.")
parser.add_option("--nocorrections",action="store_true",dest="nocorrections",default=False,
                  help="Do not apply isotope corrections even if available. Default: False.")
parser.add_option("--correction",type="choice",dest="correction",default="NNLS+Check",
                  help="Strategy for computing corrected isotopic-label intensity values. One of None, LS, LS+Check, NNLS, NNLS+Check. Default: NNLS+Check.",
                  choices=["None","LS","LS+Check","NNLS","NNLS+Check"])
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

progress.stage("Retrieve proteins, peptides, and samples")

if not opts.bygene:
    proteins = psmdb.proteins()
    proids = psmdb.proteinids()
else:
    proteins = psmdb.genes()
    proids = psmdb.geneids()

pepids = psmdb.peptideids()
if not opts.bygene:
    unshpepids = psmdb.unshared_peptideids()
else:
    unshpepids = psmdb.unshared_peptideids_bygene()
unshpepids.difference_update(psmdb.grouped_peptideids())

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
    asord = sf.getdata('analyticalsampleordinal',0)
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
            ratioorder[ratio] = (asord,asname,ordinal)

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

nratios=len(ratios)

def reporter_ion_extract(pepid,spectra,labelnames,corrections):
    bymz = defaultdict(lambda:tuple([0]+[None]*len(labelnames)))
    for s in spectra:
        correction = corrections.get(s.fileID)
        for psm in s.psms:
            ion = psm.peptideIon
            if ion.peptideID == pepid:
                key = "%.2f"%(ion.mz,)
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
                bymz[key] = max(bymz[key],tuple([sumrep]+reporters))
    return list(bymz.values())

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

ratios1 = defaultdict(list)
ratios2 = defaultdict(list)
seenlogratios = set()
progress.stage("Average peptide ion ratios over protein and sample",len(ionvalues))
for (fid,pepid),values in ionvalues.items():
    # print fid,pepid,values
    for ionintensities in values:
        weight = ionintensities[0]
        ionint = dict(list(zip(labelnames,ionintensities[1:])))
        pep = psmdb.peptide(pepid)
        if not opts.bygene:
            aggrto = list(pep.proteins())
        else:
            aggrto = list(pep.genes())
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
                for it in aggrto:
                    ratios2[(it.id,ratio)].append((weight,logratio))
                    if pep.id in unshpepids:
                        ratios1[(it.id,ratio)].append((weight,logratio))

    progress.update()
progress.done()

def cleanup(w,x):
    # print >>sys.stderr,w,x
    u = np.average(x)
    s = np.std(x)
    for i in range(len(x)-1,-1,-1):
        if x[i] > u+2*s or x[i] < u-2*s:
            del x[i]
            del w[i]
    return np.average(x)

def cleanup1(w,x):
    # print >>sys.stderr,w,x
    return np.sum(x)

for k,v in list(ratios2.items()):
    w = list(map(itemgetter(0),v))
    x = list(map(itemgetter(1),v))
    if isratio[k[1]]:
        ratios2[k] = cleanup(w,x)
    else:
        ratios2[k] = cleanup1(w,x)

for k,v in list(ratios1.items()):
    w = list(map(itemgetter(0),v))
    x = list(map(itemgetter(1),v))
    if isratio[k[1]]:
        ratios1[k] = cleanup(w,x)
    else:
        ratios1[k] = cleanup1(w,x)

for rid in seenlogratios:
    x = []
    for prid in proids:
        if (prid,rid) in ratios2:
            x.append(ratios2[(prid,rid)])
    ratios2[("Mean",rid)] = np.mean(x)
    ratios2[("StdDev",rid)] = np.std(x)
    ratios2[("Median",rid)] = np.median(x)
    for prid in proids:
        if (prid,rid) in ratios2:
            if isratio[rid]:
                ratios2[(prid,rid)] -= ratios2[("Median",rid)]
    x = []
    for prid in proids:
        if (prid,rid) in ratios1:
            x.append(ratios1[(prid,rid)])
    ratios1[("Mean",rid)] = np.mean(x)
    ratios1[("StdDev",rid)] = np.std(x)
    ratios1[("Median",rid)] = np.median(x)
    for prid in proids:
        if (prid,rid) in ratios1:
            if isratio[rid]:
                ratios1[(prid,rid)] -= ratios1[("Median",rid)]

ratios = sorted(seenlogratios,key=ratioorder.get)

def ratio2header(s):
    if '/' in rationame[s]:
        return ["Log %s"%(rationame[s],),"Unshared Log %s"%(rationame[s],)]
    if not isratio[s]:
        return [(rationame[s]+" Summed Intensity"),(rationame[s]+" Unshared Summed Intensity")]
    return [(rationame[s]+" Log Ratio"),(rationame[s]+" Unshared Log Ratio")]

pr = next(iter(proteins))
acchdr='Protein'
acckey='accession'
if pr.hasdata('geneid') or opts.bygene:
    if opts.bygene:
        acckey = 'name'
    acchdr='Gene'
headers = [acchdr]
for s in ratios:
    headers.extend(ratio2header(s))
if pr.hasdata('geneid'):
    headers.append('NCBIGeneID')
if pr.hasdata('hgncid') or pr.hasdata('authid'):
    headers.append('Authority')
for key in ('description','organism','chromosome','locus','source'):
    if pr.hasdata(key):
        headers.append(key.title())
# headers.insert(len(headers)-1,'Gene')

def rows(ratios,proteins):
    row = {acchdr:"Mean"}
    for s in ratios:
        lrhdr,unshlrhrd = ratio2header(s)
        row[unshlrhrd] = ratios1.get(("Mean",s),"")
        row[lrhdr] = ratios2.get(("Mean",s),"")
    yield row
    row = {acchdr:"Median"}
    for s in ratios:
        lrhdr,unshlrhrd = ratio2header(s)
        row[unshlrhrd] = ratios1.get(("Median",s),"")
        row[lrhdr] = ratios2.get(("Median",s),"")
    yield row
    row = {acchdr:"StdDev"}
    for s in ratios:
        lrhdr,unshlrhdr = ratio2header(s)
        row[unshlrhdr] = ratios1.get(("StdDev",s),"")
        row[lrhdr] = ratios2.get(("StdDev",s),"")
    yield row
    for pr in sorted(proteins,key=attrgetter(acckey)):
        if isdecoy(pr):
            continue
        row = {acchdr:getattr(pr,acckey)}
        anyvalues = False
        for s in ratios:
            lrhdr,unshlrhdr = ratio2header(s)
            row[unshlrhdr] = ratios1.get((pr.id,s),"")
            row[lrhdr] = ratios2.get((pr.id,s),"")
            if not anyvalues and (row[lrhdr] or row[unshlrhdr]):
                anyvalues = True
        if not anyvalues:
            continue
        if pr.metadata:
            for k in pr.metadata:
                if k in ('gene',):
                    continue
                row[k.title()] = pr.getdata(k,"")
        row['NCBIGeneID'] = pr.getdata('geneid',"")
        if pr.hasdata('hgncid'):
            row['Authority'] = pr.getdata('hgncid',"")
        elif pr.hasdata('authid'):
            row['Authority'] = pr.getdata('authid',"")
        else:
            row['Authority'] = ""
        if acchdr == 'Protein':
            row['Gene'] = ';'.join(sorted([pg.name for pg in pr.groups() if pg.type == 'Gene']))
        yield row

if opts.output:

    from dataset import CSVFileTable, XLSXFileTable, XLSFileTable, TSVFileTable
    if opts.output.endswith('.csv'):
        t = CSVFileTable(filename=opts.output,headers=headers)
        t.from_rows(rows(ratios,proteins))
    elif opts.output.endswith('.tsv'):
        t = TSVFileTable(filename=opts.output,headers=headers)
        t.from_rows(rows(ratios,proteins))
    elif opts.output.endswith('.xls'):
        t = XLSFileTable(filename=opts.output,headers=headers,sheet='Results')
        t.from_rows(rows(ratios,proteins))
    elif opts.output.endswith('.xlsx'):
        t = XLSXFileTable(filename=opts.output,headers=headers,sheet='Results')
        t.from_rows(rows(ratios,proteins))
