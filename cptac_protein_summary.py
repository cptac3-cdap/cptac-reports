#!/usr/bin/env python27

from dataset import CSVFileTable, XLSFileTable, XLSFileDataset, ZIPFileDataset, ParsimonyCSV
import sys, os, os.path, copy, re, math, gzip, bz2
from collections import defaultdict
from operator import itemgetter, attrgetter
from functools import partial
from psmdb.PepXML import PepXMLReader, ParsimonyPepXML
from peptidescan.PeptideRemapper import *
from peptidescan.OutOfCoreTable import *

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

from psmdb.PSMDb import PSMDb, model

parser.add_option("-d","--database",type="file",dest="database",remember=True,
                  help="PSM database file",default="",notNone=True,
                  name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("--bygene",action="store_true",dest="bygene",default=False,remember=True,
                  name="By gene",help="Gene-based summary. Default: False.")
parser.add_option("-q","--quiet",action="store_true",dest="quiet",default=False,remember=True,
                  name="Quiet",help="Quiet. Default: False.")
parser.add_option("-o","--output",type="savefile",dest="output",default="",remember=True,
                  help="Output file.",
                  name="Output",
                  filetypes=[("Summary Table","*.csv;*.tsv;*.xlsx;*.xls"),
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

sampids = set()
fileid2sampid = dict()
for sf in psmdb.spectrumfiles():
    sid = sf.getdata('analyticalsample')
    if not sid:
        sid = ":".join(sf.sample_names())
    sampids.add(sid)
    fileid2sampid[sf.id] = sid

# print sampids

progress.stage("Retrieve PSM counts by file and peptide")
filepepcnts = psmdb.fileid_peptideid_psmcounts()
pepcnts = psmdb.peptideid_psmcounts()

counts1 = defaultdict(int) # Distinct peptides
counts2 = defaultdict(int) # Spectral counts per peptide
counts3 = defaultdict(int) # Unshared peptides
progress.stage("Compute PSM, and distinct and unshared peptide counts by sample and peptide",
               len(filepepcnts)+len(pepcnts))
for (fid,pepid),cnt in filepepcnts.items():
    if fid in fileid2sampid:
        counts1[(pepid,fileid2sampid[fid])] = 1
        counts2[(pepid,fileid2sampid[fid])] += cnt
        if pepid in unshpepids:
            counts3[(pepid,fileid2sampid[fid])] = 1
    progress.update()
for pepid,cnt in pepcnts.items():
    counts1[(pepid,None)] = 1
    counts2[(pepid,None)] += cnt
    if pepid in unshpepids:
        counts3[(pepid,None)] = 1
    progress.update()
progress.done()


## progress.stage("Retrieve min FDR by file and peptide")
## filepepminfdr = psmdb.fileid_peptideid_psmminvalue()
## progress.stage("Compute per sample min peptide FDR",len(filepepminfdr))
## minfdr = defaultdict(lambda:1.0)
## minnzfdr = 1e+20
## for (fid,pepid),minval in filepepminfdr.iteritems():
##     if fid in fileid2sampid:
##         minfdr[(pepid,fileid2sampid[fid])] = min(minfdr[(pepid,fileid2sampid[fid])],minval)
##     minfdr[(pepid,None)] = min(minfdr[(pepid,None)],minval)
##     if minval > 0:
##         minnzfdr = min(minnzfdr,minval)
##     progress.update()
## progress.done()

## def protfdr(pepfdrs,pseudofdr):
##     bestlogfdr = 0
##     for i,fdr in enumerate(sorted(pepfdrs)):
##         if fdr == 0:
##          fdr = pseudofdr
##         bestlogfdr = min(bestlogfdr,(i+1)*math.log(fdr))
##     return math.exp(bestlogfdr)

progress.stage("Compute protein metrics",proteins.count())
pcounts1 = defaultdict(int)
pcounts2 = defaultdict(int)
pcounts3 = defaultdict(int)
pcov = defaultdict(float)
pfdr = defaultdict(float)
for pr in proteins:
    prid = pr.id
    prpep = set(pr.peptides())
    for sid in sampids:
        pcounts1[(prid,sid)] = sum([counts1[(pep.id,sid)] for pep in prpep])
        pcounts2[(prid,sid)] = sum([counts2[(pep.id,sid)] for pep in prpep])
        pcounts3[(prid,sid)] = sum([counts3[(pep.id,sid)] for pep in prpep])
        thepepids = set(pep.id for pep in prpep if counts1[(pep.id,sid)]>0)
        theunshpepids = set(pep.id for pep in prpep if counts3[(pep.id,sid)]>0)
        # pfdr[(prid,sid)] = protfdr(map(lambda pepid: minfdr[(pepid,sid)],theunshpepids),minnzfdr/2.0)
        if len(thepepids) > 0 and not opts.bygene:
            try:
                pcov[(prid,sid)] = psmdb.proteinid_coverage(prid,thepepids)
            except RuntimeError:
                pcov[(prid,sid)] = 0
        else:
            pcov[(prid,sid)] = 0
    pcounts1[(prid,None)] = sum([counts1[(pep.id,None)] for pep in prpep])
    pcounts2[(prid,None)] = sum([counts2[(pep.id,None)] for pep in prpep])
    pcounts3[(prid,None)] = sum([counts3[(pep.id,None)] for pep in prpep])
    # pfdr[(prid,None)] = protfdr(map(lambda pep: minfdr[(pep.id,None)],
    #                               filter(lambda pep: pep.id in unshpepids,prpep)),
    #                           minnzfdr/2.0)
    if not opts.bygene:
        try:
            pcov[(prid,None)] = psmdb.proteinid_coverage(prid)
        except RuntimeError:
            pcov[(prid,None)] = 0
    progress.update()
progress.done()

sampids = sorted(sampids)

# headers = ['Gene','GroupID']
if not opts.bygene:
    acchdr = "Protein"
else:
    acchdr = "Gene"
headers = [acchdr]

for s in sampids:
    headers.extend([s+" Spectral Counts",
                    s+" Distinct Peptides",
                    s+" Unshared Peptides"])

headers.extend(["Spectral Counts",
                "Distinct Peptides",
                "Unshared Peptides"])

headers.append('NCBIGeneID')
headers.append('Authority')

pr = next(iter(proteins))
for key in ('description','organism','chromosome','locus'):
    if pr.hasdata(key):
        headers.append(key.title())
if opts.bygene:
    headers.append('Proteins')
headers.append('Assays')

acckey = 'accession'
if opts.bygene:
    acckey = 'name'

def isdecoy(it):
    if hasattr(it,'decoy'):
        return (it.decoy == True)
    return (it.getdata('decoy',False) == True)

def rows(sampids,proteins):
    for pr in sorted(proteins,key=attrgetter(acckey)):
        if isdecoy(pr):
            continue
        row = dict()
        row[acchdr] = getattr(pr,acckey)
        for s in sampids:
            row[s+" Spectral Counts"] = pcounts2[(pr.id,s)]
            row[s+" Distinct Peptides"] = pcounts1[(pr.id,s)]
            row[s+" Unshared Peptides"] = pcounts3[(pr.id,s)]
            # row[s+" Gene FDR"] = pfdr[(pr.id,s)]
        row["Spectral Counts"] = pcounts2[(pr.id,None)]
        row["Distinct Peptides"] = pcounts1[(pr.id,None)]
        row["Unshared Peptides"] = pcounts3[(pr.id,None)]
        # row["Gene FDR"] = pfdr[(pr.id,None)]
        if pr.metadata:
            for k in pr.metadata:
                if k:
                    row[k.title()] = pr.getdata(k,"")
        if not opts.bygene:
            row['Proteins'] = row['Mapped']
        else:
            row['Proteins'] = ";".join(sorted(map(attrgetter("accession"),pr.proteins())))

        row['NCBIGeneID'] = pr.getdata('geneid',"")
        if pr.hasdata('hgncid'):
            row['Authority'] = pr.getdata('hgncid',"")
        elif pr.hasdata('authid'):
            row['Authority'] = pr.getdata('authid',"")
        else:
            row['Authority'] = ""
        assays = []
        for pep in pr.peptides():
            assays.extend(pep.getdata('Assay',[]))
        row['Assays'] = ';'.join(sorted(assays,key=lambda a: int(a.rsplit('-',1)[1])))
        # print row
        yield row

from dataset import CSVFileTable, XLSXFileTable, XLSFileTable, TSVFileTable
if opts.output.endswith('.csv'):
    t = CSVFileTable(filename=opts.output,headers=headers)
elif opts.output.endswith('.tsv'):
    t = TSVFileTable(filename=opts.output,headers=headers)
elif opts.output.endswith('.xls'):
    t = XLSFileTable(filename=opts.output,headers=headers,sheet=acchdr+'s')
elif opts.output.endswith('.xlsx'):
    t = XLSXFileTable(filename=opts.output,headers=headers,sheet=acchdr+'s')
t.from_rows(rows(sampids,proteins))
