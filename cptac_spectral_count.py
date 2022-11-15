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

from psmdb.PSMDb import PSMDb

parser.add_option("-d","--database",type="file",dest="database",remember=True,
                  help="PSM database file",default="",notNone=True,
                  name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("--bygene",action="store_true",dest="bygene",default=False,remember=True,
                  name="By gene",help="Gene-based summary. Default: False.")
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
    sid = ":".join(sf.sample_names())
    sampids.add(sid)
    fileid2sampid[sf.id] = sid

nsamp = len(sampids)

progress.stage("Retrieve PSM counts by file and peptide")
pepcnts = psmdb.fileid_peptideid_psmcounts()

counts1 = defaultdict(int)
counts2 = defaultdict(int)
progress.stage("Sum PSM counts of peptides by protein and sample",len(pepcnts))
for (fid,pepid),cnt in pepcnts.items():
    if fid in fileid2sampid:
        if opts.bygene:
            it = psmdb.peptide(pepid).genes()
        else:
            it = psmdb.peptide(pepid).proteins()
        for pr in it:
            counts2[(pr.id,fileid2sampid[fid])] += cnt
            counts2[(pr.id,None)] += cnt
            if pepid in unshpepids:
                counts1[(pr.id,fileid2sampid[fid])] += cnt
                counts1[(pr.id,None)] += cnt
        counts2[(None,fileid2sampid[fid])] += cnt
        counts2[(None,None)] += cnt
        if pepid in unshpepids:
            counts1[(None,fileid2sampid[fid])] += cnt
            counts1[(None,None)] += cnt
    progress.update()
progress.done()

sampids = sorted(sampids)

if not opts.bygene:
    acchdr = "Protein"
else:
    acchdr = "Gene"
headers = [acchdr]

for s in sampids:
    headers.extend([s+" Spectral Counts",s+" Unshared Spectral Counts"])

headers.extend(["Total Spectral Counts","Total Unshared Spectral Counts"])

headers.append('NCBIGeneID')
headers.append('Authority')
pr = next(iter(proteins))
for key in ('description','organism','chromosome','locus'):
    if pr.hasdata(key):
        headers.append(key.title())

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
            row[s+" Unshared Spectral Counts"] = counts1[(pr.id,s)]
            row[s+" Spectral Counts"] = counts2[(pr.id,s)]
        row["Total Unshared Spectral Counts"] = counts1[(pr.id,None)]
        row["Total Spectral Counts"] = counts2[(pr.id,None)]
        for k in pr.metadata:
            row[k.title()] = pr.getdata(k,"")
        row['NCBIGeneID'] = pr.getdata('geneid',"")
        if pr.hasdata('hgncid'):
            row['Authority'] = pr.getdata('hgncid',"")
        elif pr.hasdata('authid'):
            row['Authority'] = pr.getdata('authid',"")
        else:
            row['Authority'] = ""
        yield row
    row = dict(Gene="Total")
    for s in sampids:
        row[s+" Spectral Counts"] = counts2[(None,s)]
        row[s+" Unshared Spectral Counts"] = counts1[(None,s)]
    row["Total Spectral Counts"] = counts2[(None,None)]
    row["Total Unshared Spectral Counts"] = counts1[(None,None)]
    yield row

if opts.output:

    from dataset import CSVFileTable, XLSXFileTable, XLSFileTable, TSVFileTable
    if opts.output.endswith('.csv'):
        t = CSVFileTable(filename=opts.output,headers=headers)
        t.from_rows(rows(sampids,proteins))
    elif opts.output.endswith('.tsv'):
        t = TSVFileTable(filename=opts.output,headers=headers)
        t.from_rows(rows(sampids,proteins))
    elif opts.output.endswith('.xls'):
        t = XLSFileTable(filename=opts.output,headers=headers,sheet='Spectral Counts')
        t.from_rows(rows(sampids,proteins))
    elif opts.output.endswith('.xlsx'):
        t = XLSXFileTable(filename=opts.output,headers=headers,sheet='Spectral Counts')
        t.from_rows(rows(sampids,proteins))
