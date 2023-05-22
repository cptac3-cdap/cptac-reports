#!/usr/bin/env python3

import sys, os
from collections import defaultdict
from operator import attrgetter

import atexit
def onexit():
    print >>sys.stderr, "Type <Enter> to Exit...",
    sys.stderr.flush()
    raw_input()

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
# proteins = psmdb.proteins()
# proids = set(psmdb.proteinids())
pepids = set(psmdb.peptideids())
# unshpepids = psmdb.unshared_peptideids()
# unshpepids.remove(psmdb.grouped_peptideids())
progress.done()

# sampids = set()
samplekey = dict()
fileid2sampid = dict()
for sf in psmdb.spectrumfiles():
    sid = sf.getdata('analyticalsample')
    sidord = sf.getdata('analyticalsampleordinal',0)
    if not sid:
        sid = ":".join(sf.sample_names())
    fileid2sampid[sf.id] = sid
    samplekey[sid] = (sidord,sid)

pep = next(iter(psmdb.peptides()))
pr = pep.proteins()[0]
prisgene=False
if pr.hasdata('geneid'):
    headers = ['Peptide','Charge','Mods','MinFDR','SpectralCount','AmbigSpectralCount','Sample','Gene','Assay']
    prisgene=True
else:
    headers = ['Peptide','Charge','Mods','MinFDR','SpectralCount','AmbigSpectralCount','Sample','Protein','Gene','Assay']

def isdecoy(it):
    if hasattr(it,'decoy'):
        return (it.decoy == True)
    return (it.getdata('decoy',False) == True)

def rows(fileid2sampid):
    progress.stage("Output peptide summary",len(pepids))
    for pep in psmdb.peptides():
        if isdecoy(pep):
            continue
        charges = set()
        mods = set()
        samples = set()
        spectra = set()
        genes = set()
        proteins = set()
        minfdr = 1.0
        for ion in pep.ions:
            charges.add(ion.charge)
            mods.add("-" if not ion.modString else ion.modString)
            for psm in ion.psms:
                samples.add(fileid2sampid[psm.spectrum.fileID])
                spectra.add(psm.spectrum)
                minfdr = min(minfdr,psm.value)
        nties = 0
        for spec in spectra:
            if len(spec.distinct_peptides()) > 1:
                nties += 1
        for pr in pep.proteins():
            if isdecoy(pr):
                continue
            proteins.add(pr.accession)
            for prg in pr.groups():
                if prg.type != 'Gene':
                    continue
                if isdecoy(prg):
                    continue
                genes.add(prg.name)
        row = dict(Peptide=pep.sequence)
        row['Charge'] = ';'.join(sorted(map(str,charges)))
        row['Mods'] = ';'.join(sorted(mods))
        row['MinFDR'] = minfdr
        row['Sample'] = ';'.join(sorted(samples,key=samplekey.get))
        row['SpectralCount'] = len(spectra)
        row['AmbigSpectralCount'] = nties
        if prisgene:
            row['Gene'] = ';'.join(sorted(proteins))
        else:
            row['Protein'] = ';'.join(sorted(proteins))
            row['Gene'] = ';'.join(sorted(genes))
        row['Assay'] = ';'.join(sorted(pep.getdata('Assay',[]),key=lambda a: int(a.rsplit('-',1)[1])))
        # for key in ("Protein","Gene","Assay"):
        #     if row[key] == "":
        #       row[key] = "-"
        yield row
        progress.update()
    progress.done()

from dataset import CSVFileTable, TSVFileTable, XLSXFileTable, XLSFileTable
if opts.output.endswith('.csv'):
    t = CSVFileTable(filename=opts.output,headers=headers,compression=None)
elif opts.output.endswith('.tsv'):
    t = TSVFileTable(filename=opts.output,headers=headers,compression=None)
elif opts.output.endswith('.xls'):
    t = XLSFileTable(filename=opts.output,headers=headers,sheet='Peptides')
elif opts.output.endswith('.xlsx'):
    t = XLSXFileTable(filename=opts.output,headers=headers,sheet='Peptides')
t.from_rows(rows(fileid2sampid))
