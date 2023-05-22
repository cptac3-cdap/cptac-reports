#!/usr/bin/env python3

from dataset import CSVFileTable, XLSFileTable, XLSFileDataset, ZIPFileDataset, ParsimonyCSV
import sys, os, os.path, copy, re, math, gzip, bz2
from collections import defaultdict
from operator import itemgetter
from functools import partial
from psmdb.PepXML import PepXMLReader, ParsimonyPepXML
from peptidescan.PeptideRemapper import *
from peptidescan.OutOfCoreTable import *

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

from psmdb.filter import SpectralFDR,PeptideFDR,PeptideProphet,EValue,NoFilter
from psmdb.PSMDb import PSMDb

psmmetric_choices = ["Estimated PSM FDR % (PepArML, CPTAC-CAP, MSGF+)",
                     "Estimated Peptide FDR % (PepArML, CPTAC-CAP, MSGF+)",
                     "PSM Probability (Peptide Prophet)",
                     "PSM E-value (PepArML, CPTAC-CAP, MSGF+)",
		     "None"]
psmmetric_default = 0
psmmetric_map = dict(zip(psmmetric_choices,
                         [SpectralFDR,PeptideFDR,PeptideProphet,EValue,NoFilter]))

parser.add_option("-d","--database",type="savefile",dest="database",remember=True,
	          help="PSM database file",default="",
	          name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-i","--input",type="files",dest="input",default="",
                  help="Peptide identifications.", 
                  name="Peptide IDs",remember=True, 
                  filetypes=[("pepXML","*.pep.xml;*.pep.XML;*.pepXML"),
                             ("mzIdentML","*.mzid;*.mzIdentML"),
                             ("CSV","*.csv"),
			     ("Peptide List","*.txt"),
			     ("X!Tandem","*.tandem.*.gz"),
			     ("k-Score","*.kscore.*.gz"),
			     ("s-Score","*.sscore.*.gz"),
			     ("Mascot","*.mascot.*.gz"),
			     ("OMSSA","*.omssa.*.gz"),
			     ("MyriMatch","*.myrimatch.*.gz"),
			     ("InsPecT","*.inspect.*.gz"),
			     ("MSGF+","*.msgfplus.*.gz"),
			     ("IDPicker Database","*.idpDB")
			    ])
parser.add_option("--psmmetric",type="choice",dest="psmmetric",default=psmmetric_choices[psmmetric_default],
		  help="Metric for filtering peptide identifications. Default: %s."%(psmmetric_choices[psmmetric_default],),
		  choices=psmmetric_choices,remember=True,name="PSM Metric")
parser.add_option("-t","--psmthresh",type="float",dest="psmthresh",default=30,
                  help="Load peptide identifications with PSM metric at or better than this threshold. Default: 30%.",remember=True,
	          name="PSM Threshold")
parser.add_option("-l","--minpeplength",type="int",dest="minpeplen",default=7,
                  help="Load peptides of length at least this threshold. Default: 7.",remember=True,
	          name="Min. Peptide Length")
parser.add_option("-r","--maxrank",type="int",dest="maxrank",default=1,
                  help="Load peptide identifications of rank at most this threshold. Default: 1.",remember=True,
	          name="Max. Rank")
parser.add_option("--ties",action="store_true",dest="ties",default=False,
		  help="Use tied-scores as the same rank",remember=True,
		  name="Accept Ties")
parser.add_option("--noties",action="store_true",dest="noties",default=False,
		  help="Exactly one PSM at each rank",remember=True,
		  name="No Ties")
parser.add_option("--alignments", type="choice", choices=["Load", "None"],
                  default="Load", dest="loadalign",
                  help="Whether to load peptide to protein alignments. Default: Load",
                  name="Peptide Alignments")
parser.add_option("--pracc",type="choice",dest="pracc",default='Auto',remember=True,
                 help="Accession and description extraction rule(s).",
                 name="Accessions",
                 choices=['Auto'] + sorted(AccessionRules))
parser.add_option("--root",type="string",dest="root",default="",
                  help="Analysis names with respect to this root directory . Default: local directory.",remember=True,
	          name="Analysis Root")
parser.add_option("--analysisname",type="string",dest="analysisname",default="",
                  help="Analysis name. Default: Derived from peptide identification filename(s).",remember=True,
	          name="Analysis Name")
parser.add_option("-o","--output",type="savefile",dest="output",default="",remember=True,
	          help="Output file.",
	          name="Output",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-q","--quiet",action="store_true",dest="quiet",default=False,remember=True,
                  name="Quiet",help="Quiet. Default: False.")


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

accrule = None
if opts.loadalign == "Load":
    opts.noprotein = False
    if opts.pracc != 'Auto':
        accrule = AccessionRuleFactory(opts.pracc)
else:
    opts.noprotein = True

if psmmetric_choices.index(opts.psmmetric) in (0,1):
    opts.psmthresh /= 100.0

if not opts.root:
    opts.root = None
else:
    assert os.path.isdir(opts.root)
    opts.root = os.path.abspath(opts.root)

if not opts.analysisname:
    opts.analysisname = None

progress = ProgressText(quiet=opts.quiet)

if len(opts.input) == 0:
    opts.input = args

if opts.ties:
    trank = opts.maxrank
    ntrank = 1e+20
    rank = 1e+20
elif opts.noties:
    trank = 1e+20
    ntrank = opts.maxrank
    rank = 1e+20
else:
    trank = 1e+20
    ntrank = 1e+20
    rank = opts.maxrank

psmfilter = psmmetric_map[opts.psmmetric](opts.psmthresh,maxrank=rank,maxntrank=ntrank,maxtrank=trank,minpeplen=opts.minpeplen)

progress.stage('Initialize PSM database...')
if not opts.output:
    opts.output = opts.database
else:
    psmdb = PSMDb(filename=opts.output,clone=opts.database,quiet=opts.quiet)
    del psmdb

# if False and psmdb.psm_count() > 0:
#     progress.stage('Remove extraneous PSMs...')
#     for a in psmdb.analyses():
#         if (a.getdata('Filter') != psmfilter.__class__.__name__):
#             raise RuntimeError("Analysis filter does not match. Existing PSMs filtered using: %s"%(a.getdata('Filter'),))
#         if (a.getdata('Threshold') < psmfilter.loose):
#             raise RuntimeError("Analysis filter threshold too high. Existing PSMs filtered at: %f"%(a.getdata('Threshold'),))
#     for a in psmdb.analyses():
#         if (a.getdata('Threshold') >= psmfilter.loose):
#             a.setdata('Threshold',psmfilter.loose)
#         if (a.getdata('MaxRank') >= psmfilter.maxrank):
#             a.setdata('MaxRank',psmfilter.maxrank)
#         if (a.getdata('MinPeptideLength') < psmfilter.minpeplen):
#             a.setdata('MinPeptideLength',psmfilter.minpeplen)
#     psmdb.psm_filter(maxvalue=psmfilter.qvalue(psmfilter.loose),
#                      condition=lambda psm: (psm.getdata('rank',1e+20) <= psmfilter.maxrank))
#     psmdb.peptide_filter(condition=lambda pep: (len(pep.sequence) >= psmfilter.minpeplen))

# Load all the PSMs...
progress.stage('Load PSMs from input peptide identifications...',len(opts.input))
psmdb = PSMDb(filename=opts.output,quiet=opts.quiet)
for f in opts.input:
    f = os.path.abspath(f)
    assert(not opts.root or f.startswith(opts.root+os.sep))
    psmdb.load(f,psmfilter,noprotein=opts.noprotein,accrule=accrule,analysisname=opts.analysisname,analysisroot=opts.root,normalize=False)
    progress.update()
progress.done()

if not opts.quiet:
    psmdb.dump(sys.stdout)
