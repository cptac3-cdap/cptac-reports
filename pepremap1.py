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


from psmdb.PSMDb import PSMDb

parser.add_option("-d","--database",type="file",dest="database",remember=True,
	          help="PSM database file",default="",notNone=True,
	          name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-s","--seqdb",type="file",dest="seqdb",default="",remember=True,
                 help="Protein sequence database for remapping proteins. Optional.",
                 name="Sequence Database",
                 filetypes=[("FASTA sequence database","*.fasta")])
parser.add_option("--pracc",type="choice",dest="pracc",default='Auto',remember=True,
                 help="Accession and description extraction rule(s).",
                 name="Accessions",
                 choices=['Auto'] + sorted(AccessionRules))
parser.add_option("-o","--output",type="savefile",dest="output",default="",remember=True,
	          help="Output file.",
	          name="Output",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-U","--update",action="store_true",dest="update",default=False,remember=True,
                  help="Update (add) protein sequences. Default: Replace protein sequences.",
                  name="Update")
parser.add_option("-O","--orphans",action="store_true",dest="orphan",default=False,remember=True,
                  help="Retain orphan peptide sequences. Default: Remove orphan peptides. Default: False.",
                  name="Orphans")
parser.add_option("-D",type="choice",dest="decoy",default='None',remember=True,
                 help="Mark proteins (and their aligned peptides) decoy status. Default: Do not mark decoy status.", 
                 name="Decoy", choices=['None','True','False'])
parser.add_option("--dnamut",action="store_true",dest="dnamut",default=False,remember=True,
                  help="Consider single DNA mutation substitutions in mapped peptides. Default: False.",
                  name="DNA Mutations")
parser.add_option("--aasubst",action="store_true",dest="aasubst",default=False,remember=True,
                  help="Consider single amino-acid substitutions in mapped peptides. Default: False.",
                  name="AA Subst.")
parser.add_option("--isosubst",action="store_true",dest="isosubst",default=False,remember=True,
                  help="Consider isobaric amino-acid substitutions in mapped peptides. Default: False.",
                  name="Isobaric Subst.")
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

opts.decoy = eval(opts.decoy)

progress = ProgressText(quiet=opts.quiet)

accrule = None
if opts.seqdb or (opts.pracc != 'Auto'):
    if opts.pracc == 'Auto':
        opts.pracc = AccessionRuleAuto(opts.seqdb)
    else:
        opts.pracc = AccessionRuleFactory(opts.pracc)
    accrule = opts.pracc

progress.stage('Initialize PSM database...')
if not opts.output:
    psmdb = PSMDb(filename=opts.database,quiet=opts.quiet)
else:
    psmdb = PSMDb(filename=opts.output,clone=opts.database,quiet=opts.quiet)

# Map the peptides to a new database as necessary
if opts.seqdb:
    progress.stage("Remapping peptides to %s"%(opts.seqdb))
    if opts.dnamut:
        psmdb.peptide_remap(opts.seqdb,accrule,dnamut=1,update=opts.update,decoy=opts.decoy,noprotein=opts.orphan)
    elif opts.aasubst:
        psmdb.peptide_remap(opts.seqdb,accrule,aasubst=1,update=opts.update,decoy=opts.decoy,noprotein=opts.orphan)
    elif opts.isosubst:
        psmdb.peptide_remap(opts.seqdb,accrule,aasubst=1,validdeltas=(0.05,),update=opts.update,decoy=opts.decoy,noprotein=opts.orphan)
    else:
        psmdb.peptide_remap(opts.seqdb,accrule,update=opts.update,decoy=opts.decoy,noprotein=opts.orphan)

if not opts.quiet:
    psmdb.dump(sys.stdout)
