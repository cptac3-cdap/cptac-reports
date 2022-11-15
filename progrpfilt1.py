#!/usr/bin/env python27

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
	          help="PSM database file",default="",
	          name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-g","--groups",type="file",dest="groups",default="",remember=True,
                  help="Group name filtering. Optional.",
                  name="Group names",
                  filetypes=[("All group name files","*.xls;*.csv"),
                             ("CSV","*.csv"),
                             ("Excel","*.xls")])
parser.add_option("-t","--type",type="str",dest="type",default="",remember=True,
                  help="Group type. Default: All types.",
                  name="Group Type")
parser.add_option("-v","--complement",action="store_true",dest="complement",default=False,
                  remember=True,name="Complement",help="Complement filter list.")
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

progress = ProgressText(quiet=opts.quiet)

progress.stage('Initialize PSM database...')
if not opts.output:
    psmdb = PSMDb(filename=opts.database,quiet=opts.quiet)
else:
    psmdb = PSMDb(filename=opts.output,clone=opts.database,quiet=opts.quiet)

if opts.groups:
    progress.stage("Read group names...")
    if opts.groups.endswith(".xls"):
        t = XLSFileTable(opts.groups)
    else:
        t = CSVFileTable(opts.groups)
    groups = set()
    for r in t.rows():
        groups.add(r[t.headers()[0]])
    progress.done()
else:
    groups = None

if opts.type:
    type = opts.type
else:
    type = None

progress.message("Filter groups...")

if not opts.complement:
    psmdb.proteingroup_filter(keep=groups,type=type)
else:
    assert(groups)
    psmdb.proteingroup_filter(discard=groups,type=type)

progress.message("Remove orphan proteins...")

psmdb.protein_filter(condition=lambda pr: len(pr.groups()) > 0)

if not opts.quiet:
    psmdb.dump(sys.stdout)
