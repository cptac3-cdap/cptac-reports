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
parser.add_option("-p","--protmap",type="file",dest="protmap",default="",remember=True,
                  help="Protein accession mapping/filtering. Optional.",
                  name="Protein Mapping",
                  filetypes=[("All protein mapping files","*.xls;*.csv"),
                             ("CSV","*.csv"),
                             ("Excel","*.xls")])
parser.add_option("-v","--complement",action="store_true",dest="complement",default=False,
                  remember=True,name="Complement",help="Complement filter list.")
parser.add_option("--decoyprefix",type="str",dest="decoyprefix",default=None,remember=True,
                  name="Decoy Prefix",help="Decoy prefix. Default: None.")
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

if opts.protmap:
    progress.stage("Read accession to protein mapping file...")
    if opts.protmap.endswith(".xls"):
        t = XLSFileTable(opts.protmap)
    else:
        t = CSVFileTable(opts.protmap)
    # expect headers: <free-choice!=accession>, accession, description, url, length, defline
    assert(t.headers()[0] not in t.headers()[1:])
    assert('accession' in t.headers()[1:] or len(t.headers())==1)
    protmap = defaultdict(list)
    for r in t.rows():
        protmap[r[t.headers()[0]]].append(r)
        if opts.decoyprefix:
            key = opts.decoyprefix+r[t.headers()[0]]
            if len(t.headers())>1:
                accession = opts.decoyprefix+r.get('accession')
                protmap[key] = dict(accession=accession)
            else:
                protmap[key] = dict()
    progress.done()
    if len(t.headers()) == 1:
        progress.message("Filter proteins...")
        if opts.complement:
            psmdb.protein_filter(discard=protmap.keys())
        else:
            psmdb.protein_filter(keep=protmap.keys())
    else:
        progress.message("Map protein accessions and metadata...")
        psmdb.protein_map(protmap)
        trans = psmdb.begin()
        for pr in psmdb.proteins(connection=trans):
            if pr.accession.startswith(opts.decoyprefix):
                pr.decoy = True
            else:
                pr.decoy = False
        trans.commit(close=True)

if not opts.quiet:
    psmdb.dump(sys.stdout)
