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
parser.add_option("-g","--genemap",type="file",dest="genemap",default="",remember=True,
                  help="Protein accession to gene mapping. Optional.",
                  name="Protein to Gene",
                  filetypes=[("All protein to gene files","*.xls;*.csv"),
                             ("CSV","*.csv"),
                             ("Excel","*.xls")])
parser.add_option("-R","--remove",action="store_true",dest="remove",default=False,remember=True,                  
                  name="Remove unassociated?",help="Remove proteins not associated with a gene? Default: False.") 
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

if opts.genemap:
    progress.stage("Read accession to gene mapping file...")
    if opts.genemap.endswith(".xls"):
        t = XLSFileTable(opts.genemap)
    else:
        t = CSVFileTable(opts.genemap)
    # expect headers: <free-choice!=accession>, accession, description, url
    assert(t.headers()[0] not in t.headers()[1:])
    assert('accession' in t.headers()[1:])
    assert(len(t.headers()) > 1)
    genemap = defaultdict(list)
    for r in t.rows():
        key = r[t.headers()[0]]
        if not opts.decoyprefix:
            genemap[key].append(r)
        else:
            r1 = dict(r.items())
            r1['decoy'] = False
            genemap[key].append(r1)
            dkey = opts.decoyprefix+key
            dacc = opts.decoyprefix+r1['accession']
            genemap[dkey].append(dict(accession=dacc,decoy=True))
    progress.done()
    progress.message("Map protein accessions to gene groups...")
    psmdb.gene_map(genemap,opts.remove)

if not opts.quiet:
    psmdb.dump(sys.stdout)
