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
parser.add_option("--orthmap",type="file",dest="orthmap",default="",remember=True,
                  help="Protein/Gene accession ortholog map. Optional.",
                  name="Protein Orthologs",
                  filetypes=[("All protein mapping files","*.xls;*.csv"),
                             ("CSV","*.csv"),
                             ("Excel","*.xls")])
parser.add_option("--bygene",action="store_true",dest="bygene",default=False,remember=True,
                  name="By Gene",help="Gene-based orthologs. Default: False.")
parser.add_option("--decoyprefix",type="str",dest="decoyprefix",default=None,remember=True,
                  name="Decoy Prefix",help="Decoy prefix. Default: None.")
parser.add_option("--undo",action="store_true",dest="undo",default=False,remember=True,
                  name="Undo",help="Undo ortholog pairs. Default: False.")
parser.add_option("-o","--output",type="savefile",
                  dest="output",default="",remember=True,
                  help="Output file.",
                  name="Output",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-q","--quiet",action="store_true",
                  dest="quiet",default=False,remember=True,
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

if opts.orthmap:
    progress.stage("Read ortholog mapping file...")
    if opts.orthmap.endswith(".xls"):
        t = XLSFileTable(opts.orthmap)
    else:
        t = CSVFileTable(opts.orthmap)
    h1 = t.headers()[0]
    h2 = t.headers()[1]
    seenh2 = set()
    orthmap = dict()
    for r in t.rows():
        assert r[h1] not in orthmap, "Duplicate accession %s"%(r[h1],)
        assert r[h2] not in seenh2, "Duplicate accession %s"%(r[h2],)
        orthmap[r[h1]] = r[h2]
        if opts.decoyprefix:
            orthmap[opts.decoyprefix+r[h1]] = opts.decoyprefix+r[h2]
        seenh2.add(r[h2])
    progress.done()
    progress.message("Construct ortholog pairs...")
    if not opts.bygene:
        psmdb.protein_orthologs(orthmap)
    else:
        psmdb.gene_orthologs(orthmap)

if opts.undo:
    progress.message("Undoing ortholog pairs...")
    if not opts.bygene:
        psmdb.undo_protein_orthologs()
    else:
        psmdb.undo_gene_orthologs()

if not opts.quiet:
    psmdb.dump(sys.stdout)
