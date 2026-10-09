#!/usr/bin/env python27

import sys, os, os.path, copy, re, math, gzip, bz2
from collections import defaultdict
from operator import itemgetter
from functools import partial
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

dam = psmdb.modification('N',+0.984)

ngly = re.compile(r'N[^P][ST]')
progress.stage("Filtering out peptides that don't have N-glycan motif...",psmdb.ion_count())
def good_ion(pi):
    progress.update()
    p = pi.peptide
    for al in p.alignments:
        rcontext = al.metadata.get('right-context','')
        if not rcontext:
            rcontext = al.metadata.get('right-amino-acid','')
        assert(rcontext)
        sites = set()
        count = 0
        for m in ngly.finditer(p.sequence+rcontext):
            if m.start() < len(p.sequence):
                count += 1
                sites.add(m.start()+1)
        if count >= 1:
            for ms in pi.modifiedSites:
                if ms.modification == dam:
                    if ms.position in sites:
                        return True
            return False
    return False
psmdb.ion_filter(condition=good_ion)
progress.done()

if not opts.quiet:
    psmdb.dump(sys.stdout)
