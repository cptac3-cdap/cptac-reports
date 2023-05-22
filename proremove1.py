#!/usr/bin/env python3

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
parser.add_option("-f","--field",type="string",dest="field",default="defline",remember=True,
                  help="Protein field to match the regular expression to.",
                  name="Protein Regex")
parser.add_option("-r","--regex",type="string",dest="protregex",default="",remember=True,
                  help="Select proteins to keep by regular expression matches to protein field.",
                  name="Reg. Exp.")
parser.add_option("-v","--nonmatch",action="store_true",dest="nonmatch",default=False,remember=True,
                  name="Non-Match",help="Keep non-matching proteins. Default: False.")
parser.add_option("-M","--missing",action="store_true",dest="missing",default=False,remember=True,
                  name="Missing",help="Consider proteins missing the field to match the regular expression. Default: False.")
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

if opts.protregex:
    opts.protregex = re.compile(opts.protregex)
    nprot = psmdb.protein_count()
    progress.stage("Filtering proteins...",nprot)
    keep = True
    if opts.nonmatch:
        keep = False
    def good_protein(pr):
        progress.update()
        if (opts.field != "accession") and (not pr.metadata or opts.field not in pr.metadata):
            if opts.missing:
                return keep
            return (not keep)
        if opts.field != "accession":
            desc = pr.metadata.get(opts.field)
        else:
            desc = pr.accession
        if (not opts.protregex) or opts.protregex.search(desc):
            return keep
        return (not keep)
    psmdb.protein_filter(condition=good_protein)
    progress.done()

if not opts.quiet:
    psmdb.dump(sys.stdout)
