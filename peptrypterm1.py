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
parser.add_option("-o","--output",type="savefile",dest="output",default="",remember=True,
                  help="Output file.",
                  name="Output",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("--annotate",action="store_true",dest="annotate",default=False,remember=True,
                  name="Annotate",help="Annotate. Default: False.")
parser.add_option("--filter",type="choice",dest="filter",
                 default="Specific",remember=True,
                 help="One of N-term, C-term, Semispecific, Specific. Default: Specfic.",
                 name="Condition",choices=["N-term","C-term","Semispecific","Specific"])
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

def sequencecontext(p):
    if p.hasdata('alignments'):
        for al in p.getdata('alignments',[]):
            yield None,al[1],al[4]
    else:
        # assert(False)
        for al in p.alignments:
            lcontext = al.metadata.get('left-context','')
            if not lcontext:
                lcontext = al.metadata.get('left-amino-acid','')
            if not lcontext:
                continue
            rcontext = al.metadata.get('right-context','')
            if not rcontext:
                rcontext = al.metadata.get('right-amino-acid','')
            if not rcontext:
                continue
            yield al,lcontext,rcontext

def annotate_peptide(p):
    progress.update()
    # print p.sequence
    annots = set()
    for al,lcontext,rcontext in sequencecontext(p):
        laas = lcontext[-1] + p.sequence[0];
        raas = p.sequence[-1] + rcontext[0];
        laasmatch = False
        if re.search('^([KR][^P]|-.)$',laas):
            laasmatch = True
        # print laas,laasmatch,
        raasmatch = False
        if re.search('^([KR][^P]|.-)$',raas):
            raasmatch = True
        # print raas,raasmatch
        if laasmatch and raasmatch:
            annots.add('specific')
            if al:
                al.setdata('tryptic-termini','specific')
        elif laasmatch:
            annots.add('n-term')
            if al:
                al.setdata('tryptic-termini','n-term-only')
        elif raasmatch:
            annots.add('c-term')
            if al:
                al.setdata('tryptic-termini','c-term-only')
        else:
            if al:
                al.setdata('tryptic-termini','nonspecific')

    if 'specific' in annots:
        p.setdata('tryptic-termini','specific')
    elif 'n-term' in annots and 'c-term' in annots:
        p.setdata('tryptic-termini','semispecific')
    elif 'n-term' in annots:
        p.setdata('tryptic-termini','n-term')
    elif 'c-term' in annots:
        p.setdata('tryptic-termini','c-term')
    else:
        p.setdata('tryptic-termini','nonspecific')
    return True

def good_peptide(p):
    progress.update()
    # print p.sequence
    for lcontext,rcontext in sequencecontext(p):
        laas = lcontext[-1] + p.sequence[0];
        raas = p.sequence[-1] + rcontext[0];
        laasmatch = False
        if re.search('^([KR][^P]|-.)$',laas):
            laasmatch = True
        # print laas,laasmatch,
        raasmatch = False
        if re.search('^([KR][^P]|.-)$',raas):
            raasmatch = True
        # print raas,raasmatch
        if laasmatch and raasmatch:
            return True
        elif laasmatch and opts.filter in ('N-term','Semispecific'):
            return True
        elif raasmatch and opts.filter in ('C-term','SemiSpecific'):
            return True
    return False
if opts.annotate:
    progress.stage("Annotating peptides...",psmdb.peptide_count())
    psmdb.peptide_filter(condition=annotate_peptide)
else:
    progress.stage("Filtering peptides...",psmdb.peptide_count())
    psmdb.peptide_filter(condition=good_peptide)
progress.done()

if not opts.quiet:
    psmdb.dump(sys.stdout)
