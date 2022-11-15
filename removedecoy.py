#!/usr/bin/env python27

import sys
import re

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

from psmdb.PSMDb import PSMDb, model

parser.add_option("-d","--database",type="file",dest="database",remember=True,
	          help="PSM database file",default="",notNone=True,
	          name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-q","--quiet",action="store_true",dest="quiet",default=False,remember=True,
                  name="Quiet",help="Quiet. Default: False.")
parser.add_option("-o","--output",type="savefile",dest="output",default="",remember=True,
                  help="Output file.",
                  name="Output",
                  filetypes=[("PSM Database","*.psm")])

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

progress.stage('Initialize PSM database...')
if not opts.output:
    psmdb = PSMDb(filename=opts.database,quiet=opts.quiet)
else:
    psmdb = PSMDb(filename=opts.output,clone=opts.database,quiet=opts.quiet)

def notdecoy(it):
   progress.update()
   if hasattr(it,'decoy'):
      return (it.decoy != True)
   return (it.getdata('decoy') != True)

progress.stage('Filter peptides...',psmdb.peptide_count())
psmdb.peptide_filter(condition=notdecoy,normalize=False)
progress.done()
progress.stage('Filter proteins...',psmdb.protein_count())
psmdb.protein_filter(condition=notdecoy,normalize=False)
progress.done()
progress.stage('Filter genes...',psmdb.gene_count())
psmdb.gene_filter(condition=notdecoy,normalize=False)
progress.done()
progress.message('Normalize...')
psmdb.normalize()

if not opts.quiet:
    psmdb.dump(sys.stdout)

