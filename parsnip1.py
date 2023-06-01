#!/usr/bin/env python3

from dataset import CSVFileTable, XLSFileTable, XLSFileDataset, ZIPFileDataset, ParsimonyCSV
import sys, os, os.path, copy, re, math, gzip, bz2, random, time
from collections import defaultdict
from operator import itemgetter, attrgetter
from functools import partial, cmp_to_key
from psmdb.PepXML import PepXMLReader, ParsimonyPepXML
from peptidescan.PeptideRemapper import *
from peptidescan.OutOfCoreTable import *

import atexit
from functools import reduce
def onexit():
    print("Type <Enter> to Exit...", end=' ', file=sys.stderr)
    sys.stderr.flush()
    input()

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

from psmdb.filter import SpectralFDR,PeptideFDR,PeptideProphet,EValue
from psmdb.PSMDb import PSMDb

psmmetric_choices = ["Estimated PSM FDR % (PepArML)",
                     "Estimated Peptide FDR % (PepArML)",
                     "PSM Probability (Peptide Prophet)",
                     "PSM E-value (PepArML)"]
psmmetric_default = 0
psmmetric_map = dict(list(zip(psmmetric_choices,
                         [SpectralFDR,PeptideFDR,PeptideProphet,EValue])))

parsimony = OptionGroup(parser, "Parsimony")
output = OptionGroup(parser, "Output")
parser.add_option("-d","--database",type="file",dest="database",remember=True,
                  help="PSM database file",default="",notNone=True,
                  name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("--psmmetric",type="choice",dest="psmmetric",default=psmmetric_choices[psmmetric_default],
                  help="Metric for filtering peptide identifications. Default: %s."%(psmmetric_choices[psmmetric_default],),
                  choices=psmmetric_choices,remember=True,name="PSM Metric")
parser.add_option("-t","--psmthresh",type="float",dest="psmthresh",default=None,
                  help="Use peptide identifications with PSM metric at or better than this threshold. Default: No constraint.",remember=True,
                  name="PSM Threshold")

parser.add_option("-T","--maxfdr",type="float",dest="maxfdr",default=None,
                  help="Consider only peptides with minfdr at most this threshold. Default: No constraint.",remember=True,
                  name="Max. Pep. FDR")

parser.add_option("-C","--minspeccount",type="int",dest="minspeccount",default=1,
                  help="Consider only peptides with spectral counts at least this threshold. Default: No constraint.",remember=True,
                  name="Min. Spec. Count")

parsimony.add_option("--parsimony",type="choice",dest="parsimonychoice",default='Parsimony',remember=True,
                  help="Parsimony strategy for eliminating redundant proteins. One of Parsimony, FixedPoint, None. Default: Parsimony",
                  choices=["Parsimony","FixedPoint","None"],
                  name="Parsimony")
parsimony.add_option("--tieresolution",type="choice",dest="tieresolution",
                     default='Discard',remember=True,
                     help="Strategy for dealing with spectra with multi-peptide PSMs. One of discard multi-peptide spectra (Discard), keep only multi-peptide consitent spectra (Consistent), keep all multi-peptide spectra (Keep), group spectra by associated peptides (Group). Default: Discard",
                     choices=["Discard","Consistent","Keep","Group"],
                     name="Multi-Peptide Spectra:")
parsimony.add_option("--bygene",action="store_true",default=False,remember=True,
                     name="By Gene",help="Gene-based parsimony. Default: False.")
parsimony.add_option("--pracc",type="choice",dest="pracc",default='Auto',remember=True,
                 help="Accession and description preference rule(s).",
                 name="Accessions",
                 choices=['Auto'] + sorted(AccessionRules))
parsimony.add_option("-U","--unshared",type="int",dest="unique",default=1,remember=True,
                  help="Proteins must have at least this many unshared peptides in the final solution. Default: 1",
                  name="Unshared")
parsimony.add_option("--pepweights",type="choice",dest="pepweights",default="Equal",remember=True,name="Peptide Weights",help="Weight peptides equally, by specral count, or by -log(FDR). One of Equal, Count, FDR. Default: Equal.", choices=["Equal","Count","FDR"])
parsimony.add_option("--ncompprot",type="str",dest="ncompprot",default="",remember=True,
                     name="Max. Proteins / Comp.",help="Max. # of proteins per component. Default: Unlimited.")
parsimony.add_option("--noalignments",action="store_true",dest="noalignments",default=False,remember=True,
                  name="No Alignments",help="Do not order peptides by protein alignment in Matrix mode output. Default: False.")
parsimony.add_option("--itersolve",type="int",dest="itersolve",default=1e+20,remember=True,
                  name="Iter. Solve Thr.",help="Heuristic iterative solve threshold. Default: 20 proteins.")
parsimony.add_option("--tinitial",type="int",dest="tinitial",default=1e+20,remember=True,
                  help="Time to look for an initial solution in minutes. Default: No limit.",
                  name="Initial Sol. Time")
parsimony.add_option("--timprove",type="int",dest="timprove",default=1e+20,remember=True,
                  help="Time to look for an improved solution in minutes. Default: 10",
                  name="Improved Sol. Time")
parsimony.add_option("--ttotal",type="int",dest="ttotal",default=1e+20,remember=True,
                  help="Total time find a solution in minutes. Default: No limit.",
                  name="Improved Sol. Time")
parser.add_option("--mode", type="multichoice", dest="mode", default=["Stats"], remember=True, name="Mode", help="Output mode - dump parismony selected protein accessions, write protein-peptide matrix for each component, filter out unselected proteins, or annotate selected proteins. One of Dump, Matrix, Filter, Annotate, Stats.", multichoices=["Dump","Matrix","Filter","Annotate","Stats"])
parser.add_option("-q","--quiet",action="store_true",dest="quiet",default=False,remember=True,
                  name="Quiet",help="Quiet. Default: False.")
parser.add_option("-E","--extension",type="str",dest="extn",default=None,
                  help="Extension for automatically named files. Default: \"pars\".",remember=True,
                  name="Extension")
parser.add_option("-o","--output",type="savefile",dest="output",default="",remember=True,
                  help="Output file.",
                  name="Output",
                  filetypes=[("PSM Database","*.psm"),
                             ("Accessions","*.txt"),
                             ("Matrix","*.csv")])

parser.add_option_group(parsimony)

opts = None
while True:
    if 'exit' in error_kwargs:
        try:
            opts,args = parser.parse_args(opts=opts)
        except UserCancelledError:
            os._exit(0);
    else:
        opts,args = parser.parse_args()

    if isinstance(opts.mode,str):
        opts.mode = list(map(str.strip,opts.mode.split(',')))

    if len(opts.mode) > 1:
        if opts.output:
            parser.error('Output file must be implicit if more than one mode is chosen...',**error_kwargs)

    if opts.ncompprot:
        try:
            opts.ncompprot = int(opts.ncompprot)
            if opts.ncompprot <= 0:
                parser.error('Parismony:Max. # of Proteins per Component (--ncompprot) must be greater than zero.',**error_kwargs)
                continue
        except:
            parser.error('Parismony:Max. # of Proteins per Component (--ncompprot) must be an integer.',**error_kwargs)
            continue

    # various tests for argument consistency...

    break


psmmetricpercent = False
if opts.psmthresh != None and psmmetric_choices.index(opts.psmmetric) in (0,1):
    opts.psmthresh /= 100.0
    psmmetricpercent = True
if opts.maxfdr != None and psmmetric_choices.index(opts.psmmetric) in (0,1):
    opts.maxfdr /= 100.0
opts.parsimony = False
opts.fixedpoint = False
if opts.parsimonychoice == 'Parsimony':
    opts.parsimony = True
elif opts.parsimonychoice == 'FixedPoint':
    opts.fixedpoint = True
else:
    pass
opts.specweights = False
opts.fdrweights = False
if opts.pepweights == 'Count':
    opts.specweights = True
elif opts.pepweights == 'FDR':
    opts.fdrweights = True

opts.pdump = False
opts.pmatrix = False
opts.pfilter = False
opts.pannotate = False
opts.pstats = False
if 'Dump' in opts.mode:
    opts.pdump = True
if 'Matrix' in opts.mode:
    opts.pmatrix = True
if 'Filter' in opts.mode:
    opts.pfilter = True
if 'Stats' in opts.mode:
    opts.pstats = True
if 'Annotate' in opts.mode:
    opts.pannotate = True

opts.discardmultipep = False
opts.discardambigmultipep = False
opts.peptidegroups = False
if opts.tieresolution == 'Discard':
    opts.discardmultipep = True
elif opts.tieresolution == 'Keep':
    pass
elif opts.tieresolution == 'Group':
    opts.peptidegroups = True
elif opts.tieresolution == 'Consistent':
    opts.peptidegroups = True
    opts.discardambigmultipep = True

if not opts.extn:
    opts.extn = "pars"

# This takes care of text-based outputs...
if opts.output and not opts.pfilter and not opts.pannotate:
    opts.outhandle = open(opts.output,'w')
else:
    opts.outhandle = sys.stdout

progress = ProgressText(quiet=opts.quiet)

accrule = None
if opts.pracc != 'Auto':
    accrule = AccessionRuleFactory(opts.pracc)

psmfilter = psmmetric_map[opts.psmmetric](opts.psmthresh if opts.psmthresh != None else 1e+20)


# Check the input database has been pre-filtered appropraitely...
progress.message("Initialize PSM database...")
tbpsmdb = PSMDb(opts.database,quiet=opts.quiet)

needtb=False
needrmtie=False
needrmtie2=False
needpg=False
needclone=False
clonefn=None

if opts.psmthresh != None:
    for a in tbpsmdb.analyses():
        if (a.getdata('Filter') != psmfilter.__class__.__name__):
            print("PSM filter does not match. Input PSMs filtered using: %s"%(a.getdata('Filter'),), file=sys.stderr)
            sys.exit(1)
        elif (a.getdata('Threshold') > psmfilter.tight):
            needtb = True
            needclone = True
        elif (a.getdata('Threshold') < psmfilter.tight):
            print("PSMs already filtered more stringently than requested. Input PSMs filtered at: %s"%(a.getdata('Threshold'),), file=sys.stderr)
            sys.exit(1)

if opts.peptidegroups:
    needpg = True
    needclone = True

if opts.discardmultipep:
    needrmtie = True
    needclone = True

if opts.discardambigmultipep:
    needrmtie2 = True
    needclone = True

if opts.pannotate or opts.pfilter:
    needclone = True
    if opts.output:
        clonefn = opts.output
    else:
        clonefn = opts.database[:-4]+'.%s.psm'%opts.extn

pgtype = None
if needclone:
    progress.stage("Clone PSM database...")
    if clonefn:
        psmdb = tbpsmdb.clone(filename=clonefn,quiet=opts.quiet)
    else:
        psmdb = tbpsmdb.clone(quiet=opts.quiet)
    progress.done()
    if needtb:
        progress.stage("Remove PSMs above threshold from PSM database...")
        psmdb.psm_filter(maxvalue=psmfilter.qvalue(opts.psmthresh))
        progress.done()
    if needpg or needrmtie or needrmtie2:
        progress.stage("Making peptide groups by spectra...")
        psmdb.peptideGroupsBySpectra()
        pgtype = 'SpectraEquivClass'
        progress.done()
    if needrmtie:
        progress.stage("Remove multi-peptide spectra from PSM database...")
        psmdb.removeMultiPeptideSpectra()
        progress.done()
    elif needrmtie2:
        progress.stage("Keep only multi-peptide consistent spectra from PSM database...")
        psmdb.removeAmbiguousMultiPeptideSpectra()
        psmdb.peptideGroupsBySpectra()
        progress.done()
else:
    psmdb = tbpsmdb

progress.stage("Extract protein/peptide relationships and statistics...")
progress.stage("Peptides")
peptides = set(psmdb.peptideids())
progress.done()

if not opts.bygene:
    progress.stage("Proteins")
    proteins = set(psmdb.proteinids())
    progress.done()
else:
    progress.stage("Genes")
    proteins = set(psmdb.geneids())
    progress.done()

if opts.peptidegroups:
    progress.stage("PeptideGroups")
    pepgroups = set(psmdb.peptidegroupids(pgtype))
    progress.done()

if not opts.bygene:

    progress.stage("Proteins2Peptides")
    pr2pep = dict(psmdb.proteinid2peptideids())
    progress.done()
    progress.stage("Peptides2Proteins")
    pep2pr = dict(psmdb.peptideid2proteinids())
    progress.done()

    if opts.peptidegroups:
        progress.stage("Proteins2PeptideGroups")
        pr2pepgrp = psmdb.proteinid2peptidegroupids(pgtype)
        progress.done()
        progress.stage("PeptideGroups2Proteins")
        pepgrp2pr = psmdb.peptidegroupid2proteinids(pgtype)
        progress.done()
        progress.stage("PeptideGroups2Peptides")
        pepgrp2pep = psmdb.peptidegroupid2peptideids(pgtype)
        progress.done()
        progress.stage("PeptideGroupByProteins2Peptides")
        pepgrpprot2pep = psmdb.peptidegroupidproteinid2peptideids(pgtype)
        progress.done()

else:

    progress.stage("Genes2Peptides")
    pr2pep = dict(psmdb.geneid2peptideids())
    progress.done()
    progress.stage("Peptides2Genes")
    pep2pr = dict(psmdb.peptideid2geneids())
    progress.done()

    if opts.peptidegroups:

        progress.stage("Genes2PeptideGroups")
        pr2pepgrp = psmdb.geneid2peptidegroupids(pgtype)
        progress.done()
        progress.stage("PeptideGroups2Genes")
        pepgrp2pr = psmdb.peptidegroupid2geneids(pgtype)
        progress.done()
        progress.stage("PeptideGroups2Peptides")
        pepgrp2pep = psmdb.peptidegroupid2peptideids(pgtype)
        progress.done()
        progress.stage("PeptideGroupByGenes2Peptides")
        pepgrpprot2pep = psmdb.peptidegroupidgeneid2peptideids(pgtype)
        progress.done()


progress.stage("Peptides2PSMCount")
pep2psmcount =  psmdb.peptideid_psmcounts()
progress.done()
progress.stage("Peptides2FDR")
pep2fdr = psmdb.peptideid_psmminvalue()
progress.done()
if opts.peptidegroups:
    progress.stage("PeptideGroups2PSMCount")
    pepgrp2psmcount =  psmdb.peptidegroupid_psmcounts(pgtype)
    progress.done()
    progress.stage("PeptideGroups2FDR")
    pepgrp2fdr = psmdb.peptidegroupid_psmminvalue(pgtype)
    progress.done()

progress.stage("Min FDR value")
try:
    minqvalue = min([v for v in iter(pep2fdr.values()) if v>0])
except ValueError:
    minqvalue = 0
progress.done()

if opts.peptidegroups:
    try:
        minpepgrpqvalue = min([v for v in iter(pepgrp2fdr.values()) if v>0])
    except ValueError:
        minpepgrpqvalue = 0


prpep2align = dict()
if opts.pmatrix and not opts.noalignments and not opts.bygene:
    progress.stage("Alignments")
    prpep2align = dict(psmdb.alignmentids2pos())
    if opts.peptidegroups:
        prpepgrp2align = dict()
        for (pepgrp,prot),peps in pepgrpprot2pep.items():
            thepep = next(iter(peps))
            prpepgrp2align[(prot,pepgrp)] = prpep2align.get((prot,thepep),(-1,-1))
        prpep2align = prpepgrp2align
    progress.done()

# progress.stage("Substitutions")
# peppr2subst = dict(map(lambda t: (t[0],"%s%d->%s:%+.2f"%t[1]),psmdb.alignmentids2subst()))
peppr2subst = {}

if needtb:
    assert(False)
    progress.stage("Tie-breaker Proteins2Peptides")
    pr2tbpep = tbpsmdb.proteinid2peptideids()
    progress.done()
    progress.stage("Tie-breaker Peptides2FDR")
    tbpep2fdr = tbpsmdb.peptideid_psmminvalue()
    progress.done()
else:
    pr2tbpep = defaultdict(set)
    tbpep2fdr = defaultdict(float)

progress.stage("Counts...")
nspec = psmdb.spectrum_count()
npsm = psmdb.psm_count()
npep = psmdb.peptide_count()
nprot = psmdb.protein_count()
ngene = psmdb.gene_count()
npepgrp = psmdb.peptidegroup_count(pgtype)
progress.done()

## print >>sys.stderr, "PSMDb:"
## psmdb.dump()
## print >>sys.stderr, "Tie Breaker PSMDb:"
## tbpsmdb.dump()

progress.stage('Spectra:   %s'%nspec)
progress.stage('PSMs:      %s'%npsm)
progress.stage('Peptides:  %s'%npep)
progress.stage('PepGroups: %s'%npepgrp)
progress.stage('Proteins:  %s'%nprot)
if opts.bygene:
    progress.stage('Genes:     %s'%ngene)

from parsimony import Parsimony, Components, MinUncovered, Dominator, GraphDumper, FixedPoint

if opts.peptidegroups:
    edges = pr2pepgrp
    elements = pepgroups
    psmcount = pepgrp2psmcount
    psmfdr   = pepgrp2fdr
    minqvalue = minpepgrpqvalue
else:
    edges = pr2pep
    elements = peptides
    psmcount = pep2psmcount
    psmfdr   = pep2fdr

if opts.minspeccount > 1:
    st = len(elements)
    elements.difference_update(set(map(itemgetter(0),[t for t in iter(psmcount.items()) if t[1] < opts.minspeccount])))
    rmed = st-len(elements)
    progress.message("Removed %d peptides with spectral count < %d"%(rmed,opts.minspeccount))

if opts.maxfdr != None:
    st = len(elements)
    elements.difference_update(set(map(itemgetter(0),[t for t in iter(psmfdr.items()) if t[1] > opts.maxfdr])))
    rmed = st-len(elements)
    progress.message("Removed %d peptides with min FDR > %f"%(rmed,opts.maxfdr))

def tbcmp(pri,prj):
    tbi = sorted(map(tbpep2fdr.get,pr2tbpep[pri]-pr2tbpep[prj]),key=cmp_to_key(psmfilter.cmp))
    tbj = sorted(map(tbpep2fdr.get,pr2tbpep[prj]-pr2tbpep[pri]),key=cmp_to_key(psmfilter.cmp))
    tbi.extend([1e+20]*max(len(tbj)-len(tbi),0))
    tbj.extend([1e+20]*max(len(tbi)-len(tbj),0))
    for c in [psmfilter.cmp(*t) for t in zip(tbi,tbj)]:
        if c != 0:
            return c
    return 0

pr = None
if not opts.bygene:
    for pepid,cnt in sorted(iter(pep2psmcount.items()),key=itemgetter(1),reverse=True):
        prs = list(psmdb.peptide(pepid).proteins())
        if len(prs) > 0:
            pr = prs[0]
            break
else:
    for pepid,cnt in sorted(iter(pep2psmcount.items()),key=itemgetter(1),reverse=True):
        gns = list(psmdb.peptide(pepid).genes())
        if len(gns) > 0:
            pr = gns[0]
            break

def addaccdesc(pri):
    acc = psmdb.proteinid_accession(pri)
    desc = psmdb.proteinid_description(pri)
    if re.search(r'\b%s\b'%re.escape(acc),desc.split()[0]):
        return desc
    return " ".join([acc,desc])

def addnamedesc(gni):
    acc = psmdb.geneid_name(gni)
    desc = psmdb.geneid_description(gni)
    if desc and re.search(r'\b%s\b'%re.escape(acc),desc.split()[0]):
        return desc
    return " ".join([acc,desc])

if opts.bygene:
    addaccdesc = addnamedesc

praccfn = psmdb.proteinid_accession
if opts.bygene:
    praccfn = psmdb.geneid_name

accfn = None
if pr.metadata and pr.metadata.get('defline'):
    if not accrule:
        accrule = TestAccessionRules(pr.metadata.get('defline'))
    accfn = psmdb.proteinid_defline
    # print "ACCFN: proteinid_defline"
if (not accrule or not accfn) and pr.metadata and pr.metadata.get('description'):
    if not accrule:
        accrule = TestAccessionRules(pr.metadata.get('description'))
    accfn = addaccdesc
    # print "ACCFN: proteinid_description"
if not accrule:
    accrule = AccessionRuleFactory('FirstWord')
    if pr.metadata and pr.metadata.get('defline'):
        accrule = psmdb.proteinid_defline
        # print "ACCFN: proteinid_defline"
    elif pr.metadata and pr.metadata.get('description'):
        accrule = addaccdesc
        # print "ACCFN: proteinid_description"
    else:
        accfn = (lambda x: None)
        # print "ACCFN: None"
# print "ACCRULE:",str(accrule)

from functools import cmp_to_key

prsortkey = dict(list(zip(list(map(itemgetter(0),
                         sorted([(pri,accfn(pri)) for pri in proteins],
                                key=cmp_to_key(accrule.prefer)))),
                     list(range(len(proteins))))))
prsortkey1 = dict(list(zip(sorted(proteins),list(range(len(proteins))))))

# proteinlist = sorted(proteins,key=prsortkey.get)
# for i in range(len(proteinlist)):
#    pri = proteinlist[i]
#    for j in range(i+1,len(proteinlist)):
#        prj = proteinlist[j]
#       if accrule.prefer((pri,accfn(pri)),(prj,accfn(prj))) > 0:
#            print pri,prsortkey.get(pri),accfn(pri)
#            print prj,prsortkey.get(prj),accfn(prj)

fp = None

if opts.parsimony:

    # print len(elements),len(proteins),len(edges)

    dom = Dominator(peptides=elements,edges=edges,proteins=proteins)
    dom.dominate(tiebreakers=(pr2tbpep,tbpep2fdr,psmfilter.cmp),prsortkey=prsortkey.get,progress=progress)

    parsimony = set()

    weights = None
    if opts.specweights:
        weights = psmcount
    elif opts.fdrweights:
        weights = {}
        for k,v in list(psmfdr.items()):
            if v > 0.0:
                weights[k] = -math.log(v)
            else:
                weights[k] = -math.log(minqvalue/2)

    uitems = None
    if opts.peptidegroups:
        uitems = pepgrpprot2pep

    comp = Components(peptides=dom.peptides,edges=dom.edges,proteins=dom.proteins,progress=progress)

    for i,(pr,pep) in enumerate(comp):

        # if len(pr) > 10:
        #     gd = GraphDumper("Component%d"%i,dom.proteins,dom.peptides,dom.edges,pr,pep,
        #                      weights,psmdb.proteinid_accession,pepseq)

        if opts.peptidegroups and uitems:
            maxncpep = 0
            for pri in pr:
                pepset = reduce(set.union,[uitems.get((pgi,pri),set()) for pgi in pep],set())
                maxncpep = max(maxncpep,len(pepset))
            ncpep = maxncpep
        else:
            ncpep = len(pep)

        if ncpep >= opts.unique:

            if len(pr) == 1:

                parsimony |= pr

            else:

                kw = {}

                if opts.ncompprot:
                    kw['nprot'] = opts.ncompprot

                # for i1,pri in enumerate(sorted(pr,key=prsortkey.get)):
                #     print i1,pri,prsortkey.get(pri),accfn(pri)

                pars = MinUncovered(peptides=pep,edges=dom.edges,proteins=pr,
                                    minpepper=opts.unique,weights=weights,uitems=uitems,
                                    prsortkey=prsortkey.get,quiet=opts.quiet,**kw)

                if len(pr) >= opts.itersolve:

                    upep = defaultdict(set)
                    for pri in pr:
                        for pepi in dom.edges[pri]:
                            if opts.peptidegroups:
                                for pepi in pepgrp2pep[pepi]:
                                    upep[pepi].add(pri)
                            else:
                                upep[pepi].add(pri)
                    upepcnt=defaultdict(int)
                    for pri in pr:
                        upepcnti = 0
                        for pepi in dom.edges[pri]:
                            if opts.peptidegroups and uitems:
                                good = True
                                for pepi in uitems[(pepi,pri)]:
                                    if len(upep[pepi]) > 1:
                                        good = False
                                        break
                                if good:
                                    upepcnti += 1
                            else:
                                if len(upep[pepi]) == 1:
                                    upepcnti += 1
                        upepcnt[pri] = upepcnti

                    solution0 = ((1e+20,1e+20),None)
                    first = True
                    uthresh = opts.unique
                    nforcein = 1e+20
                    ttotal = opts.ttotal*60
                    while True:
                        forcein = list(map(itemgetter(0),[t for t in list(upepcnt.items()) if t[1] >= uthresh]))
                        if len(forcein) < nforcein:
                            nforcein = len(forcein)
                            solvestart = time.time()
                            solution = pars.solve(exub=solution0[0],exlb=(0,0),
                                                  forcein=forcein,strategy='most',
                                                  tinitial=(opts.tinitial if first else opts.timprove)*60,
                                                  timprove=opts.timprove*60,
                                                  ttotal=ttotal)
                            ttotal -= (time.time()-solvestart)
                            progress.message("Component %d: ITERSOLVE proteins %d peptides %d forced %d solution %s"%(i,len(pr),len(pep),len(forcein),solution[0]))
                            first = False
                            if solution[1] != None:
                                solution0 = solution
                            elif solution[2] == False:
                                # timeout!!!
                                break
                        uthresh += 1
                        if len(forcein) == 0:
                            break
                    solution = solution0
                    progress.message("Component %d: proteins %d peptides %d solution %s"%(i,len(pr),len(pep),solution[0]))

                else:

                    solution = pars.solve(exub=(1e+20,1e+20),exlb=(0,0),strategy='most',
                                          tinitial=opts.tinitial*60,timprove=opts.timprove*60,
                                          ttotal=opts.ttotal*60)
                    progress.message("Component %d: proteins %d peptides %d solution %s"%(i,len(pr),len(pep),solution[0]))

                if solution[1]:
                    parsimony |= solution[1]

        else:
            # progress.stage("Component %d: proteins %d peptides %d solution %s"%(i,len(pr),len(pep),(len(pep),0)))
            pass

    if opts.pmatrix or opts.pannotate:
        fp = FixedPoint(elements,edges,parsimony,pepqval=psmfdr)
        fp.iterate(100)

elif opts.fixedpoint:

    assert not opts.peptidegroups

    dom = Dominator(peptides=elements,edges=edges,proteins=proteins)
    dom.dominate(tiebreakers=(pr2tbpep,tbpep2fdr,psmfilter.cmp),progress=progress)

    fp = FixedPoint(dom.peptides,dom.edges,dom.proteins,pepqval=psmfdr)
    fp.iterate(100)
    parsimony = set([pr for pr in dom.proteins if fp.proteinprob[pr]>=(1.0-opts.psmthresh)])

else:

    dom = Dominator(peptides=elements,edges=edges,proteins=proteins,pepqval=psmfdr)
    parsimony = set(proteins)

    if opts.pmatrix:
        fp = FixedPoint(elements,edges,parsimony,pepqval=psmfdr)
        fp.iterate(100)

newpars = set()
for pr in parsimony:
    equiv = dom.equivalentto[pr]
    if len(equiv) == 1:
        newpars.add(pr)
    else:
        tbequiv = set([pr])
        tbpr = pr
        for epr in equiv:
            if tbcmp(tbpr,epr) == 0:
                tbequiv.add(epr)
            if tbcmp(tbpr,epr) > 0:
                tbpr = epr
                tbequiv = set([epr])
        if len(tbequiv) == 1:
            newpars.add(tbpr)
            if tbpr != pr:
                dom.equivalentto[tbpr] = dom.equivalentto[pr]
                dom.containedby[tbpr] = dom.containedby[pr]
                del dom.equivalentto[pr]
                del dom.containedby[pr]
                if fp:
                    fp.proteinprob[tbpr]=fp.proteinprob[pr]
                    del fp.proteinprob[pr]
        else:
            # preference = sorted([(pri,accfn(pri)) for pri in tbequiv],key=cmp_to_key(accrule.prefer))
            # if min(list(map(itemgetter(0),preference))) != preference[0][0]:
            #    for i,(pri,desc) in enumerate(preference):
            #        print(i,pri,desc)
            #     print()
            # 
            # for i1,pri in enumerate(sorted(tbequiv,key=prsortkey.get)):
            #     print(i1,pri,prsortkey.get(pri),accfn(pri))
            # print()

            # tbequiv = list(tbequiv)
            # for i1 in range(len(tbequiv)):
            #     for i2 in range(i1+1,len(tbequiv)):
            #         print(tbequiv[i1],tbequiv[i2],prsortkey.get(tbequiv[i1]),prsortkey.get(tbequiv[i2]),\
            #               accrule.prefer((tbequiv[i1],accfn(tbequiv[i1])),(tbequiv[i2],accfn(tbequiv[i2]))))

            preferred = sorted([(pri,accfn(pri)) for pri in tbequiv],key=cmp_to_key(accrule.prefer))[0][0]
            newpars.add(preferred)
            # assert(preferred == pr)
            if preferred != pr:
                print("tb equiv switch! %d in %d out"%(preferred,pr))
                for i1,(pri,desc) in sorted([(pri,accfn(pri)) for pri in tbequiv],key=cmp_to_key(accrule.prefer)):
                     print(i1+1,pri,prsortkey.get(pri),desc,('+' if preferred == pri else ('-' if pr == pri else "")))
                dom.equivalentto[preferred] = dom.equivalentto[pr]
                dom.containedby[preferred] = dom.containedby[pr]
                del dom.equivalentto[pr]
                del dom.containedby[pr]
                if fp:
                    fp.proteinprob[preferred]=fp.proteinprob[pr]
                    del fp.proteinprob[pr]

            # print

parsimony = newpars

upep = defaultdict(int)
for pep in pep2pr:
    if len(pep2pr[pep]&parsimony) == 1:
        pr = next(iter(pep2pr[pep]&parsimony))
        upep[pr] += 1
if opts.parsimony:
    assert(set(upep.keys()) == parsimony)
minupep = defaultdict(list)
for k in parsimony:
    minupep[upep[k]].append(k)
minminupep = min(minupep)


if opts.pstats or (not opts.quiet):

    progress.message("Output parsimony statistics...")

    covered = set()
    for pri in parsimony:
        covered.update(edges[pri])
    uncovered = set(elements) - covered
    wuncovered = sum(map(psmcount.get,uncovered))

    print("Min unique peptides: %d (%d proteins)"%(minminupep,len(minupep[minminupep])), file=opts.outhandle)
    if opts.unique > minminupep:
        print("Proteins with < %d unique peptides:"%(opts.unique,), end=' ', file=opts.outhandle)
        for pri in minupep[minminupep]:
            print(praccfn(pri), end=' ', file=opts.outhandle)
        print("", file=opts.outhandle)
    print("%s uncovered: %d (%d spec) [ %.2f%% (%.2f%% spec) ]"%("Peptides" if not opts.peptidegroups else "Peptide groups",len(uncovered),wuncovered,100*float(len(uncovered))/len(elements),100*float(wuncovered)/nspec), file=opts.outhandle)
    print("Parsimony set: %d proteins (%.2f%%)"%(len(parsimony),100*float(len(parsimony))/(ngene if opts.bygene else nprot)), file=opts.outhandle)
    opts.outhandle.flush()

if opts.pdump:

    progress.message("Output parsimony accessions...")
    if not opts.output:
        opts.outhandle = open(opts.database[:-4]+'.%s.txt'%opts.extn,'w')
    for k,v in sorted(list(minupep.items()),reverse=True):
        for p in v:
            equiv = dom.equivalentto[p]
            equiv.remove(p)
            print(' '.join(map(praccfn,[p] + sorted(equiv))), file=opts.outhandle)
    if opts.outhandle != sys.stdout:
        opts.outhandle.close()

def pepinprotmark(pri,pepi):
    if pepi in edges[pri]:
        return str(psmcount[pepi])
    return " "

def coverage(pri):
    try:
        if not opts.bygene:
            return psmdb.proteinid_coverage(pri)
        return 0.0
    except RuntimeError:
        return 0.0

def pepseq(pepi):
    if opts.peptidegroups:
        return '|'.join([p.sequence for p in psmdb.peptidegroupid_peptides(pepi)])
    else:
        return psmdb.peptideid_sequence(pepi)

def protdesc(pri):
    if not opts.bygene:
        desc = psmdb.proteinid_description(pri)
    else:
        desc = psmdb.geneid_description(pri)
    if not desc:
        return ""
    if ',' in desc:
        return '"%s"'%desc
    return desc

def protorg(pri):
    if not opts.bygene:
        pr = psmdb.protein(pri)
    else:
        pr = psmdb.protein_group(pri)
    org = pr.getdata('organism')
    if not org:
        return ""
    if ',' in org:
        return '"%s"'%org
    return org

if opts.pmatrix:
    progress.message("Output parsimony matrix...")
    if not opts.output:
        opts.outhandle = open(opts.database[:-4]+'.%s.csv'%opts.extn,'w')
    comp = Components(peptides=elements, edges=edges, proteins=proteins)
    for pr,pep in sorted(comp,key=lambda t: (-len(t[1]),len(t[0]),min(list(map(praccfn,t[0]))))):
        prshown = copy.copy(pr)
        if len(prpep2align) > 0:
            from partial_order import partial_order
            peporder = partial_order(pr,elements,edges,prpep2align)
            peplist = sorted(pep,key=lambda pi: peporder.get(pi,1e+20))
        else:
            peplist = sorted(pep,key=pepseq)

        # first=True
        # for pri in sorted(pr,key=praccfn):
        # if psmdb.proteinid_description(pri):
        #     if first:
        #         print >>opts.outhandle, "Protein,Description"
        #     first = False
        #     print >>opts.outhandle, ">%s,%s"%(praccfn(pri),psmdb.proteinid_description(pri))

        print("Protein,Prob,Coverage,Status,"+",".join(map(pepseq,peplist)), file=opts.outhandle)
        for pri in sorted(parsimony&pr,key=lambda pri: (upep.get(pri),praccfn(pri)),reverse=True):
            print(','.join([praccfn(pri),"%.3f"%fp.proteinprob[pri],"%.2f"%coverage(pri),'P ']+list(map(partial(pepinprotmark,pri),peplist))+[_f for _f in [peppr2subst.get((p,pri)) for p in peplist] if _f]+[protorg(pri),protdesc(pri)]), file=opts.outhandle)
            prshown.remove(pri)
            for prj in sorted(dom.equivalentto[pri]&prshown,key=praccfn):
                if prj == pri:
                    continue
                print(','.join([praccfn(prj),"%.3f"%fp.proteinprob[prj],"%.2f"%coverage(prj),'PE']+list(map(partial(pepinprotmark,prj),peplist))+[_f for _f in [peppr2subst.get((p,prj)) for p in peplist] if _f]+[protorg(prj),protdesc(prj)]), file=opts.outhandle)
                prshown.remove(prj)
            for prj in sorted(dom.containedby[pri]&prshown,key=lambda p: (-len(pr2pep[p]),praccfn)):
                print(','.join([praccfn(prj),"%.3f"%fp.proteinprob[prj],"%.2f"%coverage(prj),'C ']+list(map(partial(pepinprotmark,prj),peplist))+[_f for _f in [peppr2subst.get((p,prj)) for p in peplist] if _f]+[protorg(prj),protdesc(prj)]), file=opts.outhandle)
                prshown.remove(prj)
        for pri in sorted((set(dom.equivalentto.keys())-parsimony)&prshown,key=praccfn):
            print(','.join([praccfn(pri),"%.3f"%fp.proteinprob[pri],"%.2f"%coverage(pri),'D ']+list(map(partial(pepinprotmark,pri),peplist))+[_f for _f in [peppr2subst.get((p,pri)) for p in peplist] if _f]+[protorg(pri),protdesc(pri)]), file=opts.outhandle)
            prshown.remove(pri)
            for prj in sorted(dom.equivalentto[pri]&prshown,key=praccfn):
                print(','.join([praccfn(prj),"%.3f"%fp.proteinprob[prj],"%.2f"%coverage(prj),'DE']+list(map(partial(pepinprotmark,prj),peplist))+[_f for _f in [peppr2subst.get((p,prj)) for p in peplist] if _f]+[protorg(prj),protdesc(prj)]), file=opts.outhandle)
                prshown.remove(prj)
            for prj in sorted(dom.containedby[pri]&prshown,key=praccfn):
                print(','.join([praccfn(prj),"%.3f"%fp.proteinprob[prj],"%.2f"%coverage(prj),'C ']+list(map(partial(pepinprotmark,prj),peplist))+[_f for _f in [peppr2subst.get((p,prj)) for p in peplist] if _f]+[protorg(prj),protdesc(prj)]), file=opts.outhandle)
                prshown.remove(prj)
        for pri in sorted(prshown,key=praccfn):
            print(','.join([praccfn(pri),"%.3f"%fp.proteinprob[pri],"%.2f"%coverage(pri),'  ']+list(map(partial(pepinprotmark,pri),peplist))+[_f for _f in [peppr2subst.get((p,pri)) for p in peplist] if _f]+[protorg(pri),protdesc(pri)]), file=opts.outhandle)
        print("", file=opts.outhandle)

    if opts.outhandle != sys.stdout:
        opts.outhandle.close()

if opts.pannotate:
    progress.message("Output parsimony annotations...")
    psmdbout = psmdb
    dominant=set(dom.equivalentto)
    comp = Components(peptides=peptides, edges=pr2pep, proteins=proteins)
    i = 0
    for pr,pep in sorted(comp,key=lambda t: -max([len(pr2pep[prid]) for prid in t[0]])):
        comppr = set()
        for pri in pr:
            if pri in dominant:
                if pri in parsimony:
                    prmd = {'parsimony:status':'P',
                            'component':i+1,
                            'fixedpoint:probability':fp.proteinprob[pri]}
                    psmdbout.protein(pri).updatedata(prmd)
                    comppr.add(pri)
                    for prj in dom.equivalentto[pri]:
                        if prj == pri:
                            continue
                        prmd = {'parsimony:equivalentto':psmdbout.protein(pri).accession,
                                'parsimony:status':'PE',
                                'component':i+1,
                                'fixedpoint:probability':fp.proteinprob[prj]}
                        psmdbout.protein(prj).updatedata(prmd)
                        comppr.add(prj)
                    for prj in dom.containedby[pri]:
                        prmd = {'parsimony:containedby':psmdbout.protein(pri).accession,
                                'parsimony:status':'C',
                                'component':i+1,
                                'fixedpoint:probability':fp.proteinprob[prj]}
                        psmdbout.protein(prj).updatedata(prmd)
                        comppr.add(prj)
                else:
                    prmd = {'parsimony:status':'D',
                            'component':i+1,
                            'fixedpoint:probability':fp.proteinprob[pri]}
                    psmdbout.protein(pri).updatedata(prmd)
                    comppr.add(pri)
                    for prj in dom.equivalentto[pri]:
                        if prj == pri:
                            continue
                        prmd = {'parsimony:equivalentto':psmdbout.protein(pri).accession,
                                'parsimony:status':'DE',
                                'component':i+1,
                                'fixedpoint:probability':fp.proteinprob[prj]}
                        psmdbout.protein(prj).updatedata(prmd)
                        comppr.add(prj)
                    for prj in dom.containedby[pri]:
                        prmd = {'parsimony:containedby':psmdbout.protein(pri).accession,
                                'parsimony:status':'C',
                                'component':i+1,
                                'fixedpoint:probability':fp.proteinprob[prj]}
                        psmdbout.protein(prj).updatedata(prmd)
                        comppr.add(prj)
        for pepi in pep:
            psmdbout.peptide(pepi).setdata('component',i+1)
        psmdbout.newProteinGroup('Component:%d'%(i+1,),comppr,'Component')
        i += 1
    psmdbout.newProteinGroup('Parsimony',parsimony)

if opts.pfilter:
    progress.message("Parsimony filter ...")
    if not opts.bygene:
        psmdb.proteinid_filter(keep=parsimony)
    else:
        psmdb.geneid_filter(keep=parsimony)
