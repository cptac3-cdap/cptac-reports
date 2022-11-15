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


valid_labels = """
NoLabel 1.0
iTRAQ-113 iTRAQ-114 iTRAQ-115 iTRAQ-116 iTRAQ-117 iTRAQ-118 iTRAQ-119 iTRAQ-121
TMT-126 TMT-127 TMT-128 TMT-129 TMT-130 TMT-131
TMT-126C TMT-127C TMT-127N TMT-128C TMT-128N TMT-129C TMT-129N TMT-130C TMT-130N TMT-131C TMT-131N
TMT6-126 TMT6-127 TMT6-128 TMT6-129 TMT6-130 TMT6-131
TMT10-126 TMT10-131
TMT10-127C TMT10-127N TMT10-128C TMT10-128N TMT10-129C TMT10-129N TMT10-130C TMT10-130N
TMT11-126 TMT11-131
TMT11-126C TMT11-127C TMT11-127N TMT11-128C TMT11-128N TMT11-129C TMT11-129N TMT11-130C TMT11-130N TMT11-131C TMT11-131N
TMT16-126
TMT16-126C TMT16-127C TMT16-127N TMT16-128C TMT16-128N TMT16-129C TMT16-129N TMT16-130C
TMT16-130N TMT16-131C TMT16-131N TMT16-132C TMT16-132N TMT16-133C TMT16-133N TMT16-134N
TMT18-126
TMT18-126C TMT18-127C TMT18-127N TMT18-128C TMT18-128N TMT18-129C TMT18-129N TMT18-130C
TMT18-130N TMT18-131C TMT18-131N TMT18-132C TMT18-132N TMT18-133C TMT18-133N TMT18-134N TMT18-134C TMT18-135N
""".split()

from psmdb.PSMDb import PSMDb

parser.add_option("-d","--database",type="file",dest="database",remember=True,
                  help="PSM database file",default="",
                  name="PSM Database",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-s","--samples",type="file",dest="samples",default="",
                  help="File grouping spectral files into samples.",
                  name="Spectra Groups",remember=True,
                  filetypes=[("All group files","*.xls;*.csv"),
                             ("CSV","*.csv"),
                             ("Excel","*.xls")])
parser.add_option("-S","--samplename",type="string",dest="samplename",default="All",
                  help="Sample name for all spectral files.",
                  name="Sample Name",remember=True)
parser.add_option("-L","--labelname",type="choice",dest="labelname",default="NoLabel",
                  help="Label for all spectral files.",choices=valid_labels,
                  name="Label Name",remember=True)
parser.add_option("--self",action="store_true",dest="self",default=False,remember=True,
                  name="Spectra file as sample name",help="Use spectra file as sample name. Default: False.")
parser.add_option("-R","--remove",action="store_true",dest="remove",default=False,remember=True,
                  name="Remove Unassigned?",help="Remove spectral files not assigned to a sample. Default: False.")
parser.add_option("-A","--append",action="store_true",dest="append",default=False,remember=True,
                  name="Append Samples?",help="Append, rather than replace samples. Default: False.")
parser.add_option("-o","--output",type="savefile",dest="output",default="",remember=True,
                  help="Output file.",
                  name="Output",
                  filetypes=[("PSM Database","*.psm")])
parser.add_option("-t","--template",action="store_true",dest="dump",default=False,remember=True,
                  name="Dump Filenames",help="Output template sample file. Default: False.")
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

    if not opts.database:
        parser.error('PSM Database (-d|--database) option is required.',**error_kwargs)
        continue

    break

progress = ProgressText(quiet=opts.quiet)

if opts.dump:
    psmdb = PSMDb(filename=opts.database,quiet=opts.quiet)
    wh = open(opts.database[:-4]+'.sample.csv','w')
    print("spectrafile,sample", file=wh)
    for sf in psmdb.spectrumfiles():
        print("%s,%s"%(sf.name,""), file=wh)
    wh.close()
    sys.exit(2)

progress.stage('Initialize PSM database...')
if not opts.output:
    psmdb = PSMDb(filename=opts.database,quiet=opts.quiet)
else:
    psmdb = PSMDb(filename=opts.output,clone=opts.database,quiet=opts.quiet)

reporters = set()
isotopemap = dict()
denominatormap = dict()
samplemap = defaultdict(set)
if opts.samples:
    progress.message("Read spectra file to sample mapping file...")
    if opts.samples.endswith(".xls"):
        t = XLSFileTable(opts.samples)
    else:
        t = CSVFileTable(opts.samples)
    # expect: filename, sample
    assert 'spectrafile' in t.headers()
    assert 'sample' in t.headers()
    counts = defaultdict(int)
    for r in t.rows():
        sf = psmdb.spectrumfileid(r['spectrafile'])
        # assert sf, "Bad spectra file: %s"%r['spectrafile']
        if not sf:
            continue
        ri = r.get('label','NoLabel')
        assert ri in valid_labels
        reporters.add(ri)
        try:
            offsets = ["-2:13C+13C","-2:13C+15N","-1:13C","-1:15N","+1:15N","+1:13C","+2:15N+13C","+2:13C+13C"]
            isotopes = list(map(float,list(map(r.get,offsets))))
            if (sf.id,ri) not in isotopemap:
                isotopemap[(sf.id,ri)] = isotopes
        except (ValueError,TypeError):
            pass
        for denominator in r.get("denominator","").split(';'):
            if not denominator:
                continue
            ordinal,label = denominator.split(":")
            assert label in valid_labels
            if (sf.id,ri) not in denominatormap:
                denominatormap[(sf.id,ri)] = set([denominator])
            else:
                denominatormap[(sf.id,ri)].add(denominator)

        samplemap[r['sample']].add((sf,ri))
        counts[sf.id] += 1
        counts[(sf.id,ri)] += 1
        if r.get('analyticalsample'):
            sf.setdata('analyticalsample',r.get('analyticalsample'))

    if not opts.remove and not opts.append:
        for sf in psmdb.spectrumfiles():
            if not opts.remove:
                assert counts[sf.id] > 0, "Spectrum file %s not in any sample"%sf.name
            for ri in reporters:
                if ri == 'NoLabel':
                    assert counts[(sf.id,ri)] < 2, "Spectrum file %s assigned to more than one sample"%sf.name
                else:
                    assert counts[(sf.id,ri)] < 2, "Reporter %s in spectrum file %s assigned to more than one sample"%(ri,sf.name)


elif opts.samplename and not opts.self:

    for sf in psmdb.spectrumfiles():
        samplemap[opts.samplename].add((sf,opts.labelname))

elif opts.self:

    for sf in psmdb.spectrumfiles():
        samplemap[sf.name].add((sf,opts.labelname))

if not opts.append:
    psmdb.deleteSamples()


for name,files in list(samplemap.items()):
    files = list(files)
    s = psmdb.newSample(name,list(map(itemgetter(0),files)),list(map(itemgetter(1),files)))
    if len(isotopemap) > 0:
        for f in s.files:
            isotopes = isotopemap[(f.spectrumfile.id,f.label.name)]
            f.setdata('isotope_corrections',isotopes)
    if len(denominatormap) > 0:
        for f in s.files:
            denominators = denominatormap.get((f.spectrumfile.id,f.label.name))
            if denominators and len(denominators) > 0:
                denominators = sorted([d.split(':') for d in denominators])
                denominators = list(zip(list(map(itemgetter(1),denominators)),list(map(int,list(map(itemgetter(0),denominators))))))
                f.setdata('ratio_denominator',denominators)

if opts.remove:
    psmdb.deleteNoSampleSpectrumFiles()

if not opts.quiet:
    psmdb.dump(sys.stdout)
