#!/usr/bin/env python27

import sys
import re
from psmdb.PSMDb import PSMDb
import csv
from collections import defaultdict

assays = defaultdict(list)
for f in sys.argv[2:]:
    for r in csv.DictReader(open(f)):
        genewords = r['Gene'].split()
        gene = genewords[0]
        upacc = genewords[-1]
        peptide = r['Proteins and peptides with assays']
        assayid = r['CPTAC ID']
        d = dict()
        d['Gene'] = gene
        d['Protein'] = upacc
        d['AssayID'] = assayid
        assays[peptide].append(d)
# print assays

db = PSMDb(filename=sys.argv[1])
trans=db.begin()
for p,l in assays.items():
    pep = db.peptide(sequence=p,connection=trans)
    if not pep:
        continue

    for d in l:

        anypr = False

        for pr in pep.proteins():

            if pr.accession == d['Protein'] or \
                   pr.accession == d['Gene']:

                if d['AssayID'] not in pr.getdata('Assay',[]):
                    pr.appenddata('Assay',d['AssayID'])
                    anypr = True

            for pg in pr.groups():
                if pg.type == 'Gene' and pg.name == d['Gene']:
                    if d['AssayID'] not in pg.getdata('Assay',[]):
                        pg.appenddata('Assay',d['AssayID'])
                        anypr = True

        pep.appenddata('Assay',d['AssayID'])

        if not anypr:
            print >>sys.stderr, "Could not match assay %s on peptide %s to any protein or gene"%(d['AssayID'],p)

trans.commit(close=True)
