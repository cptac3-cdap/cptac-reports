#!/bin/env python27
import sys, csv
from scipy.stats import hypergeom

from scipy import special
from math import exp, floor, log
from psmdb.PSMDb import PSMDb
from collections import Counter, defaultdict
from peptidescan.fastalib import fasta
from fragger.MassTable import AverageMassTable

def countpr2(f):
    c = defaultdict(lambda: defaultdict(int))
    dvalue = set()
    db = PSMDb(filename=f)
    pep2fdr      = db.peptideid_psmminvalue()
    minqvalue = min([v for v in iter(pep2fdr.values()) if v>0])
    for pr in db.proteins():
        maxminfdr = -1
        count = 0
        score = 0
        for pep in pr.unshared_peptides():
            maxminfdr = max(maxminfdr,pep2fdr[pep.id])
            count += 1
            if pep2fdr[pep.id] == 0:
                score += -log(minqvalue/2.0,2.0)
            else:
                score += -log(pep2fdr[pep.id],2.0)
        if maxminfdr == 0:
            maxminfdr = minqvalue/2.0
        score1 = -log(maxminfdr,2.0)*log(count,2.0)
        c[score1][pr.decoy] += 1
        dvalue.add(pr.decoy)
    lastscore = max(c.keys())
    for score in sorted(list(c.keys()),reverse=True)[1:]:
        for decoy in dvalue:
            c[score][decoy] += c[lastscore][decoy]
        lastscore = score
    return c

def countpr3(f,protein_group):
    c = defaultdict(lambda: defaultdict(int))
    db = PSMDb(filename=f)
    for pr in db.proteins():
        c[protein_group[pr.accession]][pr.decoy] += 1
        c[0][pr.decoy] += 1
    return c

def countpr4(f):
    c = defaultdict(lambda: defaultdict(int))
    dvalue = set()
    db = PSMDb(filename=f)
    pep2fdr      = db.peptideid_psmminvalue()
    for pr in db.proteins():
        unshared = sorted([pep2fdr[p.id] for p in pr.unshared_peptides()])
        c[unshared[1]][pr.decoy] += 1
        dvalue.add(pr.decoy)
    for decoy in dvalue:
        lastval = None
        for unsh in sorted(c.keys()):
            if lastval in c:
                c[unsh][decoy] += c[lastval][decoy]
            lastval = unsh
    return c

def countpr(f):
    c = defaultdict(lambda: defaultdict(int))
    dvalue = set()
    db = PSMDb(filename=f)
    for pr in db.proteins():
        unshared = len(pr.unshared_peptides())
        # print unshared,pr.accession
        c[unshared][pr.decoy] += 1
        dvalue.add(pr.decoy)
    for unsh in range(max(c.keys()),min(c.keys())-1,-1):
        for decoy in dvalue:
            c[unsh][decoy] += c[unsh+1][decoy]
    return c

def countpr1(f,accs=None):
    c = defaultdict(lambda: defaultdict(int))
    db = PSMDb(filename=f)
    for pr in db.proteins():
        if accs == None or pr.accession in accs:
            c[-1][pr.decoy] += 1
    return c

def countpr1g(f,names=None):
    c = defaultdict(lambda: defaultdict(int))
    db = PSMDb(filename=f)
    for g in db.genes():
        if names == None or g.name in names:
            c[-1][g.getdata('decoy')] += 1
    return c

class ProteinByLength:
    def __init__(self,size=100,start=0):
        self.start=start
        self.size=size
    def group(self,l,mw):
        return (l-self.start)/self.size

class ProteinByMolWt:
    def __init__(self,size=10000.0,start=0.0):
        self.start=start
        self.size=size
    def group(self,l,mw):
        return int(math.floor((mw-self.start)/self.size))

class ProteinEqualBinByLength:
    def __init__(self,prlen_file,bins=20):
        self.prgroup = {}
        prlen = Counter()
        total = 0
        for row in csv.DictReader(prlen_file):
            prlen[row['length']] += 1
            total += 1
        perbin = total/bins
        curprs = 0
        curbin = 1
        for k,v in sorted(prlen.items()):
            self.prgroup[k] = curbin
            curprs += v
            if curprs > perbin:
                curprs = 0
                curbin += 1
    def group(self,l,mw):
        return self.prgroup[l]

class ProteinEqualBinByMW:
    def __init__(self,prlen_file,bins=20):
        self.prgroup = {}
        prlen = Counter()
        total = 0
        for row in csv.DictReader(open(prlen_file)):
            prlen[int(float(row['molwt']))] += 1
            total += 1
        perbin = total/bins
        curprs = 0
        curbin = 1
        for k,v in sorted(prlen.items()):
            self.prgroup[k] = curbin
            curprs += v
            if curprs > perbin:
                curprs = 0
                curbin += 1
    def group(self,l,mw):
        return self.prgroup[int(float(mw))]

def prgroups(lenmw_file,group_function):
    prgroup = Counter()
    for row in csv.DictReader(open(lenmw_file)):
        l = row['length']; mw = row['molwt']
        if row['accession'].startswith('XXX_'):
            continue
        prgroup[group_function(l,mw)] += 1
        prgroup[0] += 1
    return prgroup

def logchoose(n, k):
    lgn1 = special.gammaln(n+1)
    lgk1 = special.gammaln(k+1)
    lgnk1 = special.gammaln(n-k+1)
    # print lgn1,lgk1,lgnk1
    return lgn1 - (lgnk1 + lgk1)

def gauss_hypergeom(x, r, b, n):
    return exp(logchoose(r, x) +
               logchoose(b, n-x) -
               logchoose(r+b, n))

def mayu(nprot,ntarget,ndecoy):

    if ntarget < ndecoy:
        return ntarget,100.0
    prob = {}
    for fp in range(0,min(ntarget,ndecoy)+1):
        tp = ntarget-fp
        nottp = nprot-tp
        prob[fp] = gauss_hypergeom(fp,nottp,tp,ndecoy)
    totalprob = sum(prob.values())
    expectedfp = sum(fp*p/totalprob for fp,p in list(prob.items()))
    fdr = 100*expectedfp/ntarget
    return expectedfp,fdr

if __name__ == '__main__':

    bygene = False
    if sys.argv[1] == '--bygene':
        bygene = True
        sys.argv.pop(1)
    praccs = None; f1 = None
    if sys.argv[1] == '--accessions':
        praccs = set()
        for l in open(sys.argv[2]):
            sl = l.split()
            praccs.add(sl[0])
        f1 = sys.argv[2]
        sys.argv.pop(1)
        sys.argv.pop(1)
    summary = False
    if sys.argv[1] == '--summary':
        summary = True
        sys.argv.pop(1)
    score = False
    if sys.argv[1] == '--score':
        score = True
        sys.argv.pop(1)
    groups = False
    if sys.argv[1] == '--groups':
        groups = True
        sys.argv.pop(1)
        lenmw_file = sys.argv[1]
        prgroup = ProteinEqualBinByMW(lenmw_file)
        prcounts = prgroups(lenmw_file,prgroup.group)
        pr2group = {}
        for row in csv.DictReader(open(lenmw_file)):
            pr2group[row['accession']] = prgroup.group(row['length'],row['molwt'])
    byfdr = False
    if sys.argv[1] == '--byfdr':
        byfdr = True
        sys.argv.pop(1)

    nprot = int(sys.argv[1])
    for f in sys.argv[2:]:
        mayuresults = {}
        if summary:
            if bygene:
                counts = countpr1g(f,praccs)
            else:
                counts = countpr1(f,praccs)
        elif score:
            counts = countpr2(f)
        elif groups:
            counts = countpr3(f,pr2group)
        elif byfdr:
            counts = countpr4(f)
        else:
            counts = countpr(f)
        for unsh in sorted(list(counts.keys()),reverse=(True if not byfdr else False)):
            ntarget = counts[unsh][False]
            ndecoy = counts[unsh][True]
            if ntarget == 0:
                continue
            if groups:
                nprot = prcounts[unsh]
            expfp,fdr = mayu(nprot,ntarget,ndecoy)
            mayuresults[unsh] = [nprot,ntarget,ndecoy,fdr]
            # print(unsh,nprot,ntarget,ndecoy,fdr)
        if not groups:
            lastfdr = 100.0
            for unsh in sorted(list(counts.keys()),reverse=(False if not byfdr else True)):
                if unsh not in mayuresults:
                    continue
                if mayuresults[unsh][3] >= lastfdr:
                    del mayuresults[unsh]
                else:
                    lastfdr=mayuresults[unsh][3]
        totexptp=0
        for unsh in sorted(list(counts.keys()),reverse=True):
            if unsh not in mayuresults:
                continue
            nprot,ntarget,ndecoy,fdr = mayuresults[unsh]
            if unsh > 0:
                totexptp += ntarget*(1.0-(fdr/100.0))
            print((f if f1 == None else f1),"%s"%unsh,nprot,ntarget,ndecoy,"%.4f%%"%(100*float(ndecoy)/ntarget),"%.4f%%"%fdr,"%.0f"%floor(ntarget*(1-(fdr/100.0))))
        # fdr = 100.0*(ntarget-totexptp)/ntarget
        # print f,"%.4f"%-1,nprot,ntarget,ndecoy,"%.2f%%"%(100*float(ndecoy)/ntarget),"%.2f%%"%fdr,"%.0f"%floor(ntarget*(1-fdr/100.0))
