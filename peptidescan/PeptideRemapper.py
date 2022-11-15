
from __future__ import print_function

import sys,re
from .fastalib import fasta
from .ProteinHitDb import ProteinDatabase
from .OutOfCoreTable import *
from collections import defaultdict

__all__ = [ 'AccessionRules', 'AccessionRuleFactory',
            'AccessionRuleAuto', 'TestAccessionRules',
            'PeptideRemapper' ]

def cmp(a,b):
    if a < b:
        return -1
    elif a > b:
        return 1
    return 0

def AccessionRuleFactory(name):
    for cls in AccessionRuleClasses:
        if cls.name == name:
            return cls()
    return None

def AccessionRuleAuto(seqdb):
    rules = []
    for cls in AccessionRuleClasses:
        if cls.name == 'FirstWord':
            continue
        rules.append(cls())
    for i,(d,s) in enumerate(fasta(seqdb)):
        if i >= 10:
            break
        r = TestAccessionRules(d,rules)
        if r != None:
            return r
    return AccessionRuleFactory('FirstWord')

def TestAccessionRules(description,rules=None):
    if description == None:
        return None
    if rules == None:
        rules = []
        for cls in AccessionRuleClasses:
            if cls.noauto:
                continue
            rules.append(cls())
        rules.append(FirstWord())
    for r in rules:
        acc,dsc = r(description)
        # print acc,dsc,r
        if acc:
            return r
    return None

class RegExAcc:
    noauto = False
    def __call__(self,description):
        description = self.normalize(description)
        acc = self.accrex.search(description)
        # print self.__class__.__name__,description,acc.group(self.accgrp) if acc else "-"
        if acc == None:
            return None,None
        if self.dscrex:
            dsc = self.dscrex.search(description)
        else:
            return acc.group(self.accgrp),None
        if dsc == None:
            raise RuntimeError("Description regex doesn't match:\n  "+description)
        return acc.group(self.accgrp),dsc.group(self.dscgrp)
    def normalize(self,description):
        return description.lstrip('>')
    def split(self,description):
        return self(description)
    def shortacc(self,acc):
        return None
    def url(self,description):
        # print >>sys.stderr, self
        if not hasattr(self,'urlrex') or not hasattr(self,'urltmpl'):
            return None
        # print >>sys.stderr, description
        description = self.normalize(description)
        urlfields = self.urlrex.search(description)
        # print urlfields, urlfields.groupdict()
        if not urlfields:
            return None
        return self.urltmpl%urlfields.groupdict()
    def org(self,description):
        if not hasattr(self,'orgrex') or not hasattr(self,'orggrp'):
            return None
        description = self.normalize(description)
        dsc = self.orgrex.search(description)
        if not dsc:
            return None
        return dsc.group(self.orggrp)
    def gene(self,description):
        if not hasattr(self,'generex') or not hasattr(self,'genegrp'):
            return None
        description = self.normalize(description)
        dsc = self.generex.search(description)
        if not dsc:
            return None
        return dsc.group(self.genegrp)
    def entry(self,description):
        if not hasattr(self,'entryrex') or not hasattr(self,'entrygrp'):
            return None
        description = self.normalize(description)
        dsc = self.entryrex.search(description)
        if not dsc:
            return None
        return dsc.group(self.entrygrp)
    @staticmethod
    def prefer(pr1,pr2):
        pr1id,pr1defline = pr1
        pr2id,pr2defline = pr2
        cmps = set()
        for badword in "predicted hypothetical isoform uncharacterized putative cDNA homolog pseudogene LOC\\d+ like \\w+-like readthrough R000\\w+ S000\\w+ XXX_\\w+".split():
            c = cmp(1*(re.search(r'\b%s\b'%badword,pr1defline,re.I)!=None),
                    1*(re.search(r'\b%s\b'%badword,pr2defline,re.I)!=None))
            cmps.add(c)
        if -1 in cmps and 1 not in cmps:
            return -1
        if 1 in cmps and -1 not in cmps:
            return 1
        if pr1defline in pr2defline:
            return -1
        elif pr2defline in pr1defline:
            return 1
        words1 = pr1defline.split()
        words2 = pr2defline.split()
        if words1[1:-1] == words2[1:-1]:
            c = cmp(words1[-1],words2[-1])
            if c != 0:
                return c
        if words1[0] != words2[0]:
            return cmp(words1[0],words2[0])
        return cmp(pr1,pr2)
    def prefermatrix(self,prlist):
        for i,pr in enumerate(prlist):
            print(i,pr[0],pr[1])
        print(" ", end=' ')
        for i,pr in enumerate(prlist):
            print(i, end=' ')
        print()
        for i,pr in enumerate(prlist):
            print(i, end=' ')
            for pr1 in prlist:
                print("." if pr == pr1 else ("<" if (self.prefer(pr,pr1)<0) else "^"), end=' ')
            print()

class FirstWord(RegExAcc):
    noauto = True
    name = 'FirstWord'
    accrex = re.compile(r'^(\S*)\s+')
    accgrp = 1
    dscrex = re.compile(r'^(\S*)\s+?(.*)$')
    dscgrp = 2

class SecondAccFirstWord(RegExAcc):
    noauto = True
    name = 'SecondAccFirstWord'
    accrex = re.compile(r'^(([A-Za-z0-9_]+)\|([A-Za-z0-9_]+)\S*)\s+')
    accgrp = 3
    dscrex = re.compile(r'^(\S*)\s+?(.*)$')
    dscgrp = 2

class IPIShortAcc(RegExAcc):
    name = 'IPI'
    accrex = re.compile(r'^IPI:(IPI[^.\s|]*)')
    accgrp = 1
    dscrex = re.compile(r'^([^\s]*\s)(Tax_Id=\d+\s)(Gene_Symbol=(\S+)\s)(.*)$')
    dscgrp = 5
    urlrex = re.compile(r'^IPI:(?P<acc>IPI[^.\s|]*)')
    urltmpl = 'http://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=IPI&id=%(acc)s'
    generex = dscrex
    genegrp = 4

class IPIGene(RegExAcc):
    name = 'IPIGene'
    accrex = re.compile(r'\sGene_Symbol=([^-\s]+)')
    accgrp = 1
    dscrex = None
    dscgrp = None
    urlrex = re.compile(r'\s+Tax_Id=(?P<taxid>\d+)\s+Gene_Symbol=(?P<sym>[^-\s]+)')
    urltmpl = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&term=%(sym)s[sym] AND %(taxid)s[taxid]'

class UniProtAcc(RegExAcc):
    name = 'UniProt'
    accrex = re.compile(r'^(sp|tr)\|([^|-]+)(-\d+)?\|')
    accgrp = 2
    dscrex = re.compile(r'^(sp|tr)\|[^|-]+(-\d+)?\|\S*\s+(Isoform \S+( \S+)* of )?(.*)\s+OS=.*$')
    dscgrp = 5
    urlrex = re.compile(r'^(sp|tr)\|(?P<acc>[^|-]+)(-\d+)?\|')
    urltmpl = 'http://uniprot.org/uniprot/%(acc)s'
    # entryrex = re.compile(r'^(sp|tr)\|([^|-]+)(-\d+)?\|')
    # entrygrp = 2
    generex = re.compile(r'\s+GN=(\S+)')
    genegrp = 1
    orgrex = re.compile(r'\s+OS=([^=]*)( [A-Z][A-Z]=.*)?$')
    orggrp = 1
    shacc = re.compile(r'^(.*)-\d+$')
    @staticmethod
    def prefer(pr1,pr2):
        return UniProtIsoformAcc.prefer(pr1,pr2)
    def shortacc(self,acc):
        m = UniProtAcc.shacc.search(acc)
        if m:
            return m.group(1)
        return None

class UniProtID(RegExAcc):
    name = 'UniProtID'
    accrex = re.compile(r'^(sp|tr)\|([^|-]+)(-\d+)?\|(\S+)')
    accgrp = 4
    dscrex = re.compile(r'^(sp|tr)\|[^|-]+(-\d+)?\|\S*\s+(.*)\s+OS=.*$')
    dscgrp = 3
    urlrex = re.compile(r'^(sp|tr)\|(?P<acc>[^|-]+)(-\d+)?\|')
    urltmpl = 'http://uniprot.org/uniprot/%(acc)s'
    generex = re.compile(r'\s+GN=(\S+)')
    genegrp = 1
    orgrex = re.compile(r'\s+OS=([^=]*)( [A-Z][A-Z]=.*)?$')
    orggrp = 1
    @staticmethod
    def prefer(pr1,pr2):
        return UniProtIsoformAcc.prefer(pr1,pr2)

class UniProtIsoformAcc(RegExAcc):
    name = 'UniProtIsoform'
    accrex = re.compile(r'^(sp|tr)\|([^|-]+(-\d+)?)\|.*\s+OS=.*$')
    accgrp = 2
    dscrex = re.compile(r'^(sp|tr)\|[^|-]+(-\d+)?\|\S*\s+(.*)\s+OS=.*$')
    dscgrp = 3
    urlrex = re.compile(r'^(sp|tr)\|(?P<acc>[^|-]+)(-\d+)?\|')
    urltmpl = 'http://uniprot.org/uniprot/%(acc)s'
    generex = re.compile(r'\s+GN=(\S+)')
    genegrp = 1
    entryrex = re.compile(r'^(sp|tr)\|([^|-]+)(-\d+)?\|')
    entrygrp = 2
    pe = re.compile(r' PE=(\d) ')
    gn = re.compile(r' GN=(\S+) ')
    orgrex = re.compile(r'\s+OS=([^=]*)( [A-Z][A-Z]=.*)?$')
    orggrp = 1
    shacc = re.compile(r'^(.*)-\d+$')
    def shortacc(self,acc):
        m = UniProtIsoformAcc.shacc.search(acc)
        if m:
            return m.group(1)
        return None
    @staticmethod
    def prefer(pr1,pr2):
        pr1id,pr1defline = pr1
        pr2id,pr2defline = pr2
        if pr1defline == pr2defline:
            return cmp(pr1id,pr2id)
        # print pr1,pr2,
        m1 = UniProtIsoformAcc.accrex.search(pr1defline.lstrip('>'))
        m2 = UniProtIsoformAcc.accrex.search(pr2defline.lstrip('>'))
        if not m1:
            pr1db = ""; pr1acc = ""; pr1iso = ""
        else:
            pr1db = m1.group(1); pr1acc = m1.group(2); pr1iso = m1.group(3)
        if not m2:
            pr2db = ""; pr2acc = ""; pr2iso = ""
        else:
            pr2db = m2.group(1); pr2acc = m2.group(2); pr2iso = m2.group(3)
        if pr1db != pr2db:
            # print "db",-1 if pr1db == 'sp' else 1
            return -1 if pr1db == 'sp' else 1
        if pr1db == pr2db and pr1acc == pr2acc:
            if pr1iso == "":
                # print "iso",-1
                return -1
            elif pr2iso == "":
                # print "iso",1
                return 1
        pe1 = UniProtIsoformAcc.pe.search(pr1defline)
        pe2 = UniProtIsoformAcc.pe.search(pr2defline)
        pe1 = (int(pe1.group(1)) if pe1 else 1e+20)
        pe2 = (int(pe2.group(1)) if pe2 else 1e+20)
        if pe1 != pe2:
            # print "pe",cmp(pe1,pe2)
            return cmp(pe1,pe2)
        gn1 = UniProtIsoformAcc.gn.search(pr1defline)
        gn2 = UniProtIsoformAcc.gn.search(pr2defline)
        if gn1 and not gn2:
            # print "gn",-1
            return -1
        if not gn1 and gn2:
            # print "gn",1
            return 1
        pr1defline = pr1defline.lower()
        pr2defline = pr2defline.lower()
        if 'isoform' not in pr1defline and 'isoform' in pr2defline:
            # print "isoform",-1
            return -1
        if 'fragment' not in pr1defline and 'fragment' in pr2defline:
            # print "isoform",-1
            return -1
        if 'isoform' in pr1defline and 'isoform' not in pr2defline:
            # print "isoform",1
            return 1
        if 'fragment' in pr1defline and 'fragment' not in pr2defline:
            # print "isoform",1
            return 1
        # print "id",cmp(pr1id,pr2id)
        if pr1acc != pr2acc:
            if pr1acc and not pr2acc:
                return -1
            if pr2acc and not pr1acc:
                return 1
            return cmp(pr1acc,pr2acc)
        return cmp(pr1id,pr2id)

class BroadArtifact(RegExAcc):
    name = 'BroadArtifact'
    accrex = re.compile(r'^sp\|(B99\d+)\|')
    accgrp = 1
    dscrex = re.compile(r'^sp\|B99\d+\|\S*\s+(.*?)(\s+\[([^]]*)\]\s*)?$')
    dscgrp = 1
    orgrex = re.compile(r'^.* +\[([^]]*)\] *$')
    orggrp = 1

class UniProtGene(RegExAcc):
    name = 'UniProtGene'
    accrex = re.compile(r'\sGN=(\S+)')
    accgrp = 1
    dscrex = None
    dscgrp = None
    orgrex = re.compile(r'\s+OS=([^=]*)( [A-Z][A-Z]=.*)?$')
    orggrp = 1
    urlrex = re.compile(r'\sOS=(?P<species>.*)\sGN=(?P<sym>\S+)')
    urltmpl = 'http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&term=%(sym)s[sym] AND %(species)s[orgn]'

class RefSeqAcc(RegExAcc):
    name = 'RefSeqAcc'
    accrex = re.compile(r'^((gi\|(\d+)\|)?ref\|)?([NYXZA]P_[^| ]+)\|?')
    accgrp = 4
    dscrex = re.compile(r'^((gi\|(\d+)\|)?ref\|)?([NYXZA]P_[^| ]+)\|? *(.*?) +(GN=(\S+) +)?\[[^]]*\] *$')
    dscgrp = 5
    orgrex = re.compile(r'^((gi\|(\d+)\|)?ref\|)?([NYXZA]P_[^| ]+)\|? *(.*) +\[([^]]*)\] *$')
    orggrp = 6
    generex = re.compile(r'\s+GN=(\S+)')
    genegrp = 1
    shacc = re.compile(r'^(.*)\.\d+$')
    def shortacc(self,acc):
        m = RefSeqAcc.shacc.search(acc)
        if m:
            return m.group(1)
        return None
    @staticmethod
    def prefer(pr1,pr2):
        pr1id,pr1defline = pr1
        pr2id,pr2defline = pr2
        if pr1defline == pr2defline:
            return cmp(pr1id,pr2id)
        # print pr1,pr2,
        m1 = RefSeqAcc.accrex.search(pr1defline.lstrip('>'))
        m2 = RefSeqAcc.accrex.search(pr2defline.lstrip('>'))
        if not m1:
            pr1acc = "";
        else:
            pr1acc = m1.group(RefSeqAcc.accgrp);
            pr1db = pr1acc[:2]
        if not m2:
            pr2acc = "";
        else:
            pr2acc = m2.group(RefSeqAcc.accgrp);
            pr2db = pr2acc[:2]
        if pr1acc == "" or pr2acc == "":
            if pr1acc != pr2acc:
                return cmp(1*(pr1acc == ""),1*(pr2acc == ""))
            return cmp(db1id,db2id)
        if pr1db != pr2db:
            return cmp(1*(pr1db!='NP'),1*(pr2db!='NP'))
        for badword in "predicted hypothetical isoform uncharacterized putative".split():
            c = cmp(1*(re.search(r'\b%s\b',pr1defline,re.I)!=None),1*(re.search(r'\b%s\b',pr2defline,re.I)!=None))
            if c != 0:
                return c
        return cmp(pr1id,pr2id)

class RefSeqGI(RegExAcc):
    name = 'RefSeqGI'
    accrex = re.compile(r'^gi\|(\d+)\|ref\|([^|]+)\|')
    accgrp = 1
    dscrex = re.compile(r'^gi\|(\d+)\|ref\|([^|]+)\| (.*) [\w+]$')
    dscgrp = 3
    orgrex = re.compile(r'^gi\|(\d+)\|ref\|([^|]+)\| (.*) \[([^]]*)\] *$')
    orggrp = 4

class SGDAcc(RegExAcc):
    name = 'SGD'
    accrex = re.compile(r'^([A-Z0-9]+) (\w+) SGDID:(S\d+).*"(.*)"$')
    accgrp = 1
    dscrex = re.compile(r'^([A-Z0-9]+) (\w+) SGDID:(S\d+).*"(.*)"$')
    dscgrp = 4
    generex = re.compile(r'^([A-Z0-9]+) (\w+) SGDID:(S\d+).*"(.*)"$')
    genegrp = 2
    urlrex = re.compile(r'^(?P<acc>[A-Z0-9]+) (?P<gene>\w+) SGDID:(S\d+) .* "(.*)"$')
    urltmpl = 'http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=%(acc)s'

class UCSCKGAcc(FirstWord):
    noauto = True
    name = 'UCSCKGAcc'

class Gene(FirstWord):
    noauto = True
    name = 'Gene'

class OrthoGene(FirstWord):
    noauto = True
    name = 'OrthoGene'
    shacc = re.compile(r'^(.*)\(\d+\)$')
    def shortacc(self,acc):
        m = OrthoGene.shacc.search(acc)
        if not m:
            return None
        return m.group(1)
    @staticmethod
    def prefer(pr1,pr2):
        c = cmp(('|' in pr1[1])*1,('|' in pr2[1])*1)
        if c != 0:
            return -c
        return Gene.prefer(pr1,pr2)

class MultiAcc(RegExAcc):
    classes = []

    def __init__(self):
        self.rules = [c() for c in self.classes]
        self.acc2rule = dict()

    def __call__(self,description):
        for r in self.rules:
            acc,desc = r(description)
            if acc:
                if acc not in self.acc2rule:
                    self.acc2rule[acc] = r
                return acc,desc
        return None,None

    def getacc(self,description):
        for r in self.rules:
            acc,desc = r(description)
            if acc:
                assert acc in self.acc2rule
                assert self.acc2rule[acc] == r
                return acc
        raise RuntimeError('Can\'t get accession')

    def org(self,description):
        acc = self.getacc(description)
        return self.acc2rule[acc].org(description)

    def gene(self,description):
        acc = self.getacc(description)
        return self.acc2rule[acc].gene(description)

    def entry(self,description):
        acc = self.getacc(description)
        return self.acc2rule[acc].entry(description)

    def url(self,description):
        acc = self.getacc(description)
        return self.acc2rule[acc].url(description)

    def shortacc(self,acc):
        return self.acc2rule[acc].shortacc(acc)

    def prefer(self,pr1,pr2):
        for r in self.rules:
            acc1,desc1 = r(pr1[1])
            acc2,desc2 = r(pr2[1])
            if acc1 and not acc2:
                return -1
            elif acc2 and not acc1:
                return 1
            elif acc1 and acc2:
                return r.prefer(pr1,pr2)
        return cmp(pr1,pr2)

class UniProtIsoformRefSeq(MultiAcc):
    name = 'UniProtIsoformRefSeq'
    classes = [ UniProtIsoformAcc, RefSeqAcc ]
class UniProtRefSeq(MultiAcc):
    name = 'UniProtRefSeq'
    classes = [ UniProtAcc, RefSeqAcc ]
class RefSeqUniProt(MultiAcc):
    name = 'RefSeqUniProt'
    classes = [ RefSeqAcc, UniProtAcc ]
class RefSeqUniProtIsoform(MultiAcc):
    name = 'RefSeqUniProtIsoform'
    classes = [ RefSeqAcc, UniProtIsoformAcc ]
class RefSeqUniProtIsoformBroadArtifact(MultiAcc):
    name = 'RefSeqUniProtIsoformBroadArtifact'
    classes = [ RefSeqAcc, UniProtIsoformAcc, BroadArtifact ]

class UniProtIsoformDecoy(MultiAcc):
    name = 'UniProtIsoform+Decoy'
    classes = [ UniProtIsoformAcc, FirstWord ]

class OrthoRefSeqUniProtIsoform(RefSeqUniProtIsoform):
    noauto = True
    name = 'OrthoRefSeqUniProtIsoform'
    def prefer(self,pr1,pr2):
        c = cmp(('; gi|' in pr1[1])*1,('; gi|' in pr2[1])*1)
        if c != 0:
            return -c
        return RefSeqUniProtIsoform.prefer(self,pr1,pr2)

class OrthoRefSeqUniProtIsoformBroadArtifact(RefSeqUniProtIsoformBroadArtifact):
    noauto = True
    name = 'OrthoRefSeqUniProtIsoformBroadArtifact'
    def prefer(self,pr1,pr2):
        c = cmp(('; gi|' in pr1[1])*1,('; gi|' in pr2[1])*1)
        if c != 0:
            return -c
        return RefSeqUniProtIsoformBroadArtifact.prefer(self,pr1,pr2)

AccessionRuleClasses = [ FirstWord,
                         IPIShortAcc, UniProtIsoformAcc, RefSeqAcc, SGDAcc,
                         UniProtAcc, UniProtID, RefSeqGI,
                         IPIGene, UniProtGene,
                         UCSCKGAcc, Gene, OrthoGene,
                         UniProtIsoformRefSeq,
                         UniProtRefSeq,
                         RefSeqUniProt,
                         RefSeqUniProtIsoform,
                         RefSeqUniProtIsoformBroadArtifact,
                         OrthoRefSeqUniProtIsoform,
                         OrthoRefSeqUniProtIsoformBroadArtifact,
                         SecondAccFirstWord,
                         UniProtIsoformDecoy,
                         ]

AccessionRules = [c.name for c in AccessionRuleClasses]

from collections import defaultdict
from operator import itemgetter

class PeptideRemapper:

    # ptmdiffs = [ 0.0, -1.0, 1.0, -14.0, -16.0, -28.0, -32.0, -42.0 ]

    def __init__(self,peptides,seqdb,accfn=FirstWord(),verbose=False,blocksize=None,preprocess=False,translation=None,aasubst=0,dnamut=0):
        self.accfn = accfn
        self.seqdb = seqdb
        self.verbose = verbose
        self.blocksize = blocksize
        self.preprocess = preprocess
        self.peptides = peptides
        self.dnamut = dnamut
        self.aasubst = aasubst
        self.protdb = ProteinDatabase()
        self.trans = translation
        if self.trans:
            self.trans = str(self.trans)
            assert self.trans in 'FA'
        self.domap()

    def domap(self):
        from .Options import Options
        from .Run import PeptideScanCommandLine
        opt = Options()
        opt.set('i',self.seqdb)
        opt.set('C',5)
        if self.verbose:
            opt.set('v')
        if self.trans:
            opt.set('T',self.trans)
        if self.dnamut:
            opt.set('K',self.dnamut)
            opt.set('x',3)
        if self.aasubst:
            opt.set('K',self.aasubst*3)
            opt.set('x',3)
        pscmd = PeptideScanCommandLine(opt)
        kwargs = {'preprocess': self.preprocess}

        peptides = OutOfCoreDistinct(self.peptides)
        # print '\n'.join(map(repr,peptides))
        self.protdb.load_peptides(peptides)
        npep = self.protdb.peptide_count()
        if not self.blocksize:
            self.blocksize = max(npep,1)
            kwargs['blocksize'] = self.blocksize
        hit_headers = ('pracc','prdesc','pepid','start','end','leftaa','rightaa','deltamass','substitutions','prdefline')
        obssubs = defaultdict(int)
        def hit_generator():
            for offset in range(0,npep,self.blocksize):
                pepidmap = dict((p,i) for i,p in self.protdb.get_peptides(count=self.blocksize,offset=offset))
                # print '\n'.join(pepidmap.keys()[:20])
                for r in pscmd.run(list(pepidmap.keys()),**kwargs):
                    # print >>sys.stderr, r
                    pracc,prdesc = self.accfn(r.seqdescr)
                    if not pracc:
                        continue
                    if len(r.substitutions) > 1:
                        continue
                    pep = r.peptide
                    if not self.trans:
                        start = r.start
                        end = r.end
                    else:
                        if r.nucdir == 'F':
                            start = r.nucstart
                            end = r.nucend
                        else:
                            start = r.nucend
                            end = r.nucstart
                    leftaa = r.leftaa
                    rightaa = r.rightaa
                    deltamass = r.deltamass
                    if len(r.substitutions) > 0:
                        sub = r.substitutions[0]
                        subs = "%s%d->%s(%d)"%(sub.to,sub.pos,sub.frm,sub.mut)
                        obssubs[(sub.to,sub.frm,int(round(deltamass)))] += 1
                    else:
                        subs = ""
                    yield dict(list(zip(hit_headers,(pracc,prdesc,pepidmap[pep],start,end,leftaa,rightaa,deltamass,subs,r.seqdescr))))
        oochits = OutOfCoreTable(rows=hit_generator(),headers=hit_headers)
        # print >>sys.stderr, '\n'.join(map(lambda t: "%s->%s (%+d): %d"%(t[0][0],t[0][1],t[0][2],t[1]), sorted(obssubs.items(),key=itemgetter(1))))
        praccs = OutOfCoreDistinct(map(lambda r: (r.get('pracc'),),oochits))
        self.protdb.load_proteins(map(itemgetter(0),praccs))
        proidmap = self.protdb.protein_idmap()

        hititer = map(lambda r: (proidmap[r['pracc']],r['pepid'],r['start'],
                                            r['end'],r['leftaa'],r['rightaa'],r['deltamass'],r['substitutions']), oochits)
        self.protdb.load_hits(hititer)

        prots = OutOfCoreDistinct(map(lambda r: (r.get('pracc'),r.get('prdesc'),r.get('prdefline')),oochits))
        lastacc = None
        for r in prots:
            if r[0] == lastacc:
                continue
            self.protdb.set_protein_data(proidmap[r[0]],description=r[1],defline=r[2])
            lastacc = r[0]
        for d,s in fasta(self.seqdb):
            pracc,prdesc = self.accfn(d)
            if pracc in proidmap:
                self.protdb.update_protein_data(proidmap[pracc],length=len(s))
        self.protdb.commit()

    def proteins(self,peptide):
        for h in self.protdb.get_peptide_hits(peptide):
            seqdbid,pracc,prdata = self.protdb.get_protein_data(h[0])
            yield pracc,prdata['description'],h[3],h[1],h[2],h[4],prdata['defline'],h[5],h[6],prdata['length']

if __name__ == '__main__':
    import sys

    peptides = """
AAAEGAPCVGR
AAAFEQLQK
AAALEFLNR
AAFACPTSCK
AAFVQLTQK
AASVGSCSVPK
ACRPPGSPGR
ADTIDTIEK
ADTVIVGTVK
ADVDHAIEK
AEHVTDVVK
AEPASPDSPK
AEQEASLQK
AEQLSNTLK
AFLDAQNPK
AFTCLNSLK
AFTCSSSFR
AGAESILANR
AGFAVNFFK
AGIPVYAWK
AGTPVSELTK
AHAFAVQQK
AHLLLPTLK
AIDEPNNCK
AILTESENK
ALDETDSPR
ALENLLPTK
ALIEILATR
ALIEMVVAR
ALPLGDFDR
ALRPDLADK
ALSSQHQAR
ALVQCSSHR
AMISSPPFR
ANFSGMSER
ANLTVVLLR
APATTWPVR
APEIMLNSK
AQAEVEGLGK
AQAGELWVK
ASDAVQMQR
ASGPACASPSR
ASLLCHLSR
ASPLYASYK
ASPPLFQSR
ASQGPASHSR
ASSELFSQK
ASTLLAHQR
ATILISHSR
ATTNSFHVK
ATVNSQEQK
AVDYSTPPR
AVEEDGHLK
AVGLLTVISK
AVLVDMEPK
CAECGQSFR
CATDPNSYK
CCITAAPYR
CFNSTVSSR
CGAFAGPGAPR
CGQHGEPFK
CHIAVGSPGR
CLNTLGSYK
CNCHLLGLK
CRPLTPEGK
CSISCVPPAK
CSSWDGCSR
CVVVGDGAVGK
DCNPAKPEK
DCPDPCIGW
DDVLCSSMK
DDYCEAGTK
DFLPLVLGK
DGLLTVNLR
DGTATPAPIR
DIGLSVTHR
DLAAILEEK
DLEALVTDK
DLVGYSSTR
DNTAQQISK
DPLADLNIK
DPSVTQVTR
DQEALEAVK
DSLLKPGLR
DSPINANLR
DSPMNPAHK
DTLTQLNAK
DVAAHFLAR
DVAEVDTVR
DVGSLASAQR
DVPPADQEK
EALVSQLSR
EANYIGSDK
EAVTGNGIGGK
EDGILVLSR
EDKPVTGPR
EDLCNHSGK
EDPVPNGLR
EDSNSTESK
EDTVVATLR
""".split()

    r = PeptideRemapper(peptides, sys.argv[1], verbose=True, blocksize=50)
