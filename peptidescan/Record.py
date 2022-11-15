import re

class Substitution(object):
    def __init__(self,*args):
        if len(args) == 4:
            self.frm  = args[0]
            self.pos  = int(args[1])
            self.to   = args[2]
            self.mut  = int(args[3])

class SAV(object):
    def __init__(self,*args):
        self.frm = args[0]
        self.pos = int(args[1])
        self.to = args[2]
        self.dbsnp = None
        self.maf = None
        if len(args) > 3:
            self.dbsnp = args[3]
        if len(args) > 4:
            try:
                self.dbsnp = float(args[4])
            except ValueError:
                self.dbsnp = float(args[4].rstrip('%'))/100.0

class Result(object):
    def __init__(self,*sl):
        if len(sl) == 13:
            self.pepindex = int(sl[0])
            self.start = int(sl[1])
            self.end   = int(sl[2])
            self.leftaa = sl[3]
            self.peptide = sl[4]
            self.rightaa = sl[5]
            self.nucstart = int(sl[6])
            self.nucend = int(sl[7])
            self.nucframe = int(sl[8])
            self.nucdir = sl[9]
            self.nucleotides = sl[10]
            self.seqindex = sl[11]
            self.seqdescr = sl[12]
            self.deltamass = 0.0
            self.substitutions = []
            self.savs = []
            self.query = sl[4]
            if ' /delta=' in self.seqdescr and ' /sub1=' in self.seqdescr:
                m = re.search(r' /delta=([\d.+-]+)',self.seqdescr)
                assert m
                self.deltamass = float(m.group(1))
                self.query = list(self.query)
                for m in re.finditer(r' /sub(\d+)=([A-Z])(\d+)->([A-Z])[(](\d+)[)]',self.seqdescr):
                    sub = Substitution(*[m.group(i) for i in [2,3,4,5]])
                    assert self.query[sub.pos-1] == sub.to
                    self.query[sub.pos-1] = sub.frm
                    self.substitutions.append(sub)
                self.query = ''.join(self.query)
                self.seqdescr = self.seqdescr.rsplit(' /sub1=',1)[0]

    @staticmethod
    def parse(l):
        sl = l.strip().split(None,12)
        r = Result(*sl)
        return r

    def __str__(self):
        return "\t".join(map(str,[self.pepindex,self.start,self.end,
                                 self.leftaa,self.peptide,self.rightaa,
                                 self.nucstart,self.nucend,self.nucframe,
                                 self.nucdir,self.nucleotides,
                                 self.seqindex,self.seqdescr]))
