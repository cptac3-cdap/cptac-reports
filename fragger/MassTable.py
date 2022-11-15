#!/usr/bin/env python
import sys
import math

class MassTable:
    def sum(self,seq):
        return sum(map(self.w.get,seq))
    def aa(self,s):
        return self.w[s]
    def nterm(self):
        return self.w['n']
    def cterm(self):
        return self.w['c']
    def poscharge(self):
        return self.w['+']
    def mod(self,s,mod):
        return self.w[s+':'+mod]
    def chaa(self,s,mod):
        self.w[s] = self.mod(s,mod)
    def addaa(self,s,mass):
        self.w[s] = mass
    def addmod(self,s,mod,mass,relative=True):
        self.w[s+':'+mod] = mass
        if relative:
            self.w[s+':'+mod] += self.w[s]
    def valid_aa(self,s):
        return s in self.w
    def valid_mod(self,s,mod):
        return s+':'+mod in self.w
    def addmods(self,pep):
        if type(pep) == str:
            pep = self.str2pep(pep)
        for aa in pep:
            if aa in self.w:
                continue
            a,m = aa.split(':')
            if ":"+m in self.w:
                self.addmod(a,m,self.w[":"+m])
    @staticmethod
    def str2pep(pepstr):
        if not pepstr.startswith('n'):
            pepstr = 'n' + pepstr
        if not pepstr.endswith('c'):
            pepstr = pepstr + 'c'
        aalist = []
        pos = 0
        while pos < len(pepstr):
            if pos + 1 < len(pepstr) and pepstr[pos+1] == ':':
                aalist.append(pepstr[pos:(pos+3)])
                pos += 3
            else:
                aalist.append(pepstr[pos])
                pos += 1
        return aalist
    def peptide_mz(self,pep,z=1):
        return self.mz(self.peptide_mr(pep),z)
    def peptide_mr(self,pep):
        if type(pep) == str:
            pep = self.str2pep(pep)
        m = 0;
        for a in pep:
            m += self.aa(a)
        return m
    def peptide_mhp(self,pep):
        return self.peptide_mz(pep,1)
    def mzadduct(self,mr,add='+',z=1):
        return (mr+z*self.w[add])/z
    def mz(self,mr,z=1):
        return (mr+z*self.poscharge())/z
    def mr(self,mz,z=1):
        return (mz-self.poscharge())*z
    def mhp(self,m,z=None):
        if z == None:
            return self.mz(m,1)
        else:
            return self.mr(m,z)+self.poscharge()
    def aion(self):
        return self.w['-aion-']
    def bion(self):
        return self.w['-bion-']
    def yion(self):
        return self.w['-yion-']
    def immonium_ion(self,aa):
        m = self.aa(aa)-self.aion()-self.yion()
        return self.mz(m)
    def bions(self,pep,z=1):
        if type(pep) == str:
            pep = self.str2pep(pep)
        ions = []
        m = 0;
        for i,a in enumerate(pep[:-2]):
            m += self.aa(a)
            if i != 0:
                ions.append(self.mz(m,z))
        return ions
    def yions(self,pep,z=1):
        if type(pep) == str:
            pep = self.str2pep(pep)
        ions = []
        m = 0;
        for i,a in enumerate(reversed(pep[2:])):
            m += self.aa(a)
            if i != 0:
                ions.append(self.mz(m,z))
        return ions

class MonoisotopicMassTable(MassTable):
    def __init__(self):
        self.w = {\
            'A': 71.037113848,\
            'B': 114.53494,\
            'C': 103.009185648,\
            'D': 115.026943128,\
            'E': 129.042593208,\
            'F': 147.068414008,\
            'G': 57.021463768,\
            'H': 137.058911944,\
            'I': 113.084064088,\
            'K': 128.094963136,\
            'L': 113.084064088,\
            'M': 131.040485808,\
            'N': 114.042927536,\
            'P': 97.052763928,\
            'Q': 128.058577616,\
            'R': 156.101111152,\
            'S': 87.032028488,\
            'T': 101.047678568,\
            'U': 150.95363,\
            'V': 99.068414008,\
            'W': 186.079313056,\
            'X': 111.0,\
            'Y': 163.063328648,\
            'Z': 128.55059,\
            'n': 0.000000,\
            'c': 18.010560,\
            '+': 1.007825040,\
            'Na+': 22.989768,\
            '$': 0.00000000, \
            'C:m': 160.043225648,\
            'C:x': 161.014665648,\
            'M:o': 147.035395808,\
            'S:p': 87.032028488+79.966331,\
            'T:p': 101.047678568+79.966331,\
            'Y:p': 163.063328648+79.966331,\
            'n:a': 0.000000+42.01056,\
            ':o': 15.9994,\
            ':i': -131.1926,\
            ':m': 14.0266,\
            ':p': 79.9799,\
            ':e': 14.015650,\
            ':d': 28.031300,\
            ':j': 0.984016,\
            ':f': 27.994915,\
            ':a': 42.01056,\
            ':c': -57.021,\
            ':s1': -1.007825,\
            ':s2': 2*-1.007825,\
            '-aion-': 15.9994+12.0000,\
            '-bion-': 0.0,\
            '-yion-': 0.0,
            }

class AverageMassTable(MassTable):
    def __init__(self):
        self.w = {\
        'G':57.0519,\
        'A':71.0788,\
        'S':87.0782,\
        'P':97.1167,\
        'V':99.1326,\
        'T':101.1051,\
        'C':103.1388,\
        'I':113.1594,\
        'L':113.1594,\
        'N':114.1038,\
        'D':115.0886,\
        'Q':128.1307,\
        'K':128.1741,\
        'E':129.1155,\
        'M':131.1926,\
        'H':137.1411,\
        'F':147.1766,\
        'R':156.1875,\
        'Y':163.1760,\
        'W':186.2133,\
        'X':111.0,\
        'B':114.595,\
        'Z':128.6216,\
        'U':150.0379,\
        'C:m': 160.1942,\
        'C:x': 161.1790,\
        'M:o': 147.1955,\
        ':o': 15.9994,\
        ':i': -131.1926,\
        ':m': 14.0266,\
        ':p': 79.9799,\
        ':e': 14.0266,\
        ':d': 28.0532,\
        ':j': 0.9848,\
        ':f': 27.99491,\
        ':a': 42.01056,\
        ':c': -57.021,\
        'n': 0.000000,\
        'c':18.02   ,\
        '+':1.008,\
        '-':-1.008,\
        '$': 0.00000000,
        }

if __name__ == '__main__':
    import sys,os

    mt = MonoisotopicMassTable()
    print(mt.peptide_mz(sys.argv[1],1),mt.peptide_mz(sys.argv[1],2),mt.peptide_mz(sys.argv[1],3))
