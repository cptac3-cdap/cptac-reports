#!/bin/env python

import re, math, csv
import sys
from .SearchResult import *
from .fopen import *
from .Scan import Scan

class DIATSVReader:
    hmass = 1.007825040
    def __init__(self,f,alignments=True):
        self.alignments = alignments
        self.handle = None
        self.f = f
        self.any = True
        self.reset()

    def reset(self):
        if self.any:
            if self.handle:
                self.handle.close()
                self.handle = None
            self.handle = fopen(self.f,'r')
        self.any = False

    def __iter__(self):
        return self.next()

    def next(self):
        self.any = True
        dr = csv.DictReader(self.handle,dialect='excel-tab')
        for ii,r in enumerate(dr):
            pep = r['Peptide']
            splpep = re.split(r'\[([+-]?\d*\.?\d*)\]',pep)
            # print pep,splpep
            pepseq = "".join(splpep[::2])
            currseq = splpep[0]
            mods = []
            for i in range(0,len(splpep)-2,2):
                currseq += splpep[i]
                mods.append((currseq[-1],len(currseq),float(splpep[i+1])))

            nfr = r['numFragments']
            pro = r['Protein']
            for k in dr.fieldnames[3:]:
                if float(r[k]) == 0:
                    continue
                basename = k.rsplit('.',1)[0]
                sr = SearchResult()
                sr.set('basename',basename)
                sr.set('base',basename)
                sr.set('start_scan',Scan(ii+1))
                sr.set('end_scan',Scan(ii+1))
                psm = PeptideSpectrumMatch()
                psm.set('peptide',pepseq)
                if len(mods) > 0:
                    for m in mods:
                        psm.append('mods', '%s%d:%+.3f'%m)
                else:
                    psm.set('mods','-')
                psm.sort('mods',value=lambda s: (int(s.split(':')[0][1:]),float(s.split(':')[1])))
                psm.set('rank',1)
                psm.set('numFragments',int(nfr))
                psm.set('abundance',float(r[k]))
                psm.set('estfdr',0.0)
                if self.alignments:
                    for proi in pro.split(';'):
                        psm.append('protein',proi.strip('|'))
                sr.set('scores',",".join(['numFragments','abundance']))
                sr.set('parameters',"")
                sr.addPSM(psm)
                yield sr
