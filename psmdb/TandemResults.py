#!/usr/bin/env python3

from .find_elementtree import ET
import re, math
import sys
from .SearchResult import *
from .fopen import *
from .Scan import Scan
from .pymsiosplit import splitResultFileExtn

class TandemResultsReader:
    format = 'tandem'
    hmass = 1.007825040
    def __init__(self,f):
        self.handle = None
        self.f = f
        self.md = splitResultFileExtn(f)
        if not self.md or self.format != self.md['engine']:
            raise RuntimeError("Bad filename...")
        self.md['basename'] = self.md['base']

        self.tagre = re.compile(r'^{(.*)}(.*)$')
        self.context = []
        self.ns = ''
        self.typemap = {}
        for tag in ('id',
                    'z',
                    'start',
                    'end',
                    'y_ions',
                    'b_ions',
                    'missed_cleavages',
                    ):
            self.typemap[tag] = int
        for tag in ('mh',
                    'expect',
                    'sumI',
                    'maxI',
                    'fI',
                    'delta',
                    'hyperscore',
                    'nextscore',
                    'y_score',
                    'b_score'
                    ):
            self.typemap[tag] = float
        self.scores = ['y_ions','b_ions','y_score','b_score','hyperscore','nextscore','expect','sumI','maxI','fI']
        self.params = []
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
        for (event, ele) in ET.iterparse(self.handle,('start','end')):
            # print event,self.context,ele.tag,ele.attrib

            if event == 'start':
                m = self.tagre.search(ele.tag)
                if m:
                    self.context.append(m.group(2))
                else:
                    self.context.append(ele.tag)
                if len(self.context ) == 1:
                    for (k,v) in ele.attrib.iteritems():
                        if k.endswith('schemaLocation'):
                            self.ns = '{%s}'%(v.split()[0],)
                continue

            if self.context[-1] == 'group' and ele.attrib['type'] == 'model':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                yield self.extract_query(ele)

                ele.clear()

            self.context.pop()

    def extract_query(self,ele):
        sr = SearchResult()
        for k,v in self.md.iteritems():
            sr.set(k,v)
        for (k,v) in ele.attrib.iteritems():
            if v != "":
                sr.set(k,self.typemap.get(k,str)(v))
        t = None
        for descele in ele.findall(self.ns+'group/note'):
            if descele.attrib['label'] == 'Description':
                t = descele.text
        if t:
            m,st,ed = Scan.search(t,prefix=r'^Scan +Number:? *',flags=re.IGNORECASE)
            if m != 0:
                sr.set('start_scan',st)
                sr.set('end_scan',ed)
            else:
                m = re.search(r'\.(\d+)\.(\d+)\.[1234]\.dta',t)
                if m and int(m.group(1)) <= int(m.group(2)):
                    sr.set('start_scan',Scan(m.group(1)))
                    sr.set('end_scan',Scan(m.group(2)))
                else:
                    sr.set('start_scan',Scan(sr.get('id')))
                    sr.set('end_scan',Scan(sr.get('id')))
        else:
            sr.set('start_scan',Scan(sr.get('id')))
            sr.set('end_scan',Scan(sr.get('id')))

        nscoredpeps = None
        for grelt in ele.findall(self.ns+'group'):
            if grelt.attrib['label'] != 'supporting data':
                continue
            for trelt in grelt.findall('{http://www.bioml.com/gaml/}trace'):
                if trelt.attrib['type'] != 'convolution survival function':
                    continue
                nscoredpeps = int(trelt.find('{http://www.bioml.com/gaml/}Ydata/{http://www.bioml.com/gaml/}values').text.split()[0])

        z = sr.get('z')
        mh = sr.get('mh')
        sr.set('precursor_mz',(mh+(z-1)*self.hmass)/z)
        sr.unset('z')
        sr.unset('mh')
        sr.unset('label')
        sr.unset('id')
        sr.unset('type')
        j = 0
        for prelt in ele.findall(self.ns+'protein'):
            pracc = prelt.attrib['label'].split()[0]
            for domelt in prelt.findall(self.ns+'peptide/domain'):
                psm = sr.findPSM('peptide',domelt.attrib['seq'])
                if psm:
                    psm.append('protein',pracc)
                    psm.append('peptide_offset',domelt.attrib['start'])
                    continue
                j += 1
                psm = PeptideSpectrumMatch()
                psm.append('protein',pracc)
                psm.append('peptide_offset',domelt.attrib['start'])
                for (k,v) in domelt.attrib.iteritems():
                    if k not in ('id','start','end'):
                        psm.set(k,self.typemap.get(k,str)(v))
                psm.set('hit_rank',j)
                psm.set('assumed_charge',z)
                psm.set('precursor_neutral_mass',mh-TandemResultsReader.hmass)
                psm.set('calc_neutral_pep_mass',psm.get('mh')-TandemResultsReader.hmass)
                psm.rename('delta','massdiff')
                psm.rename('missed_cleavages','num_missed_cleavages')
                psm.rename('pre','prev_aa')
                psm.rename('post','next_aa')
                psm.unset('mh')
                psm.set('eval',psm.get('expect'))
                if nscoredpeps:
                    psm.set('pval',psm.get('expect')/nscoredpeps)
                ntermmod = False
                for modelt in domelt.findall(self.ns+'aa'):
                    pos = int(modelt.attrib['at'])-int(domelt.attrib['start'])
                    mass = float(modelt.attrib['modified'])
                    sym = psm.get('seq')[pos]
                    if pos == 0 and not ntermmod:
                        psm.append('mods','%s%d:%+.3f'%('[',0,mass))
                        ntermmod = True
                    else:
                        psm.append('mods','%s%d:%+.3f'%(sym,pos+1,mass))
                if psm.get('mods') == None:
                    psm.set('mods','-')
                psm.sort('mods',value=lambda s: (int(s.split(':')[0][1:]),float(s.split(':')[1])))
                psm.set('peptide',psm.get('seq'))
                psm.unset('seq')
                psm.set('scores',','.join(self.scores))
                psm.set('parameters',','.join(self.params))
                sr.addPSM(psm)
        any = False
        for m in sr.PSMiter():
            if m.get('eval') == sr.get('expect'):
                any = True
                m.set('sumI',sr.get('sumI'))
                m.set('maxI',sr.get('maxI'))
                m.set('fI',sr.get('fI'))
        sr.unset('sumI')
        sr.unset('maxI')
        sr.unset('fI')
        sr.unset('expect')
        sr.set('ranking_key','-hyperscore')
        sr.sortPSMs()
        return sr


class KScoreResultsReader(TandemResultsReader):
    format = 'kscore'

class KMScoreResultsReader(TandemResultsReader):
    format = 'kmscore'

class SScoreResultsReader(TandemResultsReader):
    format = 'sscore'

class PDEScoreResultsReader(TandemResultsReader):
    format = 'pdescore'

class KMScoreResultsReader(TandemResultsReader):
    format = 'kmscore'
