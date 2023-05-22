#!/usr/bin/env python3

from . find_elementtree import ET

import re, math
import sys
from . SearchResult import *
from collections import defaultdict
from operator import itemgetter
from . Scan import Scan
from . pymsiosplit import splitResultFileExtn, splitResultFileMZID

def ieeefloat(s):
    if s == 'inf':
        return 1e+20
    if s == '-inf':
        return -1e+20
    if s == 'nan':
        return None
    return float(s)

def getbasename(f):
    md = splitResultFileExtn(f)
    if md:
        return md['base']
    md = splitResultFileMZID(f)
    if md:
        return md['base']
    return None

class MZIdentMLReader:
    def __init__(self,h,alignments=True,filename=None):
        self.handle = h

        self.basename = None
        if filename:
            self.basename = getbasename(filename)

        self.tagre = re.compile(r'^{(.*)}(.*)$')
        self.context = []
        self.ns = ''
        self.any = True
        self.alignments = alignments
        self.reset()

    def reset(self):
        self.any = False
        self.searchdb = {}
        self.specdata = {}
        self.proteins = {}
        self.peptides = {}
        self.pepevidence = defaultdict(list)

    def __iter__(self):
        return next(self)

    def __next__(self):
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
                    for (k,v) in list(ele.attrib.items()):
                        if k.endswith('schemaLocation'):
                            self.ns = '{%s}'%(v.split()[0],)
                continue

            if self.context[-1] == 'SearchDatabase':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                if self.alignments:
                    self.extract_searchdb(ele)

                ele.clear()

            if self.context[-1] == 'SpectraData':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                self.extract_spectradata(ele)

                ele.clear()

            if self.context[-1] == 'DBSequence':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                if self.alignments:
                    self.extract_dbsequence(ele)

                ele.clear()

            elif self.context[-1] == 'Peptide':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                self.extract_peptide(ele)

                ele.clear()

            elif self.context[-1] == 'PeptideEvidence':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                if self.alignments:
                    self.extract_pepevidence(ele)

                ele.clear()

            elif self.context[-1] == 'SpectrumIdentificationResult':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                yield self.extract_result(ele)

                ele.clear()

            self.context.pop()

    def extract_searchdb(self,ele):
        id = ele.attrib['id']
        taxa = None
        dbname = None
        for param in ele.getiterator(self.ns+'cvParam'):
            if param.attrib['name'] == 'taxonomy: common name' or \
                   param.attrib['accession'] == 'MS:1001468':
                taxa = param.attrib['value']
            if param.attrib['name'] == 'database name' or \
                   param.attrib['accession'] == 'MS:1001013':
                dbname = param.attrib['value']
        self.searchdb[id] = (dbname,taxa)

    def extract_spectradata(self,ele):
        id = ele.attrib['id']
        loc = ele.attrib.get('location')
        name = ele.attrib.get('name')
        basename = None
        if not basename and name:
            basename = getbasename(name)
        if not basename and loc:
            basename = getbasename(loc)
        if not basename:
            basename = getbasename(id)
        if not basename:
            basename = id
        self.specdata[id] = basename

    def extract_dbsequence(self,ele):
        id = ele.attrib['id']
        acc = ele.attrib.get('accession')
        if acc == None:
            acc = id
        desc = None
        for param in ele.findall(self.ns+'cvParam'):
            if param.attrib['name'] == 'protein description' or \
               param.attrib['accession'] == 'MS:1001088':
                desc = param.attrib['value']
        sdbref = ele.attrib.get('searchDatabase_ref')
        self.proteins[id] = (acc,desc,sdbref)

    def extract_peptide(self,ele):
        id = ele.attrib['id']
        seq = ele.find(self.ns+'PeptideSequence').text
        mods = []
        for mod in ele.findall(self.ns+'Modification'):
            delta = float(mod.attrib['monoisotopicMassDelta'])
            pos = int(mod.attrib['location'])
            if pos == 0:
                aa = '['
            elif pos == (len(seq)+1):
                aa = ']'
            else:
                aa = seq[pos-1]
            mods.append((aa,pos,delta))
        mods.sort(key=itemgetter(1,2,0))
        self.peptides[id] = (seq,mods)

    def extract_pepevidence(self,ele):
        id = ele.attrib['id']
        laa = ele.attrib.get('pre')
        raa = ele.attrib.get('post')
        try:
            st = int(ele.attrib['start'])
        except KeyError:
            st = None
        try:
            ed = int(ele.attrib['end'])
        except KeyError:
            ed = None
        dbseqref = ele.attrib['dBSequence_ref']
        pepref = ele.attrib['peptide_ref']
        assert pepref in self.peptides
        self.pepevidence[pepref].append((self.proteins[dbseqref],laa,st,ed,raa))

    def extract_result(self,ele):
        sr = SearchResult()
        specid = ele.attrib['spectrumID']
        if len(self.specdata) > 1 and ele.attrib['spectraData_ref'] in self.specdata:
            basename = self.specdata[ele.attrib['spectraData_ref']]
            sr.set('basename',basename)
        else:
            assert(self.basename)
            sr.set('basename',self.basename)
        index = None
        scan = None
        title = None
        rt = None
        precint = None
        m = re.search(r'^index=(\d+)$',specid)
        if m:
            index = m.group(1)
            sr.set('index',int(index))
        m = re.search(r'^controllerType=(\d+) controllerNumber=(\d+) scan=(\d+)$',specid)
        if m:
            # scan = '.'.join(map(str,[m.group(1),m.group(2),m.group(3)]))
            scan = int(m.group(3))
        # Put the others here...
        for cvelt in ele.findall(self.ns+'cvParam'):
            if cvelt.attrib['accession'] == 'MS:1000796' or \
               cvelt.attrib['name'] == 'spectrum title':
                title = cvelt.attrib['value']
            if cvelt.attrib['accession'] == 'MS:1001141' or \
               cvelt.attrib['name'] == 'intensity of precursor ion':
                precint = float(cvelt.attrib['value'])
            if cvelt.attrib['accession'] == 'MS:1000894' or \
               cvelt.attrib['name'] == 'retention time':
                rt = float(cvelt.attrib['value'])
            if cvelt.attrib['accession'] == 'MS:1001115' or \
               cvelt.attrib['name'] == 'scan number(s)':
                scan = int(cvelt.attrib['value'])
        if scan:
            sr.set('start_scan',scan)
            sr.set('end_scan',scan)
        elif title:
            m,st,ed = Scan.search(title,prefix=r'^Scan +Number:? *',flags=re.IGNORECASE)
            if m:
                sr.set('start_scan',st)
                sr.set('end_scan',ed)
            else:
                m = re.search(r'\.(\d+)\.(\d+)\.\d',title)
                if m and int(m.group(1)) <= int(m.group(2)):
                    sr.set('start_scan',Scan(m.group(1)))
                    sr.set('end_scan',Scan(m.group(2)))
                else:
                    raise RuntimeError("Can't determine scan number from spectrum title")
        elif index:
            sr.set('start_scan',index)
            sr.set('end_scan',index)

        if rt:
            sr.set('retention_time',rt)
        if precint:
            sr.set('precursor_intensity',precint)
        sii = ele.find(self.ns+'SpectrumIdentificationItem')
        precursor_mz = float(sii.attrib['experimentalMassToCharge'])
        sr.set('precursor_mz',precursor_mz)

        rankingvalues = set()
        for it,sii in enumerate(ele.findall(self.ns+'SpectrumIdentificationItem')):
            psm = PeptideSpectrumMatch()
            z = int(sii.attrib['chargeState'])
            psm.set('assumed_charge',z)
            rank = int(sii.attrib['rank'])
            psm.set('rank',rank)
            psm.set('notierank',it+1)
            psm.set('precursor_neutral_mass',precursor_mz*z-1.0078*z)
            calculatedMassToCharge = float(sii.attrib['calculatedMassToCharge'])
            psm.set('calc_neutral_pep_mass',calculatedMassToCharge*z-1.0078*z)
            psm.set('massdiff',precursor_mz*z-calculatedMassToCharge*z)
            pepref = sii.attrib['peptide_ref']
            psm.set('peptide',self.peptides[pepref][0])
            for m in self.peptides[pepref][1]:
                psm.append('mods','%s%d:%+.3f'%m)
            if psm.get('mods') == None:
                psm.set('mods','-')
            psm.sort('mods',value=lambda s: (int(s.split(':')[0][1:]),float(s.split(':')[1])))
            for pepev in self.pepevidence[pepref]:
                psm.append('protein',pepev[0][0])
                if pepev[0][1]:
                    psm.append('protein_descr',pepev[0][1],delim='\t')
                if pepev[0][2] in self.searchdb and \
                       self.searchdb[pepev[0][2]][1]:
                    psm.append('protein_org',self.searchdb[pepev[0][2]][1])
                if pepev[0][2] in self.searchdb and \
                       self.searchdb[pepev[0][2]][0]:
                    psm.append('protein_source',self.searchdb[pepev[0][2]][0])
                if pepev[1]:
                    psm.append('peptide_prev_aa',pepev[1])
                if pepev[2]:
                    psm.append('peptide_offset',pepev[2])
                if pepev[4]:
                    psm.append('peptide_next_aa',pepev[4])

            params = []
            for scoreelt in (sii.findall(self.ns+'cvParam')+sii.findall(self.ns+'userParam')):
                sname = scoreelt.attrib['name']
                try:
                    val = scoreelt.attrib['value']
                    val = float(scoreelt.attrib['value'])
                    val = int(scoreelt.attrib['value'])
                except:
                    pass
                if val != "" and val != None:
                    psm.set(sname, val)
                    params.append(sname)
            if psm.has('MS-GF:SpecEValue'):
                psm.set('eval',psm.get('MS-GF:SpecEValue'))
                params.append('eval')
                ranking_key = 'eval'
            if psm.has('IDPicker:MinQValue'):
                psm.set('estfdr',psm.get('IDPicker:MinQValue'))
                params.append('estfdr')
                ranking_key = 'estfdr'
            elif psm.has('MS-GF:QValue'):
                psm.set('estfdr',psm.get('MS-GF:QValue'))
                params.append('estfdr')
                ranking_key = 'estfdr'
            if psm.has('MS-GF:PepQValue'):
                psm.set('pepfdr',psm.get('MS-GF:PepQValue'))
                params.append('pepfdr')
            if ranking_key == 'eval':
                rankingvalues.add("%.5e"%psm.get('eval'))
                psm.set('trank',len(rankingvalues))
            elif ranking_key == 'estfdr':
                rankingvalues.add("%.5e"%psm.get('estfdr'))
                psm.set('trank',len(rankingvalues))
            # print psm
            if not sr.has('precursor_area') and psm.has('CPTAC-CDAP:PrecursorArea'):
                sr.set('precursor_area',psm.get('CPTAC-CDAP:PrecursorArea'))
            if not sr.has('precursor_area') and psm.has('CPTAC-CAP:PrecursorArea'):
                sr.set('precursor_area',psm.get('CPTAC-CAP:PrecursorArea'))

            if psm.get('CPTAC-CDAP:FullyLocalized') == "Y":
                # Check the localization is correct...
                phospos1 = set()
                if psm.get('mods') != "-":
                    for m in psm.get('mods').split(','):
                        aapos,delstr = m.split(':')
                        aa = aapos[0]
                        pos = int(aapos[1:])
                        if abs(float(delstr)-80) < 0.1:
                            phospos1.add(pos)
                phospos2 = set()
                curpos = 0
                for chunk in psm.get('CPTAC-CDAP:PhosphoRSPeptide').split(']')[:-1]:
                    aas,prob = chunk.split('[')
                    curpos += len(aas)
                    if float(prob) >= 99:
                        phospos2.add(curpos)
                if phospos1 != phospos2:
                    psm.set('CPTAC-CDAP:FullyLocalized',"N")

            iTRAQkeys = """
                iTRAQ114 iTRAQ115 iTRAQ116 iTRAQ117 iTRAQ113 iTRAQ118 iTRAQ119 iTRAQ121
                iTRAQ4-114 iTRAQ4-115 iTRAQ4-116 iTRAQ4-117
                iTRAQ8-114 iTRAQ8-115 iTRAQ8-116 iTRAQ8-117
                iTRAQ8-113 iTRAQ8-118 iTRAQ8-119 iTRAQ8-121
            """.split()
            TMTkeys = """
                TMT18-126C TMT18-127C TMT18-127N TMT18-128C TMT18-128N TMT18-129C TMT18-129N TMT18-130C
                TMT18-130N TMT18-131N TMT18-131C TMT18-132N TMT18-132C TMT18-133N TMT18-133C TMT18-134N TMT18-134C TMT18-135N
                TMT16-126C TMT16-127C TMT16-127N TMT16-128C TMT16-128N TMT16-129C TMT16-129N TMT16-130C
                TMT16-130N TMT16-131N TMT16-131C TMT16-132N TMT16-132C TMT16-133N TMT16-133C TMT16-134N
                TMTXX-134C TMTXX-135N

                TMT11-126 TMT11-126C TMT11-127C TMT11-127N TMT11-128C TMT11-128N TMT11-129C TMT11-129N TMT11-130C
                TMT11-130N TMT11-131N TMT11-131C TMT11-131
                TMTXX-132N TMTXX-132C

                TMT10-126 TMT10-127C TMT10-127N TMT10-128C TMT10-128N TMT10-129C TMT10-129N TMT10-130C TMT10-130N TMT10-131
                TMTXX-131C TMTXX-132N

                TMT6-126 TMT6-127 TMT6-128 TMT6-129 TMT6-130 TMT6-131
                TMT2-126 TMT2-127

                TMT-127C TMT-127N TMT-128C TMT-128N TMT-129C TMT-129N TMT-130C TMT-130N
                TMT-126 TMT-127 TMT-128 TMT-129 TMT-130 TMT-131
                TMT-131N TMT-131C
            """.split()
            # print "\n----",
            for k in iTRAQkeys + TMTkeys:
                for prefix in ("CPTAC-CDAP:","CPTAC-CAP:",""):
                    value = psm.get(prefix+k)
                    if value != None:
                    # print "\n",prefix,k,value,
                        if k.startswith('iTRAQ'):
                            if '-' not in k:
                                label = "iTRAQ-"+k[-3:]
                            else:
                                label = k
                            k1 = "%s_reporter_ion_intensity"%(label,)
                            k2 = "%s_reporter_ion_massdelta"%(label,)
                        elif k.startswith('TMT'):
                            k1 = "%s_reporter_ion_intensity"%(k,)
                            k2 = "%s_reporter_ion_massdelta"%(k,)
                        if sr.has(k1):
                            continue
                        if not isinstance(value,str):
                            sr.set(k1,value)
                            continue
                        v1,v2 = value.split('/')
                        # print v1,v2,
                        v1 = float(v1)
                        # print k1,v1,
                        sr.set(k1,v1)
                        try:
                            v2 = float(v2)+0.0 # avoid -0.0!
                            # print k2,v2,
                            sr.set(k2,v2)
                        except ValueError:
                            pass
            psm.set('metrics',params)
            sr.addPSM(psm)
        sr.set('ranking_key',ranking_key)
        sr.sortPSMs()
        # print sr
        return sr
