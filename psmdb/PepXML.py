#!/bin/env python

from . find_elementtree import ET
from xml.sax.saxutils import escape
import re, math, operator
from . SearchResult import *
import sys
from fragger.MassTable import MonoisotopicMassTable as MT

class PepXMLReader(object):
    def __init__(self,h):
        self.handle = h
        self.tagre = re.compile(r'^{(.*)}(.*)$')
        self.context = []
        self.ns = ''
        self.typemap = {}
        self.ranking_key = 'hit_rank'
        self.format = 'pepxml'
        self.analysis_results = []
        self.mt = MT()
        for tag in ('start_scan',
                    'end_scan',
                    'assumed_charge',
                    'index',
                    'hit_rank',
                    'peptide_offset',
                    'num_tot_proteins',
                    'num_matched_ions',
                    'tot_num_ions',
                    'num_tol_term',
                    'num_missed_cleavages',
                    ):
            self.typemap[tag] = int
        for tag in ('precursor_neutral_mass',
                    'retention_time_sec',
                    'calc_neutral_pep_mass',
                    'massdiff',
                    ):
            self.typemap[tag] = float
        self.typemap['decoy'] = lambda s: bool(int(s))

    def make_path(self,*l):
        return '/'.join(map(lambda t: self.ns+t,l))

    def __iter__(self):
        return self.next()

    def next(self):
        for (event, ele) in ET.iterparse(self.handle,('start','end')):
            # print event,self.context,ele.tag,ele.attrib

            if event == 'start':
                m = self.tagre.search(ele.tag)
                if m:
                    self.context.append(m.group(2))
                else:
                    self.context.append(ele.tag)
                if len(self.context) == 1:
                    for (k,v) in ele.attrib.iteritems():
                        if k.endswith('schemaLocation'):
                            self.ns = '{%s}'%(v.split()[0],)
                continue

            if self.context[-1] == 'spectrum_query':

                # ele is now a DOM like object with everything from
                # the highest level scan...

                yield self.extract_query(ele)

                ele.clear()

            self.context.pop()

    def extract_query(self,ele):
        sr = SearchResult()
        for (k,v) in ele.attrib.iteritems():
            sr.set(k,self.typemap.get(k,str)(v))
        z = sr.get('assumed_charge')
        sr.unset('assumed_charge')
        sr.unset('index')
        precursor_mass = sr.get('precursor_neutral_mass')
        sr.unset('precursor_neutral_mass')
        if sr.has('retention_time_sec'):
            sr.rename('retention_time_sec','retention_time')

        if sr.get('spectrum'):
            t=sr.get('spectrum')
            m = re.search(r'^(.*)\.(\d+)\.(\d+)\.\d',t)
            if m and int(m.group(2)) <= int(m.group(3)):
                sr.set('start_scan',int(m.group(2)))
                sr.set('end_scan',int(m.group(3)))
                sr.set('basename',m.group(1))
            else:
                raise "Can't determine scan number from spectrum field value"

        # print ele.tag,ele.attrib

        for hitelt in ele.findall(self.make_path('search_result','search_hit')):
            psm = PeptideSpectrumMatch()
            for (k,v) in hitelt.attrib.iteritems():
                psm.set(k,self.typemap.get(k,str)(v))
            psm.set('assumed_charge',z)
            psm.set('precursor_neutral_mass',precursor_mass)
            for ap in hitelt.findall(self.make_path('alternative_protein')):
                for key in ('protein','peptide_offset','peptide_prev_aa','peptide_next_aa','num_tol_term'):
                    if ap.attrib.has_key(key):
                        psm.append(key,self.typemap.get(key,str)(ap.attrib[key]))
                for key in ('protein_descr',):
                    if ap.attrib.has_key(key):
                        psm.append(key,self.typemap.get(key,str)(ap.attrib[key]),delim='\t')

            for modelt in hitelt.findall(self.make_path('modification_info','mod_aminoacid_mass')):
                pos = int(modelt.attrib['position'])
                mass = float(modelt.attrib['mass'])
                sym = psm.get('peptide')[pos-1]
                massdiff = float(modelt.attrib['mass']) - self.mt.aa(sym)
                psm.append('mods','%s%d:%+.3f'%(sym,pos,math.floor(1000*massdiff)/1000))
            for modelt in hitelt.findall(self.make_path('modification_info')):
                if modelt.attrib.has_key('mod_nterm_mass'):
                    massdiff = float(modelt.attrib['mod_nterm_mass'])
                    psm.append('mods','%s%d:%+.3f'%('[',0,math.floor(1000*massdiff)/1000))
            if psm.get('mods') == None:
                psm.set('mods','-')
            psm.sort('mods',value=lambda s: (int(s.split(':')[0][1:]),float(s.split(':')[1])))
            for scoreelt in hitelt.findall(self.make_path('search_score')):
                psm.set(scoreelt.attrib['name'],self.typemap.get(scoreelt.attrib['name'],str)(scoreelt.attrib['value']))
            for paramelt in hitelt.findall(self.make_path('parameter')):
                psm.set(paramelt.attrib['name'],self.typemap.get(paramelt.attrib['name'],str)(paramelt.attrib['value']))
            for analtag,analkey in self.analysis_results:
                for analelt in hitelt.getiterator(self.make_path(analtag)):
                    psm.set(analtag+'.'+analkey,self.typemap.get(analtag+'.'+analkey,str)(analelt.attrib[analkey]))
            if psm.has('expect'):
                psm.rename('expect','eval')
            sr.addPSM(psm)
        sr.set('ranking_key',self.ranking_key)
        sr.sortPSMs()
        return sr

class PepXMLRewriter:
    def __init__(self):
        pass

    def rewrite(self,input,output):
        handle = open(input,'r')
        document = ET.parse(handle)
        handle.close()
        root = document.getroot()
        for k in root.attrib:
            if k.endswith('schemaLocation'):
                xmlns = root.attrib[k].split()[0]
                break
        self.ns = '{'+xmlns+'}'
        for ele in root.getiterator(self.ns+'search_hit'):
            # print >>sys.stderr, "HERE0"
            self.fixhit(ele)
        for ele in root.getiterator(self.ns+'search_result'):
            # print >>sys.stderr, "HERE1"
            self.fixresult(ele)
        for ele in root.getiterator(self.ns+'spectrum_query'):
            # print >>sys.stderr, "HERE2"
            self.fixspectrum(ele)
        PepXMLRewriter.cleanup(root,lambda e: e.tag)

        self.out = open(output,'wb')
        self.out.write('<?xml version="1.0" encoding="ISO-8859-1"?>\n')
        document.write(self)
        self.out.write('\n')
        self.out.close()

    def write(self,s):
        # print >>sys.stderr, "Writing out: '%s'\n"%s
        if 'xmlns:ns0' in s:
            s = s.replace(r'xmlns:ns0','xmlns')
        if '<ns0:' in s:
            s = s.replace('<ns0:','<')
        if '</ns0:' in s:
            s = s.replace('</ns0:','</')
        if "'" in s:
            s = s.replace("'","&apos;")
        self.out.write(s)

    @staticmethod
    def cleanup(elem,filter):
        out = []
        for e in elem:
            PepXMLRewriter.cleanup(e, filter)
            if filter(e):
                out.append(e)
        elem[:] = out

    def fixhit(self,ele):
        return

    def fixresult(self,ele):
        if not any(sh.tag for sh in ele.findall(self.ns+'search_hit')):
            ele.tag = None
        return

    def fixspectrum(self,ele):
        if not any(sr.tag for sr in ele.findall(self.ns+'search_result')):
            ele.tag = None
        return

class ParsimonyPepXML(PepXMLRewriter):
    def __init__(self,parsimony):
        self.parsimony = set(parsimony)
        # print >>sys.stderr, self.parsimony
        self.prkeys = "protein protein_descr num_tol_term peptide_prev_aa peptide_next_aa".split()

    def fixhit(self,ele):

        # for e in ele:
        #     print >>sys.stderr, e.tag, e.text, e.tail

        nprot = int(ele.attrib['num_tot_proteins'])
        proteins = [dict([(k,ele.attrib.get(k)) for k in self.prkeys])]
        for apr in ele.getiterator(self.ns+'alternative_protein'):
            proteins.append(dict([(k,apr.attrib.get(k)) for k in self.prkeys]))
            apr.tag = None
            headtext = apr.text
            tailtext = apr.tail
        assert (len(proteins) == nprot), "%s,%s: %d != %d"%(ele.attrib['peptide'],ele.attrib['protein'],len(proteins),nprot)
        # print >>sys.stderr, proteins
        proteins = filter(lambda d: d['protein'] in self.parsimony,proteins)
        # print >>sys.stderr, proteins

        if len(proteins) == 0:
            ele.tag = None
            return

        ele.attrib['num_tot_proteins'] = str(len(proteins))
        pr = proteins[0]
        for k in self.prkeys:
            if pr.get(k):
                ele.attrib[k] = pr[k]
            else:
                if k in ele.attrib:
                    del ele.attrib[k]

        for i,apr in enumerate(proteins[1:]):
            apr_ele = ET.Element(self.ns+"alternative_protein")
            apr_ele.text = headtext
            apr_ele.tail = tailtext
            for k in self.prkeys:
                if apr.get(k):
                    apr_ele.attrib[k] = apr[k]
            ele.insert(i,apr_ele)
