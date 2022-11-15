
from __future__ import print_function

import sys,os,re

class fastalib_filter:
    def keep_(self,d,s):
        return self.keep(d,s)
    def keep(self,d,s):
        assert 0, 'bool keep(defline,sequence) not defined in ' + __name__ + '.'

class valid_aa(fastalib_filter):
    def __init__(self,validaa='ACDEFGHIKLMNPQRSTVWY'):
        self.nonvalid = re.compile('[^%s]'%(validaa,))
    def keep(self,d,s):
        return (len(self.nonvalid.findall(s))==0)

class all(fastalib_filter):
    def keep(self,d,s):
        return 1;

class none(fastalib_filter):
    def keep(self,d,s):
        return 0

class fastalib_writer:
    def write_(self,ofh,d,s):
        return self.write(ofh,d,s)
    def write(self,ofh,d,s):
        assert 0, 'write(defline,sequence) not defined in ' + __name__ + '.'
    def start(self,ofh):
        pass
    def finish(self,ofh):
        return True

class fasta_entry(fastalib_writer):
    def write(self,ofh,d,s):
        print('>%s\n%s' % (d,s), file=ofh)

def open_read(fname):
    bits = fname.split('.')
    if (bits[-1] == 'gz'):
        handle = os.popen('gunzip -c %s' % fname, 'r')
    elif (not os.path.isfile(fname) and os.path.isfile(fname+".gz")):
        handle = os.popen('gunzip -c %s.gz' % fname, 'r')
    else:
        handle = open(fname,'r')
    return handle

def open_write(fname,addgz=0):
    bits = fname.split('.')
    # print bits
    if (bits[-1] == 'gz'):
        handle = os.popen('gzip -c > %s' % fname, 'w')
    elif (addgz):
        handle = os.popen('gzip -c > %s.gz' % fname, 'w')
    else:
        if not os.path.exists(fname):
            handle = open(fname,'w')
        else:
            handle = open(fname,'a')
    return handle

class fasta:

    def __init__(self,input_file_name='-',output_file_name='-'):

        if input_file_name == '-':
            self.input = sys.stdin
        else:
            self.input = open_read(input_file_name)

        if output_file_name == '-':
            self.output = sys.stdout
        else:
            self.output = open_write(output_file_name)

        self.filter = all()
        self.writer = fasta_entry()

    def set_filter(self,f):
        self.filter = f

    def set_writer(self,w):
        self.writer = w

    def __iter__(self):
        return next(self)

    def __next__(self):
        firsttime = True
        nextchunkseq = True
        try:
            self.input.seek(0)
        except:
            pass
        import re
        recre = re.compile(r'^>',re.MULTILINE)
        wsre = re.compile(r'\s+',re.MULTILINE)
        defline = ""
        sequence = ""
        bufsize = 16*1024*1024
        buf = self.input.read(bufsize)
        while buf != '':
            records = recre.split(buf)
            # print records
            if nextchunkseq:
                sequence += wsre.sub('',records[0]);
            else:
                dlns = records[0].split('\n',1)
                defline += dlns[0].rstrip();
                sequence = wsre.sub('',dlns[1]);
                nextchunkseq = True
            for r in records[1:]:
                if not firsttime and self.filter.keep_(defline,sequence):
                    yield(defline,sequence)
                firsttime = False
                dlns = r.split('\n',1)
                if len(dlns) < 2:
                    nextchunkseq = False
                    defline = dlns[0].lstrip();
                else:
                    defline = dlns[0].strip();
                    sequence = wsre.sub('',dlns[1]);
            buf = self.input.read(bufsize)
        if not firsttime and self.filter.keep(defline,sequence):
            yield (defline,sequence)

    def process(self):
        firsttime = True
        nextchunkseq = True
        try:
            self.input.seek(0)
        except:
            pass
        import re
        recre = re.compile(r'^>',re.MULTILINE)
        wsre = re.compile(r'\s+',re.MULTILINE)
        defline = ""
        sequence = ""
        self.writer.start(self.output)
        bufsize = 16*1024*1024
        buf = self.input.read(bufsize)
        while buf != '':
            records = recre.split(buf)
            # print records
            if nextchunkseq:
                sequence += wsre.sub('',records[0]);
            else:
                dlns = records[0].split('\n',1)
                defline += dlns[0].rstrip();
                sequence = wsre.sub('',dlns[1]);
                nextchunkseq = True
            for r in records[1:]:
                if not firsttime and self.filter.keep_(defline,sequence):
                    try:
                        self.writer.write_(self.output,defline,sequence)
                    except IOError:
                        sys.exit(0)
                firsttime = False
                dlns = r.split('\n',1)
                if len(dlns) < 2:
                    nextchunkseq = False
                    defline = dlns[0].lstrip();
                else:
                    defline = dlns[0].strip();
                    sequence = wsre.sub('',dlns[1]);
            buf = self.input.read(bufsize)
        if not firsttime and self.filter.keep(defline,sequence):
            try:
                self.writer.write_(self.output,defline,sequence)
            except IOError:
                sys.exit(0)
        return self.writer.finish(self.output)
