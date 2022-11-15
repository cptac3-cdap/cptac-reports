
import os
import gzip
import bz2
import sys
import urllib
import tempfile
import re

__all__ = [ 'fopen', 'charCounter' ]

def fopen(f,spec,comp=None):
    if callable(f):
        f = f()
    if type(f) == str or type(f) == unicode:
        sf = f.split(':')
        if sf[0] in ('http','file','ftp','https'):
            if spec == 'r':
                (tf,x) = urllib.urlretrieve(f)
                if f.endswith('.gz') or comp in ('gz','tgz'):
                    h = handle_wrapper(gzip.GzipFile(tf),tf)
                if f.endswith('.bz2') or comp == 'bz2':
                    h = handle_wrapper(bz2.BZ2File(tf),tf)
                else:
                    tfh = open(tf,'rbU')
                    h = handle_wrapper(tfh,tf)
            else:
                raise "Unsupported I/O spec: '%s'"%(spec,)
        elif f.endswith('.gz') or comp in ('gz','tgz'):
            if spec == 'r':
                h = gzip.open(f,'rb')
            elif spec == 'w':
                h = gzip.open(f,'wb')
            else:
                raise "Unsupported I/O spec: '%s'"%(spec,)
        elif f.endswith('.bz2') or comp == 'bz2':
            if spec == 'r':
                h = bz2.BZ2File(f,'rb')
            elif spec == 'w':
                h = bz2.BZ2File(f,'wb')
            else:
                raise "Unsupported I/O spec: '%s'"%(spec,)
        elif f == '-':
            if spec == 'r':
                h = sys.stdin
            elif spec == 'w':
                h = sys.stdout
            else:
                raise "Unsupported I/O spec: '%s'"%(spec,)
        else:
            if spec == 'r':
                h = open(f,spec+'bU')
            else:
                h = open(f,spec+'b')
    elif type(f) == file:
        h = f
    return h

class lfwrapper:
    def __init__(self,h):
        self.h = h
        self.re = re.compile(r'(\r|\r\n|\n)')
        self.pos = -1

    def __iter__(self):
        return self.next()

    def seek(p):
        assert (p == 0)
        self.h.seek(p)
        self.pos = -1

    def next(self):
        if self.pos < 0:
            self.buf = self.h.read(8196)
            self.lastchunk = ''
        while self.buf:
            sb = self.re.split(self.buf)
            if len(sb) > 1:
                if self.pos < 0:
                    # print "0>"+(self.lastchunk+sb[0])
                    self.pos = 0
                    yield (self.lastchunk+sb[0]+'\n')
                start = self.pos
                for i in range(start+2,len(sb)-1,+2):
                    # print str(i)+">"+sb[i]
                    self.pos = i
                    yield (sb[i]+'\n')
                self.lastchunk = sb[-1]
            else:
                self.lastchunk += sb[-1]
            self.buf = self.h.read(8196)
            self.pos = -1
        # print ">"+self.lastchunk
        yield self.lastchunk+'\n'

    def __getattr__(self,attr):
        return getattr(self.h,attr)

class handle_wrapper:
    def __init__(self,h,tf=None):
        self.h = h
        self.tf = tf

    def __del__(self):
        if self.tf:
            os.unlink(self.tf)

    def __getattr__(self,attr):
        return getattr(self.h,attr)

class charCounter:
    def __init__(self,h):
        self.h = h
        self.count = 0;

    def __getattr__(self,attr):
        return getattr(self.h,attr)

    def chars(self):
        return self.count

    def write(self,str):
        self.count += len(str);
        return self.h.write(str)
