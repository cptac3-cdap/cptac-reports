
import tempfile,os.path,sys,re,itertools,random,math,heapq
import gzip
import bz2
import warnings
from pickle import load, dump, HIGHEST_PROTOCOL

class OutOfCoreTable:
    def __init__(self,rows,force=False,headers=None):
        if headers:
            self.headers_ = headers
        else:
            self.headers_ = rows.headers()
        self.oocfiles = None
        self.therows = None
        self.maxrows = 50000
        self.force = force
        self.from_rows(rows)

    def __del__(self):
        if not self.oocfiles:
            return
        for f in self.oocfiles:
            try:
                os.unlink(f)
            except IOError:
                pass

    def __iter__(self):
        return self.rows()

    def __next__(self):
        return self.rows()

    def merge(self,files):
        for f in files:
            h=open(f,'rb')
            while True:
                try:
                    r = load(h)
                except EOFError:
                    break
                yield r
            h.close()

    def rows(self):
        if self.therows != None:
            for r in self.therows:
                yield dict(list(zip(self.headers_,r)))
        else:
            for r in self.merge(self.oocfiles):
                yield dict(list(zip(self.headers_,r)))

    def from_rows(self,rows):
        rows = iter(rows)
        self.oocfiles = []
        self.therows = []
        i = 0
        while True:
            try:
                r = next(rows)
            except StopIteration:
                break
            i += 1
            self.therows.append(list(map(r.get,self.headers_)))
            if i >= self.maxrows:
                tfd,tfn = tempfile.mkstemp()
                os.close(tfd)
                # print >>sys.stderr, "Writing out file:",tfn,"...",
                wh = open(tfn,'wb')
                for r in self.therows:
                    dump(r,wh,HIGHEST_PROTOCOL)
                wh.close()
                # print >>sys.stderr, 'done'
                self.oocfiles.append(tfn)
                self.therows = []
                i = 0
        if self.force or len(self.oocfiles) > 0:
            tfd,tfn = tempfile.mkstemp()
            os.close(tfd)
            # print >>sys.stderr, "Writing out last file:",tfn,"...",
            wh = open(tfn,'wb')
            for r in self.therows:
                dump(r,wh,HIGHEST_PROTOCOL)
            wh.close()
            # print >>sys.stderr, 'done'
            self.oocfiles.append(tfn)
            self.therows = None

class OutOfCoreSortedTable:
    def __init__(self,rows,key,force=False,keyisrow=False,headers=None,nodupkeys=False):
        if headers:
            self.headers_ = headers
        else:
            self.headers_ = rows.headers()
        self.key_ = key
        self.oocfiles = None
        self.therows = None
        self.maxrows = 50000
        self.force = force
        self.keyisrow = keyisrow
        self.nodupkeys = nodupkeys
        self.from_rows(rows)

    def __del__(self):
        if not self.oocfiles:
            return
        for f in self.oocfiles:
            try:
                os.unlink(f)
            except IOError:
                pass

    def __iter__(self):
        return self.rows()

    def __next__(self):
        return self.rows()

    def merge(self,files):
        rs = []
        hs = []
        lastk = None
        for i,f in enumerate(files):
            h=open(f,'rb')
            hs.append(h)
            try:
                k,r = load(h)
            except EOFError:
                continue
            heapq.heappush(rs, (k,r,i))
        while True:
            try:
                k,r,i = heapq.heappop(rs)
            except IndexError:
                break
            if not(self.nodupkeys and lastk and k == lastk):
                yield k,r
                lastk = k
            try:
                k,r = load(hs[i])
            except EOFError:
                continue
            heapq.heappush(rs, (k,r,i))
        for h in hs:
            h.close()

    def rows(self):
        if self.therows != None:
            lastk = None
            if not self.keyisrow:
                for k,r in self.therows:
                    if not(self.nodupkeys and lastk and k == lastk):
                        yield dict(list(zip(self.headers_,r)))
                        lastk = k
            else:
                for k,r in self.therows:
                    if not(self.nodupkeys and lastk and k == lastk):
                        yield dict(list(zip(self.headers_,k)))
                        lastk = k
            return

        if not self.keyisrow:
            for k,r in self.merge(self.oocfiles):
                yield dict(list(zip(self.headers_,r)))
        else:
            for k,r in self.merge(self.oocfiles):
                yield dict(list(zip(self.headers_,k)))

    def from_rows(self,rows):
        rows = iter(rows)
        self.oocfiles = []
        self.therows = []
        thelastrow = None
        while True:
            try:
                r = next(rows)
            except StopIteration:
                break
            k = self.key_(r)
            if not self.keyisrow:
                therow = (k,list(map(r.get,self.headers_)))
            else:
                therow = (k,None)
            if self.nodupkeys and thelastrow and therow[0] == thelastrow[0]:
                continue
            self.therows.append(therow)
            thelastrow = therow
            if len(self.therows) >= self.maxrows:
                self.therows.sort()
                tfd,tfn = tempfile.mkstemp()
                os.close(tfd)
                # print >>sys.stderr, "Writing out file:",tfn,"...",
                wh = open(tfn,'wb')
                for r in self.therows:
                    dump(r,wh,HIGHEST_PROTOCOL)
                wh.close()
                # print >>sys.stderr, 'done'
                self.oocfiles.append(tfn)
                self.therows = []
                i = 0
        self.therows.sort()
        if self.force or len(self.oocfiles) > 0:
            tfd,tfn = tempfile.mkstemp()
            os.close(tfd)
            # print >>sys.stderr, "Writing out last file:",tfn,"...",
            wh = open(tfn,'wb')
            for r in self.therows:
                dump(r,wh,HIGHEST_PROTOCOL)
            wh.close()
            # print >>sys.stderr, 'done'
            self.oocfiles.append(tfn)
            self.therows = None

        bs = 5
        while len(self.oocfiles) > bs:
            mergefiles = self.oocfiles[:bs]
            tfd,tfn = tempfile.mkstemp()
            os.close(tfd)
            # print >>sys.stderr, "Writing out merged file:",tfn,"...",
            wh = open(tfn,'wb')
            lastr = None
            for r in self.merge(mergefiles):
                if r != lastr:
                    dump(r,wh,HIGHEST_PROTOCOL)
                    lastr = r
            wh.close()
            # print >>sys.stderr, 'done'
            self.oocfiles.append(tfn)
            del self.oocfiles[:bs]
            for f in mergefiles:
                try:
                    os.unlink(f)
                except IOError:
                    pass

class OutOfCoreDistinct:
    def __init__(self,values):
        self.ooc = OutOfCoreSortedTable(OutOfCoreDistinct.values(values),headers=['value'],
                                        key=lambda x: (x.get('value'),), keyisrow = True, nodupkeys = True)

    @staticmethod
    def values(it):
        for v in it:
            yield dict(value=v)

    def __iter__(self):
        return self.rows()

    def rows(self):
        first = True
        for r in self.ooc.rows():
            v = r.get('value')
            if first:
                yield v
                first = False
                lastv = v
            elif v != lastv:
                yield v
                lastv = v
