
import array

class SearchResult:
    def __init__(self):
        self.psm = []
        self.md_ = {}
    def addPSM(self,psm):
        self.psm.append(psm)
    def sortPSMs(self,keyname=None):
        if not keyname:
            keyname = self.get('ranking_key','rank')
        if keyname.startswith('-'):
            self.psm.sort(key=lambda m: m.get(keyname[1:]),reverse=True)
        else:
            self.psm.sort(key=lambda m: m.get(keyname))
    def getPSM(self,i):
        return self.psm[i-1]
    def findPSM(self,key,value):
        l = [ m for m in self.psm if m.get(key) == value ]
        if len(l) == 0:
            return None
        elif len(l) == 1:
            return l[0]
        else:
            raise RuntimeError("Non-unique key value pair in findPSM")
    def setPSMs(self,psms):
        self.psm = psms
    def getPSMs(self):
        return self.psm
    def PSMiter(self):
        for m in self.psm:
            yield m
    def __len__(self):
        return len(self.psm)
    def has(self,key):
        return (key in self.md_)
    def get(self,key,df=None):
        return self.md_.get(key,df)
    def set(self,key,value):
        self.md_[key] = value
    def rename(self,old_key,new_key):
        self.md_[new_key] = self.md_[old_key]
        del self.md_[old_key]
    def copy(self,old_key,new_key):
        self.md_[new_key] = self.md_[old_key]
    def unset(self,key):
        if key in self.md_:
            del self.md_[key]
    def metadata(self):
        return self.md_
    def __str__(self):
        s = ">%s"%self.get('spectrum',"%s.%s.%s"%(self.get('base'),self.get('start_scan'),self.get('end_scan')))
        for (k,v) in sorted(iter(list(self.md_.items())),key=lambda k_v: k_v[0]):
            if k not in ('spectrum','ranking_key'):
                s += " %s:%s"%(k,v)
        s += '\n'
        for (i,h) in enumerate(self.psm):
            s += "%d. "%(i+1,)
            s += str(h)
            s += '\n'
        return s[:-1]

class PeptideSpectrumMatch:
    def __init__(self):
        self.md_ = {}
    def get(self,key,df=None):
        return self.md_.get(key,df)
    def has(self,key):
        return (key in self.md_)
    def set(self,key,value):
        self.md_[key] = value
    def append(self,key,value,delim=','):
        if key not in self.md_:
            self.md_[key] = str(value)
        else:
            if type(self.md_[key]) != type(''):
                self.md_[key] = str(self.md_[key])
            self.md_[key] += (delim + str(value))
    def sort(self,key,value,delim=','):
        if key not in self.md_:
            return
        if type(self.md_[key]) != type(''):
            return
        l = self.md_[key].split(delim)
        if len(l) < 2:
            return
        l.sort(key=value)
        self.md_[key] = ','.join(l)
    def unset(self,key):
        if key in self.md_:
            del self.md_[key]
    def rename(self,old_key,new_key):
        self.md_[new_key] = self.md_[old_key]
        del self.md_[old_key]
    def copy(self,old_key,new_key):
        self.md_[new_key] = self.md_[old_key]
    def metadata(self):
        return self.md_
    def __str__(self):
        s = "%s"%self.get('peptide')
        for (k,v) in sorted(iter(list(self.md_.items())),key=lambda k_v1: k_v1[0]):
            if k not in ('peptide'):
                s += " %s:%s"%(k,v)
        return s

