
def filter_criteria(psm):
    r={}
    if psm.has('estfdr'):
        r['estfdr'] = float(psm.get('estfdr'))
    if psm.has('qvalue'):
        r['estfdr'] = float(psm.get('qvalue'))
    if psm.has('pepfdr'):
        r['pepfdr'] = float(psm.get('pepfdr'))
    if psm.has('peptideprophet.probability'):
        r['probability'] = float(psm.get('peptideprophet.probability'))
        r['estfdr'] = 1-r['probability']
    if psm.has('probability'):
        r['probability'] = float(psm.get('probability'))
        r['estfdr'] = 1-r['probability']
    if psm.has('eval'):
        r['evalue'] = float(psm.get('eval'))
    if psm.has('expect'):
        r['evalue'] = float(psm.get('expect'))
    if psm.has('evalue'):
        r['evalue'] = float(psm.get('evalue'))
    if psm.has('hit_rank'):
        r['rank'] = int(psm.get('hit_rank'))
    if psm.has('notierank'):
        r['ntrank'] = int(psm.get('notierank'))
    if psm.has('rank'):
        r['rank'] = int(psm.get('rank'))
    if psm.has('trank'):
        r['trank'] = int(psm.get('trank'))
    if psm.has('peptide'):
        r['length'] = len(psm.get('peptide'))
    return r

class PSMFilter(object):
    def __init__(self,field,minimize,tight,loose=None,maxrank=1,maxntrank=1e+20,maxtrank=1e+20,minpeplen=7):
        self.tight = tight
        if loose == None:
            loose = tight
        else:
            loose = 1e+20
        self.minimize = minimize
        self.field = field
        self.maxrank = maxrank
        self.maxntrank = maxntrank
        self.maxtrank = maxtrank
        self.minpeplen = minpeplen
        if self.minimize:
            self._lt = float.__lt__
            self._gt = float.__gt__
            self.loose = max(tight,loose)
        else:
            self.loose = min(tight,loose)
            self._lt = float.__gt__
            self._gt = float.__lt__
    def retain(self,psm):
        return self.loose_filter(filter_criteria(psm))
    def loose_filter(self,r):
        return self.good(r,self.loose)
    def tight_filter(self,r):
        return self.good(r,self.tight)
    def good(self,r,threshold):
        return self.lteq(self.value(r),threshold) and r['rank'] <= self.maxrank and r.get('ntrank',99) <= self.maxntrank and r.get('trank',99) <= self.maxtrank and r['length'] >= self.minpeplen
    def prefer(self,r,threshold):
        if threshold == None:
            return True
        return self.lt(self.value(r),threshold)
    def value(self,r):
        if self.field == None:
            return None
        if isinstance(r,float):
            return r
        return float(r[self.field])
    def psmqvalue(self,psm):
        return self.qvalue(filter_criteria(psm))
    def psmcriteria(self,psm):
        return filter_criteria(psm)
    def qvalue(self,r):
        if self.field == None:
            return None
        if self.minimize:
            return self.value(r)
        # For probability
        return (1-self.value(r))
    def lt(self,value1,value2):
        return self._lt(value1,value2)
    def lteq(self,value1,value2):
        if self.field == None:
            return True
        return (not self._gt(value1,value2))
    def cmp(self,value1,value2):
        if self._lt(value1,value2):
            return -1
        if self._gt(value1,value2):
            return 1
        return 0
    def tight_filter_sql(self,column):
        if self.minimize:
            return (column <= self.tight)
        return (column >= self.tight)

    def loose_filter_sql(self,column):
        if self.minimize:
            return (column <= self.loose)
        return (column >= self.loose)

class PeptideProphet(PSMFilter):
    def __init__(self,tight,**kw):
        PSMFilter.__init__(self,'probability',False,tight,None,**kw)

class SpectralFDR(PSMFilter):
    def __init__(self,tight,**kw):
        PSMFilter.__init__(self,'estfdr',True,tight,None,**kw)

class PeptideFDR(PSMFilter):
    def __init__(self,tight,**kw):
        PSMFilter.__init__(self,'pepfdr',True,tight,None,**kw)

class EValue(PSMFilter):
    def __init__(self,tight,**kw):
        PSMFilter.__init__(self,'evalue',True,tight,None,**kw)

class MZIdPassThreshold(PSMFilter):
    def __init__(self,**kw):
        PSMFilter.__init__(self,'passThreshold',True,"true",**kw)

class NoFilter(PSMFilter):
    def __init__(self,tight=1e+20,**kw):
        PSMFilter.__init__(self,None,None,tight=tight,**kw)
