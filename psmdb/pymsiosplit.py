import re, os
srre1 = re.compile(r'^(.*?)(-(\d+))?(\.d(\d+))?\.(\w[\w\d_]+)\.(\d+)(\.(gz|bz2|zip|tgz))?$')
def splitResultFileExtn(f):
    m = srre1.search(f)
    if m == None:
        return None
    (p,b) = os.path.split(m.group(1))
    decoy = m.group(5)
    if decoy != None:
        decoy=int(decoy)
    comp = m.group(9)
    if comp == None:
        comp = ''
    chunk = m.group(3)
    if chunk != None:
        chunk = int(chunk)
    fullbase = b
    if m.group(2):
        fullbase += m.group(2)
    if m.group(4):
        fullbase += m.group(4)
    return dict(path=p,base=b,fullbase=fullbase,
                chunk=chunk,decoy=decoy,
                engine=m.group(6).lower(),instance=int(m.group(7)),
                compression=comp)

srre2 = re.compile(r'^(.*?)(\.raw.*?)?\.(mzid|mgf|raw|mzml)(\.(gz|bz2))?$',re.I)
def splitResultFileMZID(f):
    m = srre2.search(f)
    if m == None:
        return None
    (p,b) = os.path.split(m.group(1))
    comp = m.group(4)
    return dict(path=p,base=b,compression=comp)
