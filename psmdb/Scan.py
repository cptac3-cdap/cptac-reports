
import re

class Scan(tuple):
    def __new__(cls,*args):
        if len(args) == 1:
            try:
                return tuple.__new__(cls,(int(args[0]),))
            except (ValueError,TypeError) as e:
                pass
            try:
                return tuple.__new__(cls,list(map(int,args[0])))
            except (ValueError,TypeError):
                pass
            try:
                if isinstance(args[0],str):
                    return tuple.__new__(cls,list(map(int,args[0].split('.'))))
            except (ValueError,TypeError):
                pass
        return tuple.__new__(cls,list(map(int,args)))
    def __str__(self):
        return '.'.join(map(str,self))
    regex = r'\s*(\s*\d+\s*\.)*\s*\d+\s*'
    pairregex = '(?P<start_scan>' + regex + ')' + '(-(?P<end_scan>' + regex + '))?'
    @staticmethod
    def search(s,prefix="",suffix="",flags=0):
        m = re.search(prefix+Scan.pairregex+suffix,s,flags)
        if not m:
            return 0,None,None
        if m.group('end_scan') != None:
            return 2,Scan(m.group('start_scan')),Scan(m.group('end_scan'))
        return 1,Scan(m.group('start_scan')),Scan(m.group('start_scan'))

if __name__=='__main__':

    import sys

    print(Scan(1,2,3))
    print(Scan([1,2,3]))
    print(Scan((1,2,3)))
    print(Scan(' 1','2 ','3'))
    print(Scan((1,)))
    print(Scan(1))
    print(Scan('1'))
    print(Scan('1.2.3'))

    print(Scan.search('Scan Number: 290', prefix='^Scan +Number:? *',suffix="$",flags=re.IGNORECASE))
    print(Scan.search('Scan Number: 290 ', prefix='^Scan +Number:? *',suffix="$",flags=re.IGNORECASE))
    print(Scan.search('Scan Number: 290', prefix="Scan Number: "))
    print(Scan.search(' 1.2.3'))
    print(Scan.search(' 1.2.3 - 3. 2'))
    print(Scan.search('Locus:1.2.3-3.2',prefix="Locus:\s*"))
    print(Scan.search('1'))
    print(Scan.search('1-2'))
