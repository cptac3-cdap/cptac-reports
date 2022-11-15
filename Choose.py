
import operator
from functools import reduce

def fact(n):
    return reduce(operator.mul,list(range(2,n+1)),1)

def combinations(n,k):
    return fact(n)/(fact(n-k)*fact(k))

class choose:
    """
    choose(l,k) -> iterator for all choices of k items from iterable l.
    choose(n,k) -> iterator for all choices of k items from range(0,n).

    The choose object yields tuples containing k list elements chosen,
    without replacement, from l or range(0,n).
    """

    def __init__(self,l,k):
        if type(l) == type(0):
            self.l = list(range(0,l))
        else:
            self.l = list(l)
        self.n = len(self.l)
        self.n1 = self.n - 1
        self.k = k
        self.k1 = self.k - 1
        self.nk = self.n - self.k

    def __iter__(self):
        return next(self)

    def __next__(self):
        if self.k == 0:
            yield ()
            return
        if len(self.l) < self.k:
            return
        self.x = list(range(0,self.k))
        yield self.listelts()
        while self.x[0] < self.nk:
            j = self.k1
            if self.x[j] < self.n1:
                self.x[j] += 1
                yield self.listelts()
            else:
                while (self.x[j] - self.x[j-1]) == 1:
                    j -= 1
                self.x[j-1] += 1
                while j < self.k:
                    self.x[j] = self.x[j-1]+1
                    j += 1
                yield self.listelts()

    def listelts(self):
        return tuple([ self.l[i] for i in self.x])

if __name__ == '__main__':

    import sys
    l = sys.argv[1]
    k = int(sys.argv[2])

    # help(choose)

    print(fact(int(l)), fact(int(k)), fact(int(l)-int(k)))
    print(combinations(int(l),int(k)))

    for t in choose(l,k):
        print(t)

    for t in choose(int(l),k):
        print(t)
