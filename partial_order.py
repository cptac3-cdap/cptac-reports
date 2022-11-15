
from collections import defaultdict

def partial_order(proteins,peptides,edges,prpep2align):
    ordinal = dict()
    inedge = defaultdict(set)
    outedge = defaultdict(set)
    allpeps = set()
    for pr in proteins:
        ordered_peps = sorted(edges[pr]&peptides,
                              key=lambda p: (prpep2align.get((pr,p)),p))
        for i in range(0,len(ordered_peps)-1):
            inedge[ordered_peps[i+1]].add(ordered_peps[i])
            outedge[ordered_peps[i]].add(ordered_peps[i+1])
        allpeps.update(ordered_peps)
    s = set()
    for pep in allpeps:
        if len(inedge[pep]) == 0:
            s.add(pep)
    if len(s) == 0:
        s.add(min(allpeps))
    nextlabel=1
    labeled = set()
    while len(s) > 0:
        u = s.pop()
        ordinal[u] = nextlabel
        labeled.add(u)
        nextlabel += 1
        for v in outedge[u]:
            inedge[v].remove(u)
            if len(inedge[v]) == 0 and v not in labeled:
                s.add(v)
    return ordinal
