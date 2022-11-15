from collections import defaultdict
from operator import itemgetter
from copy import copy, deepcopy
from Choose import choose
import time
import heapq
import math
import sys
from functools import reduce, cmp_to_key

class BranchAndBound:
    def __init__(self,**kw):
        self.heap = []
        self.quiet = kw.get('quiet',False)
        self.log_interval = kw.get('interval',5)
        self.key = self.greedy
    def greedy(self,state):
        return self.cost(state)
    def most(self,state):
        return -len(state.chosen),-state.depth
    def bfs(self,state):
        return state.depth,self.cost(state)
    def dfs(self,state):
        return -state.depth,self.cost(state)
    def bylowerbound(self,state):
        return state.lb
    def rekey(self, key=None, ub=None):
        if key:
            self.key = key
        h = []
        for k,s in self.heap:
            if ub != None and s.lb >= ub:
                continue
            heapq.heappush(h,self.row(s))
        self.heap = h
    def row(self, state):
        return self.key(state),state
    def push(self, state, depth=None, lb=None):
        if depth != None:
            state.depth=depth
        if lb != None:
            state.lb = lb
        heapq.heappush(self.heap,self.row(state))
    def pop(self):
        try:
            return heapq.heappop(self.heap)
        except IndexError:
            pass
        return None,None
    def solve(self, exub=1e+20, exlb=-1e+20,
              forcein=None, strategy='greedy',
              tinitial=1e+20, timprove=1e+20,ttotal=1e+20):
        assert strategy in ('greedy','most','bfs','dfs','bylowerbound')
        self.key = getattr(self,strategy)
        start_time = time.time()
        end_time = start_time + min(tinitial,ttotal)
        state = self.empty_solution()
        if forcein:
            for pr in forcein:
                state.select_protein(pr)
        self.push(state,depth=0,lb=self.lower_bound(state))
        best = (exub,None)
        iter = 0; last = time.time()
        while time.time() < end_time:
            key,state = self.pop()
            if state == None:
                break
            if not self.quiet and time.time() > last + self.log_interval:
                print("Iteration %7d nodes %7d best %7s current: %s"%(iter,len(self.heap),best[0],self.logstats(state)))
                sys.stdout.flush()
                last = time.time()
            iter += 1
            # print state.lb, best[0], state.lb >= best[0]
            if state.lb >= best[0]:
                continue
            if self.feasible(state):
                cost = self.cost(state)
                if cost < best[0]:
                    if best[1] == None or True:
                        end_time = min(time.time()+timprove, start_time+ttotal);
                    best = (cost,state)
                    if cost <= exlb:
                        break
                    self.newsolution(state)
                    self.rekey(ub=best[0])
                    if not self.quiet and len(self.heap) > 0:
                        print("Iteration %7d nodes %7d best %7s"%(iter,len(self.heap),best[0]))
                        sys.stdout.flush()
                continue
            if not self.explore_children(state,best[0]):
                continue
            for ch in self.children(state):
                lb = self.lower_bound(ch)
                if lb >= best[0]:
                    continue
                self.push(ch,depth=state.depth+1,lb=lb)
            if self.dorekey():
                print("Do rekey!", file=sys.stderr)
                sys.stderr.flush()
                self.rekey(ub=best[0])
        if len(self.heap) > 0:
            print("WARNING: Time limit reached. Solution not provably optimal!")
            sys.stdout.flush()
        if not self.quiet:
            print("Iteration %7d nodes %7d best %7s DONE"%(iter,len(self.heap),best[0]))
            sys.stdout.flush()
        return (best[0],best[1],len(self.heap) == 0)

class BranchAndBoundState:
    def __init__(self):
        pass

    def uitemset(self,pgsi,pri):
        return reduce(set.union,[self.uitems.get((pgi,pri),set()) for pgi in pgsi],set())

    def init(self,peps,prots,edges,weights,uitems=None,pr2key=None):

        # Preprocess to ensure that the protein set is minimal
        lprots = list(prots)
        pepedges = {}
        for pr in lprots:
            pepedges[pr] = (edges[pr]&peps)
        contained = set()
        for i in range(0,len(lprots)):
            for j in range(i+1,len(lprots)):
                if j in contained:
                    continue
                if pepedges[lprots[j]] <= pepedges[lprots[i]]:
                    contained.add(lprots[j])
        prots = copy(prots)
        for pr in contained:
            prots.remove(pr)
            del pepedges[pr]
        edges = pepedges

        self.weights = weights
        self.edges = edges
        self.uitems = uitems
        self.pr2key = pr2key
        self.invedges = defaultdict(set)
        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)
        self.chosen = set()
        self.intersection = {}
        for pr in prots:
            pepsi = (self.edges[pr]&peps)
            if uitems:
                uit = self.uitemset(pepsi,pr)
                self.intersection[pr] = (len(uit),sum(map(self.weights.get,pepsi)))
            else:
                self.intersection[pr] = (len(pepsi),sum(map(self.weights.get,pepsi)))
        self.coveredby = dict((pep,len(self.invedges[pep]&prots)) for pep in peps)
        notcovered = [t for t in list(self.coveredby.items()) if t[1]==0]
        self.uncovered = (len(notcovered),sum([self.weights.get(t[0]) for t in notcovered]))
        self.uniqpeps = defaultdict(set)

    def dump(self):
        print("I:",' '.join("%d:%d,%d"%(k,v[0],v[1]) for k,v in sorted(self.intersection.items())), file=sys.stderr)
        print("C:",' '.join("%d:%d"%(k,v) for k,v in sorted(self.coveredby.items())), file=sys.stderr)
        print(self.uncovered, file=sys.stderr)
        print(self.chosen, file=sys.stderr)

    def proteins_byintersection(self):
        return sorted(list(self.intersection.items()),key=lambda t: (-t[1][0],self.pr2key(t[0]),t[0]))

    def proteins_bywintersection(self):
        return sorted(list(self.intersection.items()),key=lambda t: (-t[1][1],self.pr2key(t[0]),t[0]))

    def max_protein_byintersection(self):
        try:
            return min(list(self.intersection.items()),key=lambda t: (-t[1][0],self.pr2key(t[0]),t[0]))
        except ValueError:
            pass
        return None,(0,0.0)

    def clone(self):
        child = BranchAndBoundState()
        child.edges = self.edges
        child.weights = self.weights
        child.uitems = self.uitems
        child.pr2key = self.pr2key
        child.invedges = self.invedges
        child.intersection = copy(self.intersection)
        child.coveredby = copy(self.coveredby)
        child.uncovered = self.uncovered
        child.chosen = copy(self.chosen)
        child.uniqpeps = copy(self.uniqpeps)
        return child

    def select_protein(child,pr):
        # print >>sys.stderr, "Selecting protein",pr
        del child.intersection[pr]
        child.chosen.add(pr)
        # These are the unqiue peptides for this new selected protein.
        if child.uitems:
            child.uniqpeps[pr] = child.uitemset(child.edges[pr]&set(child.coveredby),pr)
        else:
            child.uniqpeps[pr] = (child.edges[pr]&set(child.coveredby))
        child.minpepper = len(child.uniqpeps[pr])
        # Now update the other chosen proteins' unique lists
        for pr1 in child.chosen:
            if pr1 == pr:
                continue
            if child.uitems:
                child.uniqpeps[pr1] -= child.uitemset(child.edges[pr],pr1)
            else:
                child.uniqpeps[pr1] -= child.edges[pr]
            child.minpepper = min(len(child.uniqpeps[pr1]),child.minpepper)
        prottocheck = set()
        for pep in child.edges[pr]:
            if pep in child.coveredby:
                for pr1 in child.invedges[pep]:
                    if pr1 in child.intersection:
                        prottocheck.add(pr1)
                del child.coveredby[pep]
        prottozero = set()
        remainingpeps = set(child.coveredby)
        for pr1 in prottocheck:
            pepsi = (child.edges[pr1]&remainingpeps)
            if child.uitems:
                uit = child.uitemset(pepsi,pr1)
                child.intersection[pr1] = (len(uit),sum(map(child.weights.get,pepsi)))
            else:
                child.intersection[pr1] = (len(pepsi),sum(map(child.weights.get,pepsi)))
            if child.intersection[pr1][0] == 0:
                prottozero.add(pr1)
        return prottozero

    def reject_protein(child,pr):
        # print >>sys.stderr, "Rejecting protein",pr
        del child.intersection[pr]
        for pep in child.edges[pr]:
            if pep in child.coveredby:
                child.coveredby[pep] -= 1
                if child.coveredby[pep] == 0:
                    t = child.uncovered
                    child.uncovered = (t[0]+1,t[1]+child.weights.get(pep))
        return child

    def npeptide(self):
        return len(self.coveredby)

    def wpeptide(self):
        return sum(map(self.weights.get,self.coveredby))

    def nchosen(self):
        return len(self.chosen)

    def nprotein(self):
        return len(self.intersection)

class Parsimony(BranchAndBound):
    def __init__(self,proteins,peptides,edges,*args,**kw):
        self.left = 0
        self.minpepper = 1
        if 'left' in kw:
            self.left = kw['left']
            del kw['left']
        if 'minpepper' in kw:
            self.minpepper = kw['minpepper']
            del kw['minpepper']
        if kw.get('weights') != None:
            self.weights = kw['weights']
            del kw['weights']
        else:
            self.weights = dict((p,1) for p in peptides)
        if kw.get('uitems') != None:
            self.uitems = kw['uitems']
        else:
            self.uitems = None
        if kw.get('prsortkey') != None:
            self.pr2key = kw['prsortkey']
            del kw['prsortkey']
        else:
            self.pr2key = defaultdict(float)
        BranchAndBound.__init__(self,*args,**kw)
        self.prot = proteins
        self.pep = peptides
        self.edges = edges
        self.seensolution = False
        self.key = self.dfs

    def solve(self,*args,**kw):
        best = BranchAndBound.solve(self,*args,**kw)
        if best[1]:
            return best[0],best[1].chosen,best[2]
        return 1e+20,None,best[2]

    def empty_solution(self):
        state = BranchAndBoundState()
        state.init(self.pep,self.prot,self.edges,self.weights,self.uitems,self.pr2key)
        return state

    def greedy(self,state):
        return state.wpeptide()

    def bfs(self,state):
        return state.depth,state.wpeptide()

    def dfs(self,state):
        return -state.depth,state.wpeptide()

    def bylowerbound(self,state):
        return state.lb,self.dfs(state)

    def cost(self,state):
        return state.nchosen()

    def logstats(self,state):
        # state.dump()
        return "lb %7s proteins %7s peptides %7s chosen %7s"%(state.lb,state.nprotein(),
                                                              state.npeptide(),state.nchosen())

    def feasible(self,state):
        return state.wpeptide() <= self.left

    def lower_bound(self,state):
        wpep = state.wpeptide()
        if wpep <= self.left:
            return state.nchosen()
        bestcase = 0
        nbestcase = 0
        for pr,(intr,wintr) in state.proteins_bywintersection():
            if intr < self.minpepper:
                continue
            nbestcase += 1
            bestcase  += wintr
            if bestcase >= (wpep-self.left):
                break
        if bestcase < (wpep-self.left):
            return 1e+20
        lb1 = state.nchosen() + nbestcase
        return lb1

    def explore_children(self,state,ub):
        if state.uncovered[1] > self.left:
            return False
        return True

    def children(self,state):
        pr,(intr,wintr) = state.max_protein_byintersection()
        if intr >= self.minpepper:
            ch = state.clone()
            tozero = ch.select_protein(pr)
            if ch.minpepper >= self.minpepper:
                for pr1 in tozero:
                    ch.reject_protein(pr1)
                yield ch
            ch = state.clone()
            ch.reject_protein(pr)
            for pr1 in tozero:
                ch.reject_protein(pr1)
            yield ch

    def newsolution(self,state):
        if not self.seensolution:
            self.seensolution = True
            # self.key = self.bylowerbound

    def dorekey(self):
        return False
        if not self.seensolution:
            return False
        if len(self.heap) >= 100 and self.key != self.dfs:
            self.key = self.dfs
            return True
        if len(self.heap) < 10 and self.key != self.bylowerbound:
            self.key = self.bylowerbound
            return True
        return False

class MinUncovered(Parsimony):
    def __init__(self,proteins,peptides,edges,nprot=1e+20,*args,**kw):
        self.nprot = nprot
        Parsimony.__init__(self,proteins,peptides,edges,*args,**kw)

    def cost(self,state):
        return (state.wpeptide(),state.nchosen())

    def feasible(self,state):
        pr,(intr,wintr) = state.max_protein_byintersection()
        if intr < self.minpepper:
            return True
        if state.nchosen() >= self.nprot:
            return True
        return False

    def lower_bound(self,state):
        if self.feasible(state):
            return self.cost(state)
        nchosen = state.nchosen()
        bestcase = 0
        nbestcase = 0
        for pr,(intr,wintr) in state.proteins_bywintersection():
            if intr < self.minpepper:
                continue
            nbestcase += 1
            bestcase  += wintr
            if nbestcase >= (self.nprot-nchosen):
                break
        return (max(state.wpeptide() - bestcase,state.uncovered[1]),state.nchosen()+nbestcase)

    def explore_children(self,state,ub):
        return True

class Components(object):
    def __init__(self,proteins,peptides,edges,maxpepdeg=1e+20,progress=None):
        self.prot = proteins
        self.pep = peptides
        self.edges = edges
        self.invedges = defaultdict(set)
        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)
        badpep = set()
        for pep in list(self.invedges.keys()):
            if len(self.invedges[pep]) > maxpepdeg:
                badpep.add(pep)
        self.pr2comp = {}
        self.component = []
        if progress:
            progress.stage("Finding components...",len(self.prot))
        for pr in self.prot:
            if pr in self.pr2comp:
                continue
            thecomp = set()
            queue = set([pr])
            while len(queue) > 0:
                pr0 = queue.pop()
                thecomp.add(pr0)
                peps = ((self.edges[pr0] & self.pep) - badpep)
                pr1 = reduce(set.union,list(map(self.invedges.get,peps)),set())
                pr1 &= self.prot
                pr1 -= thecomp
                queue.update(pr1)
            for pr1 in thecomp:
                self.pr2comp[pr1] = len(self.component)-1
            self.component.append((thecomp,
                                   self.pep & reduce(set.union,
                                              list(map(self.edges.get,thecomp)))))
            if progress:
                progress.update()
        if progress:
            progress.done()
            progress.message("Components: %d"%len(self.component))

    def __iter__(self):
        return next(self)

    def __next__(self):
        for pr,pep in self.component:
            yield pr,pep

def mycmp(a,b):
    if a < b:
       return -1
    elif b < a:
       return 1
    return 0

class Dominator:
    def __init__(self,peptides,edges,proteins):
        self.peptides = set(peptides)
        self.proteins = set(proteins)
        self.edges = copy(edges)
        self.eliminated = set()
        self.equivalentto = defaultdict(set)
        self.containedby = defaultdict(set)
        self.invedges = defaultdict(set)
        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)

    def dominate(self,prsortkey=None,tiebreakers=None,eliminate=0,promiscuouslimit=1e+20,relationships=True,progress=None):
        pr2tbpep = defaultdict(set)
        tbpep2fdr = defaultdict(lambda :1e+20)
        tbcmp = mycmp
        if tiebreakers:
            pr2tbpep = tiebreakers[0]
            tbpep2fdr = tiebreakers[1]
            tbcmp = tiebreakers[2]
        if progress:
            progress.message("Proteins: %d peptides: %d"%(len(self.proteins),len(self.peptides)))
        if prsortkey != None:
            lprots = sorted(self.proteins,key=prsortkey)
        else:
            lprots = list(self.proteins)
        # print self.proteins
        peps = set()
        for pep in self.peptides:
            if len(self.invedges[pep]) <= promiscuouslimit:
                peps.add(pep)
        pepedges = {}
        for pr in lprots:
            pepedges[pr] = (self.edges[pr]&peps)
        invedges = defaultdict(set)
        for i,pr in enumerate(lprots):
            for v in pepedges[pr]:
                invedges[v].add(i)
        contained = set()
        if progress:
            progress.stage("Checking for dominated proteins",len(lprots))
        for i in range(0,len(lprots)):
            # print i,pepedges[lprots[i]]
            if progress:
                progress.update()
            if len(pepedges[lprots[i]]) <= eliminate:
                contained.add(lprots[i])
                if relationships:
                    self.eliminated.add(lprots[i])
                continue
            # print contained
            checked = set()
            for pep in pepedges[lprots[i]]:
                # print invedges[pep], checked, contained
                for j in invedges[pep]:
                    # print j, j in checked, j in contained
                    if j == i:
                        continue
                    if j in checked:
                        continue
                    checked.add(j)
                    # if lprots[j] in contained:
                    #     continue
                    # print " ",j,pepedges[lprots[j]]
                    if pepedges[lprots[i]] >= pepedges[lprots[j]]:
                        if pepedges[lprots[i]] > pepedges[lprots[j]]:
                            if relationships:
                                self.containedby[lprots[i]].add(lprots[j])
                            contained.add(lprots[j])
                        elif prsortkey != None and i < j:
                            if relationships:
                                self.equivalentto[lprots[i]].add(lprots[j])
                                self.containedby[lprots[i]].add(lprots[j])
                            contained.add(lprots[j])
                        else:
                            tbi = sorted(map(tbpep2fdr.get,pr2tbpep[lprots[i]]-pr2tbpep[lprots[j]]),key=cmp_to_key(tbcmp))
                            tbj = sorted(map(tbpep2fdr.get,pr2tbpep[lprots[j]]-pr2tbpep[lprots[i]]),key=cmp_to_key(tbcmp))
                            tbi.extend([1e+20]*max(len(tbj)-len(tbi),0))
                            tbj.extend([1e+20]*max(len(tbi)-len(tbj),0))
                            # if len(tbi) > 0:
                            #     print >>sys.stderr, (tbi,i),(tbj,j)
                            if (tbi,i) < (tbj,j):
                                # print >>sys.stderr, (tbi,i),(tbj,j)
                                if relationships:
                                    self.equivalentto[lprots[i]].add(lprots[j])
                                    self.containedby[lprots[i]].add(lprots[j])
                                contained.add(lprots[j])

        if progress:
            progress.done()
        # print contained
        self.proteins = set(lprots) - contained
        self.edges = pepedges
        try:
            self.peptides = reduce(set.union,list(map(self.edges.get,self.proteins)))
        except TypeError:
            self.peptides = set()
        self.invedges = defaultdict(set)

        # Fix equivalence and contained relationships...
        if relationships:
            for p in self.proteins:
                self.equivalentto[p] = self.collapse(p,self.equivalentto)
                self.containedby[p] = (self.collapse(p,self.containedby) - self.equivalentto[p])
            for k in list(self.equivalentto):
                if k not in self.proteins:
                    del self.equivalentto[k]
            for k in list(self.containedby):
                if k not in self.proteins:
                    del self.containedby[k]
            # print self.proteins
            # print self.equivalentto
            # print self.containedby

        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)
        if progress:
            progress.message("Proteins: %d peptides: %d"%(len(self.proteins),len(self.peptides)))

    @staticmethod
    def collapse(p0,edges):
        collapsed = set()
        collapsed.add(p0)
        # print collapsed
        tocheck = set()
        # print p0,edges[p0]
        if p0 in edges:
            tocheck.update(edges[p0])
        while len(tocheck) > 0:
            # print collapsed,tocheck
            p1 = tocheck.pop()
            collapsed.add(p1)
            if p1 in edges:
                tocheck.update(edges[p1]-collapsed)
        return collapsed

    def force(self,forced,progress=None):
        nunique = defaultdict(int)
        if progress:
            progress.stage("Counting unique peptides per protein",len(self.peptides))
        for pep in self.peptides:
            if progress:
                progress.update()
            hitset = self.invedges[pep]
            if len(hitset) == 1:
                nunique[next(iter(hitset))] += 1
        if progress:
            progress.done()
        self.forced = set()
        if progress:
            progress.stage("Finding forced proteins",len(self.proteins))
        for pr in self.proteins:
            if progress:
                progress.update()
            if nunique[pr] >= forced:
                self.forced.add(pr)
        if progress:
            progress.done()
            progress.message("Forced proteins: %d"%len(self.forced))
        self.proteins -= self.forced
        if progress:
            progress.stage("Determining uncovered peptides",len(self.forced))
        self.covered = set()
        for pr in self.forced:
            if progress:
                progress.update()
            self.covered.update(self.edges[pr])
        if progress:
            progress.done()
            progress.message("Covered peptides: %d"%len(self.covered))
        self.peptides -= self.covered
        self.dominate(relationships=False,progress=progress)

class FixedPoint:
    def __init__(self,peptides,edges,proteins,pepqval):
        self.peptides = set(peptides)
        self.proteins = set(proteins)
        # edges: pr -> peps
        self.edges = {}
        for pr in self.proteins:
            self.edges[pr] = (self.peptides & edges[pr])
        self.pepprob = dict([(t[0],1.0-t[1]) for t in iter(pepqval.items())])
        # print self.pepprob
        self.invedges = defaultdict(set)
        for k,s in self.edges.items():
            for v in s:
                self.invedges[v].add(k)
        self.edgeweight = defaultdict(float)
        for pep in self.invedges:
            n = len(self.invedges[pep])
            for pr in self.invedges[pep]:
                self.edgeweight[(pr,pep)] = 1.0/float(n)
        self.proteinprob = self.protprob()
        # print self.proteinprob

    def protprob(self):
        protprob = defaultdict(float)
        for pr in self.edges:
            falseprob = 1
            for pep in self.edges[pr]:
                # print >>sys.stderr, pr,pep,self.edgeweight[(pr,pep)],self.pepprob[pep],1.0-self.edgeweight[(pr,pep)]*self.pepprob[pep]
                falseprob *= (1.0-self.edgeweight[(pr,pep)]*self.pepprob[pep])
            protprob[pr] = 1-falseprob
        return protprob

    def adjedges(self):
        edgeweight = defaultdict(float)
        for pep in self.invedges:
            total = sum(self.proteinprob[pr] for pr in self.invedges[pep])
            for pr in self.invedges[pep]:
                try:
                    edgeweight[(pr,pep)] = self.proteinprob[pr]/total
                except ZeroDivisionError:
                    edgeweight[(pr,pep)] = 0.0
        return edgeweight

    def iterate(self,steps):
        for i in range(steps):
            self.edgeweight = self.adjedges()
            protprob = self.protprob()
            delta = 0.0
            for pr in self.proteinprob:
                delta = max(abs(self.proteinprob[pr] - protprob[pr]),delta)
            self.proteinprob = protprob
            iterations = i
            # print i,delta
            if delta < 1e-5:
                break
        if iterations == (steps-1):
            print("Warning: Fixed-point solution has not converged...", file=sys.stderr)
        return

class GraphDumper:
    def __init__(self,name,proteins,peptides,edges,prset,pepset,weights=None,protein_name=None,peptide_name=None):
        w = weights
        if not w:
            w = defaultdict(lambda x: 1)
        h = open('%s.nodes.tsv'%name,'w')
        print("NodeID\tType\tWeight"+"\tName" if (peptide_name or protein_name) else "", file=h)
        for pri in prset:
            print("PR%d\tProtein\t%s"%(pri,0)+"\t%s"%protein_name(pri) if protein_name else "", file=h)
        for pepi in pepset:
            print("PE%d\tPeptide\t%s"%(pepi,w[pepi])+"\t%s"%peptide_name(pepi) if peptide_name else "", file=h)
        h.close()
        h = open('%s.edges.tsv'%name,'w')
        print("Protein\tPeptide", file=h)
        for pri in prset:
            for pepi in edges[pri]:
                if pepi in pepset:
                    print("PR%s\tPE%s"%(pri,pepi), file=h)
        h.close()
