#!/bin/env python

from __future__ import print_function

__all__ = [ 'permutations', 'select', 'choose', 'product', 'list_accumulator',
            'tuple_accumulator', 'concat_accumulator', 'set_accumulator', 'itermatchings']

import copy, sys
from collections import defaultdict
from operator import itemgetter
import itertools
from functools import reduce

def permutations(x):
    """

    Iterator for all permutations of an iterable

    >>> from combinatorics import permultations
    >>> for p in permutations([1,4,6,10]):                         
    ...     print p                                                
    ...                                                            
    [1, 4, 6, 10]                                                  
    [1, 4, 10, 6]                                                  
    [1, 6, 4, 10]                                                  
    [1, 6, 10, 4]                                                  
    [1, 10, 4, 6]                                                  
    [1, 10, 6, 4]                                                  
    [4, 1, 6, 10]                                                  
    [4, 1, 10, 6]                                                  
    [4, 6, 1, 10]                                                  
    [4, 6, 10, 1]                                                  
    [4, 10, 1, 6]                                                  
    [4, 10, 6, 1]                                                  
    [6, 1, 4, 10]                                                  
    [6, 1, 10, 4]                                                  
    [6, 4, 1, 10]                                                  
    [6, 4, 10, 1]                                                  
    [6, 10, 1, 4]                                                  
    [6, 10, 4, 1]                                                  
    [10, 1, 4, 6]                                                  
    [10, 1, 6, 4]                                                  
    [10, 4, 1, 6]                                                  
    [10, 4, 6, 1]                                                  
    [10, 6, 1, 4]                                                  
    [10, 6, 4, 1]                                                  
                                                                   
    """
    if len(x) <= 1:
        yield x
    else:
        for i in range(len(x)):
            for p in permutations(list(x[:i]) + list(x[i + 1:])):
                yield [x[i]] + list(p)

def select(x,k):
    """

    Iterator for all selections (order is significant) of k items from
    an iterable

    >>> from combinatorics import select
    >>> for p in select([1,4,6,10],2):                         
    ...     print p                                                
    ...                                                            
    [1, 4]  
    [1, 6]  
    [1, 10] 
    [4, 1]  
    [4, 6]  
    [4, 10] 
    [6, 1]  
    [6, 4]  
    [6, 10] 
    [10, 1] 
    [10, 4] 
    [10, 6] 
                                                                   
    """
    if k > len(x):
        return
    if k == 0:
        yield []
    for i in range(len(x)):
        for s in select(list(x[:i]) + list(x[i + 1:]),k-1):
            yield [x[i]] + list(s)
    
def choose(x,k):
    """

    Iterator for all choices (order is not significant) of k items from
    an iterable

    >>> from combinatorics import choose
    >>> for p in choose([1,4,6,10],2):                         
    ...     print p                                                
    ...                                                            
    [1, 4]  
    [1, 6]  
    [1, 10] 
    [4, 1]  
    [4, 6]  
    [4, 10] 
    [6, 1]  
    [6, 4]  
    [6, 10] 
    [10, 1] 
    [10, 4] 
    [10, 6] 
                                                                   
    """
    if k > len(x):
        return
    if k == 0:
        yield []
    for i in range(len(x)):
        for s in choose(list(x[i+1:]),k-1):
            yield [x[i]] + list(s)
    
def subset(x):
    for i in range(len(x)+1):
        for s in choose(x,i):
            return s

class set_accumulator:
    @staticmethod
    def new(x):
        return set(x)
    @staticmethod
    def add(l,x):
        l.add(x)
        return l

class list_accumulator:
    @staticmethod
    def new(x):
        return [x]
    @staticmethod
    def add(l,x):
        l.append(x)
        return l

class tuple_accumulator:
    @staticmethod
    def new(x):
        return tuple([x])
    @staticmethod
    def add(l,x):
        l = list(l)
        l.append(x)
        return tuple(l)

class concat_accumulator:
    @staticmethod
    def new(x):
        return str(x)
    @staticmethod
    def add(l,x):
        return (l + str(x))

class matching_accumulator:
    @staticmethod
    def new(x):
        return list(x[0]),list(x[1])
    @staticmethod
    def add(l,x):
        l = [ list(l[0]), list(l[1]) ]
        l[0].extend(list(x[0]))
        l[1].extend(list(x[1]))
        return l

def product(*args,**kw):
    """
    Cartesian products of at least two iterables as a list

    >>> from combinatorics import product
    >>> for p in product('abc',[1,2,3],['ab','cd']):
    ....    print p
    ....
    ('a', 1, 'ab')
    ('a', 1, 'cd')
    ('a', 2, 'ab')
    ('a', 2, 'cd')
    ('a', 3, 'ab')
    ('a', 3, 'cd')
    ('b', 1, 'ab')
    ('b', 1, 'cd')
    ('b', 2, 'ab')
    ('b', 2, 'cd')
    ('b', 3, 'ab')
    ('b', 3, 'cd')
    ('c', 1, 'ab')
    ('c', 1, 'cd')
    ('c', 2, 'ab')
    ('c', 2, 'cd')
    ('c', 3, 'ab')
    ('c', 3, 'cd')
    >>> from combinatorics import product, concat_accumulator
    >>> for p in product('abc',[1,2,3],['ab','cd'],
                         accumulator=concat_accumulator):
    ....    print p
    ....
    a1ab
    a1cd
    a2ab
    a2cd
    a3ab
    a3cd
    b1ab
    b1cd
    b2ab
    b2cd
    b3ab
    b3cd
    c1ab
    c1cd
    c2ab
    c2cd
    c3ab
    c3cd
    """
    if 'accumulator' not in kw:
        accumulator = tuple_accumulator
    else:
        accumulator = kw['accumulator']
    if len(args) == 0:
        return []
    l0 = map(accumulator.new,args[0])
    for l in args[1:]:
        l1 = []
        for i0 in l0:
            for i in l:
                l1.append(accumulator.add(copy.copy(i0),i))
        l0 = l1
    return l0

def itermatchings(items1,items2,matchtest):
    eq = dict()
    list1 = list(items1)
    n1 = len(list1)
    list2 = list(items2)
    n2 = len(list2)

    if n1 != n2:
        raise StopIteration

    for inds in itertools.permutations(range(n2)):
        badmatch = False
        for ai,bi in zip(range(n1),inds):
            if (ai,bi) in eq:
                if not eq[(ai,bi)]:
                    badmatch = True
                    break
            else:
                a = list1[ai]; b = list2[bi]
                if not matchtest(a,b):
                    eq[(ai,bi)] = False
                    badmatch = True
                    break
                else:
                    eq[(ai,bi)] = True
        if not badmatch:
            yield list1,list(map(list2.__getitem__,inds))

def iterplacements(items1,items2):
    list1 = list(items1)
    list2 = list(items2)
    n1 = len(list1)
    n2 = len(list2)
    if n2 < n1:
        raise StopIteration
    if n2 == 0:
        yield [],[]
        return
    args = []
    for j in range(n2):
        args.append(range(n1))
    for pos in product(*args,**dict(accumulator=list_accumulator)):
        if len(set(pos)) != n1:
            continue
        result={}
        for i in range(n1):
            result[i] = []
        for j,i in enumerate(pos):
            result[i].append(j)
        retval = []
        for i,js in result.items():
            retval.append((list1[i],map(lambda j: list2[j],js)))
        yield list(map(itemgetter(0),retval)),list(map(itemgetter(1),retval))

def iterpairs(items1,items2):
    list1 = list(items1)
    n1 = len(list1)
    list2 = list(items2)
    n2 = len(list2)

    if n1 != n2:
        raise StopIteration
    for inds in itertools.permutations(range(n2)):
        yield list(zip(list1,map(list2.__getitem__,inds)))

def iterecmatchings(items1,items2,matchtest):
    list1 = list(items1)
    n1 = len(list1)
    list2 = list(items2)
    n2 = len(list2)

    # print list1
    # print list2
    # print n1,n2

    if n1 != n2:
        raise StopIteration

    ec1 = defaultdict(list)
    ec2 = defaultdict(list)
    
    placed = set()
    for i in range(n1):
        if i in placed:
            continue
        placed.add(i)
        ec1[i].append(i)
        for j in range(i+1,n1):
            if j in placed:
                continue
            # print i,j,
            # sys.stdout.flush()
            if not matchtest(list1[i],list1[j]):
                # print "no match"
                # sys.stdout.flush()
                continue
            # else:
            #   print "match"
            #   sys.stdout.flush()
            ec1[i].append(j)
            placed.add(j)

    for j in range(n2):
        found = False
        for i in ec1:
            # print i,j,
            sys.stdout.flush()
            if matchtest(list1[i],list2[j]):
                # print "match"
                # sys.stdout.flush()
                ec2[i].append(j)
                found = True
                break
            # else:
            #   print "no match"
            #   sys.stdout.flush()
        if not found:
            raise StopIteration

    # print ec1
    # print ec2

    anypairs = []
    args = []
    for i in ec1:
        if len(ec1[i]) != len(ec2[i]):
            raise StopIteration
        args.append(iterpairs(ec1[i],ec2[i]))
        anypairs.extend(zip(ec1[i],ec2[i]))

    # print anypairs

    # in some cases, we can short cut by just trying some equiv class
    # consistent matching, as there is no topology that matters,
    # rather than invoking the expensive product enumeration:

    yield list(map(lambda t: list1[t[0]],anypairs)),list(map(lambda t: list2[t[1]],anypairs)) 
    
    for mp in itertools.product(*args):
        pairs = reduce(lambda a,b: list(a) + list(b), mp)
        yield list(map(lambda t: list1[t[0]],pairs)),list(map(lambda t: list2[t[1]],pairs))

def itergenmatchings(items1,items2,matchtest):

    list1 = list(items1)
    n1 = len(list1)
    list2 = list(items2)
    n2 = len(list2)

    if n1 != n2:
        raise StopIteration

    edges = defaultdict(set)
    inedges = defaultdict(set)
    for i,i1 in enumerate(list1):
        for j,i2 in enumerate(list2):
            if matchtest(i1,i2):
                edges[i].add(j)
                inedges[j].add(i)

    # for i in edges:
    #     print "%s:"%i," ".join(map(str,sorted(edges[i])))

    outdegree = defaultdict(int)
    indegree = defaultdict(int)
    for i in range(n1):
        outdegree[i] = len(edges[i])
    for j in range(n2):
        indegree[j] = len(inedges[j])

    startpairs = []
    outdegreezero = 0
    for i in range(n1):
        if outdegree[i] == 0:
            raise StopIteration
        if outdegree[i] == 1:
            j = next(iter(edges[i]))
            startpairs.append((i,j))
    
    indegreezero = 0
    for j in range(n2):
        if indegree[j] == 0:
            raise StopIteration
        if indegree[j] == 1:
            i = next(iter(inedges[j]))
            startpairs.append((i,j))

    startpairs = list(set(startpairs))

    # print startpairs

    if len(startpairs) != len(set(map(itemgetter(0),startpairs))):
        raise StopIteration
    if len(startpairs) != len(set(map(itemgetter(1),startpairs))):
        raise StopIteration

    startn1set = set(range(n1)) - set(map(itemgetter(0),startpairs))
    startn2set = set(range(n2)) - set(map(itemgetter(1),startpairs))

    # print startn1set,startn2set

    if len(startn1set) == 0:
        l1 = map(lambda t: list1[t[0]],startpairs)
        l2 = map(lambda t: list2[t[1]],startpairs)
        yield list(l1),list(l2)

    start = (startpairs,startn1set,startn2set)
    partialsolutions = [start]
    while len(partialsolutions) > 0:
        pairs,tochoose1,tochoose2 = partialsolutions.pop()
        # print pairs,tochoose1,tochoose2
        if len(tochoose1) == 0:
            assert len(tochoose2) == 0
            l1 = map(lambda t: list1[t[0]],pairs)
            l2 = map(lambda t: list2[t[1]],pairs)
            yield list(l1),list(l2)
        else:
            i1 = sorted(tochoose1,key=outdegree.get)[0]
            newtochoose1 = set(filter(lambda i: i != i1,tochoose1))
            for i2 in (edges[i1]&tochoose2):
                newtochoose2 = set(filter(lambda i: i != i2,tochoose2))
                partialsolutions.append((pairs+[(i1,i2)],newtochoose1,newtochoose2))

def itergenmaximalmatchings(items1,items2,matchtest):

    list1 = list(items1)
    n1 = len(list1)
    list2 = list(items2)
    n2 = len(list2)

    edges = defaultdict(set)
    inedges = defaultdict(set)
    for i,i1 in enumerate(list1):
        for j,i2 in enumerate(list2):
            if matchtest(i1,i2):
                edges[i].add(j)
                inedges[j].add(i)

    # for i in edges:
    #     print "%s:"%i," ".join(map(str,sorted(edges[i])))
    # sys.stdout.flush()

    outdegree = defaultdict(int)
    for i in range(n1):
        outdegree[i] = len(edges[i])

    startn1set = set(range(n1))
    startn2set = set(range(n2))

    # print startn1set,startn2set

    start = ([],startn1set,startn2set)
    partialsolutions = [start]
    while len(partialsolutions) > 0:
        pairs,tochoose1,tochoose2 = partialsolutions.pop(0)
        # print pairs, tochoose1, tochoose2
        if len(tochoose1) == 0:
            maximal = True
            for i1 in map(lambda t: t[0],filter(lambda t: t[1] == None,pairs)):
                if len(edges[i1] & tochoose2) > 0:
                    maximal = False
                    break
            # print "!!",pairs,tochoose1,tochoose2,maximal
            if maximal:
                l1 = map(lambda t: list1[t[0]],filter(lambda t: t[1] != None,pairs))
                l2 = map(lambda t: list2[t[1]],filter(lambda t: t[1] != None,pairs))
                tc1 = map(lambda i: list1[i],map(lambda t: t[0],filter(lambda t: t[1] == None,pairs)))
                tc2 = map(lambda i: list2[i],tochoose2)
                yield list(l1),list(l2),list(tc1),list(tc2)
        else:
            i1 = sorted(tochoose1,outdegree.get)[0]
            newtochoose1 = set(filter(lambda i: i != i1,tochoose1))
            for i2 in (edges[i1]&tochoose2):
                newtochoose2 = set(filter(lambda i: i != i2,tochoose2))
                partialsolutions.append((pairs+[(i1,i2)],newtochoose1,newtochoose2))
            partialsolutions.append((pairs+[(i1,None)],newtochoose1,set(tochoose2)))
            partialsolutions.sort(key=lambda t: (-len(t[0]),sum(1 for _ in filter(lambda tt: tt[1] == None,t[0]))))

def testperm(l):
    print("Permutations of",','.join(map(str,l)))
    for p in permutations(l):
        print(','.join(map(str,p)))

def testselect(l,k):
    print("Selections (size %d) of"%k,','.join(map(str,l)))
    for p in select(l,k):
        print(','.join(map(str,p)))

def testchoose(l,k):
    print("Choices (size %d) of"%k,','.join(map(str,l)))
    for p in choose(l,k):
        print(','.join(map(str,p)))

def testprod(*args,**kw):
    print("Products of", ", ".join(map(str,args)))
    for p in product(*args,**kw):
        print(p)

def testiterecmatch(*args):
    print("Iter EC matching of")
    print("  ",args[0])
    print("and")
    print("  ",args[1])
    for i,p in enumerate(iterecmatchings(*args)):
        print(i+1,p)

def testitergenmatch(*args,**kw):
    print("Iter general matching of")
    print("  ",args[0])
    print("and")
    print("  ",args[1])
    for i,p in enumerate(itergenmatchings(*args,**kw)):
        print(i+1,p)

if __name__ == "__main__":

    testprod('abc','def','ijk',[1,2,3,4],
             accumulator=tuple_accumulator)
    testprod('abc','def','ijk',[1,2,3,4],
             accumulator=concat_accumulator)

    testperm([1,4,6,10])
    testperm(range(5))
    testperm('abcabc')

    testselect([1,4,6,10],4)
    testselect(range(5),1)
    testselect('abcabc',3)

    testchoose([1,4,6,10],2)
    testchoose([1,4,6,10],3)
    testchoose(range(5),1)
    testchoose('abcabc',3)

    testiterecmatch(["1.1","1.2","1.3","2.1","2.2","3.1","3.2","3.3","3.4"],
                    ["2.1","2.2","1.1","3.1","3.2","1.2","3.3","3.4","1.3"],
                    lambda x,y: int(float(x))==int(float(y)))

    testitergenmatch(["1.1","1.2","1.3","2.1","2.2","3.1","3.2","3.3","3.4"],
                     ["2.1","2.2","1.1","3.1","3.2","1.2","3.3","3.4","1.3"],
                     lambda x,y: int(float(x))==int(float(y)))

