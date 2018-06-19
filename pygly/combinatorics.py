#!/bin/env python

__all__ = [ 'permutations', 'product', 'list_accumulator',
            'tuple_accumulator', 'concat_accumulator']

import copy

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
            for p in permutations(x[:i] + x[i + 1:]):
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
    if k > x:
        return
    if k == 0:
        yield []
    for i in range(len(x)):
        for s in select(x[:i] + x[i + 1:],k-1):
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
    if k > x:
        return
    if k == 0:
        yield []
    for i in range(len(x)):
        for s in choose(x[i+1:],k-1):
            yield [x[i]] + list(s)
    
def subset(x):
    for i in range(len(x)+1):
        for s in choose(x,i):
            return s

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

def testperm(l):
    print "Permutations of",','.join(map(str,l))
    for p in permutations(l):
	print ','.join(map(str,p))

def testselect(l,k):
    print "Selections (size %d) of"%k,','.join(map(str,l))
    for p in select(l,k):
	print ','.join(map(str,p))

def testchoose(l,k):
    print "Choices (size %d) of"%k,','.join(map(str,l))
    for p in choose(l,k):
	print ','.join(map(str,p))

def testprod(*args,**kw):
    print "Products of", ", ".join(map(str,args))
    for p in product(*args,**kw):
        print p

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

