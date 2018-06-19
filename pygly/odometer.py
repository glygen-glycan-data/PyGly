#!/usr/bin/env python

import copy

class basic_odometer:
    def __init__(self):
        self.inrange_ = False
        self.total_max_set_ = False
        self.weighted_total_max_set_ = False
        self.n_ = 0
    def init(self):
        self.value_ = copy.deepcopy(self.min_)
        self.inrange_ = True
        self.ord_ = 0;
    def set_size(self,n):
        if n == 0 and self.n_ == 0:
            return
        self.min_ = [0]*n
        self.max_ = [0]*n
        self.value_ = [0]*n
        self.total_weights_ = [1]*n
        self.n_ = n
    def set_min(self,m,i=-1):
        if i == -1:
            self.min_ = [m]*self.n_
        else:
            self.min_[i] = m
    def set_max(self,m,i=-1):
        if i == -1:
            self.max_ = [m]*self.n_
        else:
            self.max_[i] = m
    def set_total_max(self,t):
        self.total_max_ = t
        self.total_max_set_ = True
    def set_weighted_total_max(self,t):
        self.weighted_total_max_ = t
        self.weighted_total_max_set_ = True
    def set_weighted_total_weights(self,w,i=-1):
        if i == -1:
            self.total_weights_ = [w]*self.n_
        else:
            self.total_weights_[i] = w;
    def get_value(self,i):
        return self.value_[i]
    def values(self):
        return self.value_
    def inrange(self):
        return self.inrange_
    def inc(self):
        i = 0
        if i < self.n_:
            self.value_[i]+=1
            weighted_total=self.weighted_sum()
            total=self.sum()
            while (self.value_[i] > self.max_[i]) or \
                      (self.total_max_set_ and total > self.total_max_) or \
                      (self.weighted_total_max_set_ and weighted_total > self.weighted_total_max_):
                self.value_[i] = self.min_[i]
                i += 1
                if i >= self.n_:
                    break
                self.value_[i] += 1
                total=self.sum()
                weighted_total=self.weighted_sum()
            self.ord_ += 1
        if i >= self.n_:
            self.inrange_ = False
    def weighted_sum(self):
        return sum([self.total_weights_[j]*self.value_[j] for j in xrange(0,self.n_)])
    def sum(self):
        return sum(self.value_)
    def write(self,h,sep=''):
        for j in xrange(0,self.n_-1):
            h.write('%d%s'%(self.value_[j],sep))
        h.write('%d'%(self.value_[self.n_-1],))

class composite_odometer:
    def __init__(self):
        self.inrange_ = False
        self.total_max_set_ = False
        self.n_ = 0
    def init(self):
        for o in self.value_:
            o.init()
        self.inrange_ = True
    def set_size(self,n):
        if n == 0 and self.n_ == 0:
            return
        self.value_ = [ basic_odometer() ]
        for i in xrange(1,n):
            self.value_.append(basic_odometer())
        self.total_weights_ = [1]*n
        self.n_ = n
    def set_min(self,m,i=-1,j=-1):
        if i == -1:
            for o in self.value_:
                o.set_min(m,j)
        else:
            self.values_[i].set_min(m,j)
    def set_max(self,m,i=-1,j=-1):
        if i == -1:
            for o in self.value_:
                o.set_max(m,j)
        else:
            self.values_[i].set_max(m,j)
    def set_weighted_total_max(self,t):
        self.total_max_ = t
        self.total_max_set_ = True
    def set_weighted_total_weights(self,w,i=-1):
        if i == -1:
            self.total_weights_ = [w]*self.n_
        else:
            self.total_weights_[i] = w;
    def get_value(self,i):
        return self.value_[i]
    def values(self):
        return self.value_
    def inrange(self):
        return self.inrange_
    def weighted_sum(self):
        return sum([self.total_weights_[j]*self.value_[j].weighted_sum() for j in xrange(0,self.n_)])
    def sum(self):
        return sum([self.value_[j].sum() for j in xrange(0,self.n_)])
    def write(self,h,sep=' ',sep1=''):
        for j in xrange(0,self.n_-1):
            self.value_[j].write(h,sep1)
            h.write(sep)
        self.value_[self.n_-1].write(h,sep1)
    def inc(self):
        i = 0
        if i < self.n_:
            if self.total_max_set_:
                total=self.weighted_sum()
                total -= self.total_weights_[i]*self.value_[i].weighted_sum()
                self.value_[i].set_weighted_total_max(self.total_max_-total)
            self.value_[i].inc()
            # print >>sys.stdout, "1: ",
            # self.write(sys.stdout)
            # print >>sys.stdout, " ",total
            while (not self.value_[i].inrange()):
                self.value_[i].init()
                # print >>sys.stdout, "2: ",
                # self.write(sys.stdout)
                # print >>sys.stdout, " ",total
                i += 1
                if i >= self.n_:
                    break
                if self.total_max_set_:
                    total=self.weighted_sum()
                    total -= self.total_weights_[i]*self.value_[i].weighted_sum()
                    self.value_[i].set_weighted_total_max(self.total_max_-total)
                self.value_[i].inc()
                # print >>sys.stdout, "3: "
                # self.write(sys.stdout)
                # print >>sys.stdout, " ",total
        if i >= self.n_:
            self.inrange_ = False

if __name__ == "__main__":
    import sys

    n = int(sys.argv[1])
    m = int(sys.argv[2])

    o = basic_odometer()
    #o = composite_odometer()
    # o.set_size(2)
    # o.get_value(0).set_size(n)
    #o.get_value(1).set_size(n)

    o.set_size(n)
    o.set_min(0)
    o.set_max(m)
    o.set_max(0,3)
    o.set_weighted_total_max(2*m)
    o.set_total_max(m)
    for i in xrange(0,n):
        o.set_weighted_total_weights(i+1,i)
    #     o.get_value(1).set_weighted_total_weights(i+1,i)
    o.init()
    while o.inrange():
        o.write(sys.stdout)
        print >>sys.stdout
        o.inc()
