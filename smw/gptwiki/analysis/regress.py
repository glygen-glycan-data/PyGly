
import numpy as np
import scipy.stats as spstat
from operator import itemgetter
import copy

class LinearRegression(object):

    def __init__(self,**kw):
	pass
    
    def normalize(self,points):
	apoints = np.array(points)
        ind = np.argsort(apoints[:,0])
        x = apoints[ind,0]
        y = apoints[ind,1]
        invind = np.argsort(ind)
        return x,y,invind
    
    def evaluate(self,params,points):
        x,y,ind = self.normalize(points)
        yfit = params['slope']*x+params['intercept']
        return dict(yfit=yfit[ind],residuals=yfit[ind]-y[ind])

class SimpleLinearRegression(LinearRegression):

    def __init__(self,**kw):
        self._corr_coeff = kw.get('corr_coeff',False)
	super(SimpleLinearRegression,self).__init__(**kw)

    def fit(self,points):
	x,y,ind = self.normalize(points)
        a, b, r, p, e = spstat.linregress(x,y)
	result = dict(slope=a, intercept=b)
	if self._corr_coeff:
	    result['r'] = r
        return result

class IteratedRobustLinearRegression(LinearRegression):

    def __init__(self,**kw):
	self._stopping_rvalue = kw.get('max_rvalue',0.99)
	self._poor_linearity_rvalue = kw.get('min_rvalue',0.5)
	self._max_removed = kw.get('max_removed',20)
	self._min_points = kw.get('min_points',5)
	super(IteratedRobustLinearRegression,self).__init__(**kw)

	self._slr = SimpleLinearRegression(corr_coeff=True)

    def fit(self,points):
	points1 = copy.copy(sorted(points))
	if len(points1) < self._min_points:
	    raise RuntimeError('Not enough points')

	removed = 0
	while True:
	    try:
	        params=self._slr.fit(points1)
	    except:
		raise RuntimeError('fit failed')

	    result=self._slr.evaluate(params,points1)

	    if params['r'] > self._stopping_rvalue:
		break
	    
	    if params['r'] < self._poor_linearity_rvalue:
		raise RuntimeError("Can't do the fit!")
	  
	    index = np.argmax(abs(result['residuals']))

	    print points1[index],' is removed!' 
	    del points1[index]

	    removed += 1
	    if removed > self._max_removed:
                raise RuntimeError('Too many points removed!')
            if len(result['residuals']) < self._min_points:
                raise RuntimeError('Not enough points!')
	
	return params

if __name__ == "__main__":
    import sys, os

    lr = IteratedRobustLinearRegression(min_points=2)
    points = [(1,2),(2,3),(3,4)]
    params = lr.fit(points)
    print params
    print lr.evaluate(params,points)

    points = [(1,2.5),(2,3.9),(3,4.3)]
    params = lr.fit(points)
    print params
    print lr.evaluate(params,points)
