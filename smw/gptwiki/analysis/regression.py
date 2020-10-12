
import numpy as np
import scipy.stats as spstat
import copy

np.seterr(all='raise')

class LinearRegression(object):

    def __init__(self,**kw):
	self._minpoints = kw.get('minpoints',2)
    
    def normalize(self,points):
	self.check(points)
	apoints = np.array(points)
	assert len(apoints.shape) == 1 or (len(apoints.shape) == 2 and apoints.shape[1] == 2)
	if len(apoints.shape) == 2 and apoints.shape[1] == 2:
            ind = np.argsort(apoints[:,0])
            x = apoints[ind,0]
            y = apoints[ind,1]
            invind = np.argsort(ind)
            return x,y,invind
	else:
            ind = np.argsort(apoints)
            x = apoints[ind,0]
            invind = np.argsort(ind)
            return x,invind

    def check(self,points):
	if len(points) < self._minpoints:
	    raise ValueError("Too few points for linear regression")
    
    def evaluate(self,params,points):
        x,y,ind = self.normalize(points)
        yfit = self.y(x)
        return dict(yfit=yfit[ind],residuals=yfit[ind]-y[ind])

    def y(self,params,x):
	return params['slope']*x+params['intercept']

class SimpleLinearRegression(LinearRegression):

    """
        slr = SimpleLinearRegression(corr_coeff=True)
        points = [(1,2),(2,3),(3,4)]
	params = slr.fit(points)
	print params
	print slr.evaluate(params,points)
    """

    def __init__(self,**kw):
	super(SimpleLinearRegression,self).__init__(**kw)

    def fit(self,points):
	x,y,ind = self.normalize(points)
        a, b, r, p, e = spstat.linregress(x,y)
	result = dict(slope=a, intercept=b, r=r)
        return result

class IteratedRobustLinearRegression(LinearRegression):

    def __init__(self,**kw):
	# set parameters, like outlier threshold, etc. see above for example
	self._stopping_rvalue = kw.get('max_rvalue',0.99)
	self._poor_linearity_rvalue = kw.get('min_rvalue',0.5)
	self._maxremoved = kw.get('max_removed',20)
	super(IteratedRobustLinearRegression,self).__init__(**kw)
	# initialize simple linear regression we will call repeatedly.
	self._slr = SimpleLinearRegression()

    def fit(self,points):
	self.check(points)
	points1 = copy.copy(points)

	removed = 0
	while True:
	    try:
	        params=self._slr.fit(points1)
	    except:
		raise RuntimeError('fit failed')

	    result=self._slr.evaluate(params,points1)
	    
	    # decide whether to stop or not...
	    # Use r, residuals, etc.
	    if params['r'] > self._stopping_rvalue:
		break
	    
	    # Check it is reasonable to keep going, if not raise RuntimeError("message")
	    if params['r'] < self._poor_linearity_rvalue:
		raise RuntimeError("Can't do the fit!")
	  
	    index = np.argmax(abs(result['residuals']))

	    # print points1[index],' is removed!' 
	    del points1[index]

	    removed += 1
	    if removed > self._maxremoved:
                raise RuntimeError('Too many points removed!')
            if len(result['residuals']) < self._minpoints:
                raise RuntimeError('Too few points!')

	# Once out of loop, result has values from final SLR fit
	params['retained_points'] = list(points1)
	return params

if __name__ == "__main__":
    import sys, os

   # lr = SimpleLinearRegression()
    lr = IteratedRobustLinearRegression()
    points = [(1,2),(2,3),(3,4)]
    print points
    print lr.normalize(points)
    points = [(2,3),(3,4),(1,2)]
    print points
    print lr.normalize(points)
    params = lr.fit(points)
    print params
    print lr.evaluate(params,points)

    points = [(1,2.5),(2,3.9),(3,4.3)]
    params = lr.fit(points)
    print params
    print lr.evaluate(params,points)

    points = [(1,2.5),(3,4.3),(2,3.9)]
    params = lr.fit(points)
    print params
    print lr.evaluate(params,points)
