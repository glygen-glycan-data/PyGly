
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
	ind = np.argsort(apoints[:,0])
	x = apoints[ind,0]
	y = apoints[ind,1]
	invind = np.argsort(ind)
        return x,y,invind

    def check(self,points):
	if len(points) < self._minpoints:
	    raise ValueError("Too few points for linear regression")

    def evaluate(self,params,points):
	x,y,ind = self.normalize(points)
	yfit = params['slope']*x+params['intercept']
	return dict(yfit=yfit[ind],residuals=yfit[ind]-y[ind])

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
	super(IteratedRobustLinearRegression,self).__init__(**kw)
	# initialize simple linear regression we will call repeatedly.
	self._slr = SimpleLinearRegression()

    def fit(self,points):
	self.check(points)
	points1 = copy.copy(points)
	while True:
	    params=self._slr.fit(points1)
	    result=self._slr.evalute(params,points1)
	    # decide whether to stop or not...
	    # Use r, residuals, etc.
	    if stop:
		break
	    # Check it is reasonable to keep going, if not raise RuntimeError("message")
	    if not resonable:
		raise RuntimeError("Can't do the fit!")
	    # manipulate points1 to remove one pair of values
	    points1.remove((None,None))
	
	# Once out of loop, result has values from final SLR fit
	return result

if __name__ == "__main__":
    import sys, os

    lr = SimpleLinearRegression()
    # lr = IteratedRobustLinearRegression(outlier_threshold=2)
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
