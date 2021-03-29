
import numpy as np
import scipy.optimize as spopt
import scipy.stats as spstat
import warnings

warnings.simplefilter("error", spopt.OptimizeWarning)

from operator import itemgetter
import math
import sys
import traceback

np.seterr(all='raise')

class Shape(object):
    """
    Abstract base class for Gaussian, paired half-gaussian, etc. peak shapes
    """
    def __call__(self,x,*p,**kw):
	return self.f(x,*p,**kw)
    def extra(self,*args,**kw):
	pass

class Gaussian(Shape):
    rtp = np.sqrt(2*np.pi)
    def f(self,x,*p,**kw):
	" params are mu, sigma, height "
        return p[2]*np.exp(-(x-p[0])**2/(2*p[1]**2))
    def width(self,*p):
	return 6.0*p[1]
    def height(self,*p):
	return p[2]
    def area(self,*p):
	return p[2]*p[1]*self.rtp
    def start(self,*p):
	return p[0]-3.0*p[1]
    def end(self,*p):
	return p[0]+3.0*p[1]
    def centroid(self,*p):
	return p[0]

class MultiGaussian(Gaussian):
    def f(self,x,*p,**kw):
	"""
          x is dimension (k,n), returns (k,n)
          params are mu, sigma, height_1, ...., height_k
        """
	k = len(self.series)-1
	y = np.zeros(x.shape)
	for i in range(k):
	    if self.series[i+1]-self.series[i] == 0:
		continue
	    if max((x[self.series[i]:self.series[i+1]]-p[0])**2/(2*p[1]**2)) > 100:
		continue
	    y[self.series[i]:self.series[i+1]] = p[i+2]*np.exp(-(x[self.series[i]:self.series[i+1]]-p[0])**2/(2*p[1]**2))
	return y
    def extra(self,*args,**kw):
	self.series = args[0]
    def width(self,*p):
	return 6.0*p[1]
    def height(self,*p):
	return tuple(p[2:])
    def area(self,*p):
	return map(lambda h: h*p[1]*self.rtp, p[2:])
    def start(self,*p):
	return p[0]-3.0*p[1]
    def end(self,*p):
	return p[0]+3.0*p[1]
    def centroid(self,*p):
	return p[0]

class PairedHalfGaussian(Gaussian):
    def f(self,x,*p):
	" params are mu, sigma1, sigma2, height "
	li = (x<p[0]); ri = (x>=p[0])
	y = np.zeros(x.shape)
	y[li] = p[3]*np.exp(-(x[li]-p[0])**2/(2*p[1]**2))
	y[ri] = p[3]*np.exp(-(x[ri]-p[0])**2/(2*p[2]**2))
	return y
    def width(self,*p):
	return 3.0*p[1]+3.0*p[2]
    def height(self,*p):
	return p[3]
    def area(self,*p):
	return (0.5*p[1]+0.5*p[2])*p[3]*self.rtp

class ParamEstimator(object):
    """
    Abstract base class for initial estimates of shape parameters
    """
    def __call__(self,x,y):
	return self.param(x,y)

class GaussianParam(ParamEstimator):
    """ Returns mean, sigma, height """
    def param(self,x,y):
        height = np.max(y[0])
        mean = np.sum(x[0]*y[0])/np.sum(y[0])
        sigma = np.sqrt(np.sum(y[0]*(x[0]-mean)**2)/np.sum(y[0]))
        return mean, sigma, height

class MultiGaussianParam(ParamEstimator):
    """ Returns mean, sigma, height_1, ..., height_k """
    def param(self,x,y):
	heights = []
	means = []
	sigmas = []
	for i in range(len(x)):
	    height = np.max(y[i])
	    heights.append(height)
	    if height <= 0:
		continue
            mean = np.sum(x[i]*y[i])/np.sum(y[i])
            sigma = np.sqrt(np.sum(y[i]*(x[i]-mean)**2)/np.sum(y[i]))
	    means.append(mean)
	    sigmas.append(sigma)
	# print >>sys.stderr, means, sigmas
	param = [np.mean(means),np.mean(sigmas)] + heights
        return tuple(param)

class PairedHalfGaussianParam(ParamEstimator):
    """ Returns mean, sigma1, sigma2, height """
    def param(self,x,y):
        height = np.max(y)
        mean = np.sum(x*y)/np.sum(y)
        sigma = np.sqrt(np.sum(y*(x-mean)**2)/np.sum(y))
        return mean, sigma, sigma, height

class RangeFinderBase(object):
    def __call__(self,x,y,*args,**kwargs):
	return self.range(x,y,*args,**kwargs)

class SlopeBasedRangeFinder(RangeFinderBase):

    def __init__(self,**kw):
	self._slopethresh = kw.get('slope_threshold',1)
	self._shoulderthresh = kw.get('shoulder_threshold',0.1)
	self._peak_radius = kw.get('peak_radius',1.5)

    def range(self,x,y,x0=None,xs=None,xe=None):
	if x0 == None:
	    if xs == None:
	        maxi = np.argmax(y)
	        x0 = x[maxi]
	    else:
		ind = ((x>=xs)&(x<=xe))
		maxy = np.max(y[ind])
		maxi = min(np.where((y==maxy)&ind))
		x0 = x[maxi]
	ind = (x>=(x0-self._peak_radius))&(x<=(x0+self._peak_radius))
	x1 = x[ind]
	y1 = y[ind]
        maxi = np.argmax(y1)
	y1max = y1[maxi]
	x1max = x1[maxi]
	s1 = np.zeros(len(x1))
	# print >>sys.stderr, x1,y1
	for i in range(2,len(x1)-2):
	    xs = x1[i-2:i+2]
	    ys = y1[i-2:i+2]
	    # print >>sys.stderr, xs,ys
	    s1[i] = ( (np.mean(xs)*np.mean(ys)-np.mean(xs*ys)) /
                      (np.mean(xs)*np.mean(xs)-np.mean(xs*xs)) )
	    # scale the slope by the order of magnitude of the max peak in the area...
	    s1[i] = round(s1[i]/math.pow(10,int(math.log10(y1max))),4)
        xs = x1[0]
	xe = x1[-1]
	for i in range(2,len(x1)-2):
	    if (x1[i] < x1max) and s1[i] < self._slopethresh and y1[i] < y1max*self._shoulderthresh:
		xs = max(x1[i],xs)
	    if (x1[i] > x1max) and s1[i] < -self._slopethresh and y1[i] < y1max*self._shoulderthresh:
		xe = min(x1[i],xe)
	return xs,xe

class PeakFitterBase(object):
    def __init__(self,**kw):
	if kw.get('shape'):
	    self._shape = shape
	if kw.get('param'):
	    self._param = param
	if kw.get('range'):
	    self._range = range
    def normalize(self,*points):
	x = []
	y = []
	for i in range(len(points)):
	    spoints = sorted(points[i])
	    x.append(np.array(map(itemgetter(0),spoints)))
	    y.append(np.array(map(itemgetter(1),spoints)))
	return x,y
    def restrict(self,x,y,xs=-1e+20,xe=1e+20):
	x1 = []
	y1 = []
	for i in range(len(x)):
	    ind = ((x[i]>=xs)&(x[i]<=xe))
	    x1.append(x[i][ind])
	    y1.append(y[i][ind])
	return x1,y1
    def flatten(self,al):
	# d is a list of value arrays
	cs = np.array([len(ai) for ai in al])
	cs = np.cumsum(cs)
	ind = [0] + list(cs)
	return np.concatenate(al),ind
    def stack(self,arr,ind):
	st = []
	for i in range(len(ind)-1):
	    st.append(np.array(arr[ind[i]:ind[i+1]]))
	return st
    def fit(self,*points,**kwargs):
	x,y = self.normalize(*points)
	if 'x0' in kwargs:
	    xs,xe = self._range(x[0],y[0],x0=kwargs.get('x0'))
	elif 'xs' in kwargs:
	    xs,xe = self._range(x[0],y[0],xs=kwargs.get('xs'),xe=kwargs.get('xe'))
	else:
	    xs,xe = self._range(x[0],y[0])
	x,y = self.restrict(x,y,xs,xe)
	p0 = self._param(x,y)
	xall,xind = self.flatten(x)
	yall,yind = self.flatten(y)
	assert xind == yind
	self._shape.extra(xind)
	try:
	    p1,pcov = spopt.curve_fit(self._shape,xall,yall,p0)
	except FloatingPointError:
	    raise
	except RuntimeError:
	    raise
	except spopt.minpack.error, e:
	    raise FloatingPointError()
	except spopt.OptimizeWarning, e:
	    raise FloatingPointError()
	p1 = tuple(p1)
	xs = self._shape.start(*p1)
	xe = self._shape.end(*p1)
	return xs,xe,p1
    def evaluate(self,xs,xe,param,*points):
	x,y = self.normalize(*points)
	x,y = self.restrict(x,y,xs,xe)
	width = self._shape.width(*param)
	r = []
	for i in range(len(x)):
	    self._shape.extra([0,len(x[i])])
	    yifit = self._shape(x[i],*param)
	    try:
	        slope, intercept, r_value, p_value, std_err = spstat.linregress(y[i],yifit)
	        r.append(r_value)
	    except FloatingPointError:
		r.append(None)
	height = self._shape.height(*param)
	area = self._shape.area(*param)
	return {'width': width, 'height': height, 'area': area, 'r': r}
    def curve(self,xs,xe,param,*points):
	x,y = self.normalize(*points)
	x,y = self.restrict(x,y,xs,xe)
	xfit = []
	for i in range(len(x)):
	    xfit.append(np.arange(min(x[i]),max(x[i]),0.01))
	xallfit,xind = self.flatten(xfit)
	self._shape.extra(xind)
	yallfit = self._shape(xallfit,*param)
	yfit = self.stack(yallfit,xind)
	return x,y,xfit,yfit

class GaussianPeakFit(PeakFitterBase):
    _shape = Gaussian()
    _param = GaussianParam()
    _range = SlopeBasedRangeFinder()

class PairedHalfGaussianPeakFit(PeakFitterBase):
    _shape = PairedHalfGaussian()
    _param = PairedHalfGaussianParam()
    _range = SlopeBasedRangeFinder()

class MultiGaussianPeakFit(PeakFitterBase):
    _shape = MultiGaussian()
    _param = MultiGaussianParam()
    _range = SlopeBasedRangeFinder()

class SimpleLinearRegression(object):
    def fit(self,points):
	spoints = sorted(points)
        x = np.array(map(itemgetter(0),spoints))
        y = np.array(map(itemgetter(1),spoints))
	a, b, r, p, e = spstat.linregress(x,y)
	return dict(slope=a, intercept=b, r=r)

if __name__ == "__main__":

    import sys, json

    series = []
    titles = []
    for fn in sys.argv[1:]:
	chrom = json.loads(open(fn).read())
        peaks = [ (p['rt'],p['int']) for p in chrom['peaks'] ]
	series.append(peaks)
	titles.append(chrom['title'].split()[0])

    ga = MultiGaussianPeakFit()
    xs,xe,pa = ga.fitpeak(*series)
    print xs,xe
    print pa

    x0,y0,xfit,yfit = ga.curve(xs,xe,pa,*series)

    # metrics = ga.evaluate(peaks,xs,xe,pa)
    # print xs,xe,pa,metrics

    from pylab import *
    colors = ['b','r','c','m','y']
    for i in range(min(6,len(x0))):
	plot(x0[i],y0[i],colors[i]+'.')
    for i in range(min(6,len(xfit))):
    	plot(xfit[i],yfit[i],colors[i]+'--',label=titles[i])
    legend()
    show()
    
