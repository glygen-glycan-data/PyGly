#!/bin/env python27

import getwiki
from analysis.regression import SimpleLinearRegression as SLR
from operator import itemgetter

slr = SLR()
points = [(1,2.5),(3,4.3),(2,3.9)]
print map(itemgetter(0),points)
print map(itemgetter(1),points)
params = slr.fit(points)
print params['slope']
print params['intercept']
result = slr.evaluate(params,points)
print result['yfit']
print result['residuals']

