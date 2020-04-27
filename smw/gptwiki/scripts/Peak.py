import math
import numpy as np

window = 3

def xic (dreader, column_name, evgrt = {},second_xic = False):
    x = []
    y = []
#    print 'column_name=',column_name
    for row in dreader:
        for pid in evgrt:
            if column_name.find(pid) != -1:
		if float(row['rt']) > (float(evgrt[pid]) - window) and float(row['rt']) < (float(evgrt[pid]) + window):
        	    x.append(float(row['rt']))
        	    y.append(float(row[column_name]))
    return x, y

def centroid (x = [], y = []):
    index = np.argmax(y)
    rt = x[index]
    start_time = rt - 1.5
    end_time = rt + 1.5
    newx = []
    newy = []
    for n, i in enumerate(x):
	if i >= start_time and i <= end_time:
	    newx.append(x[n])
	    newy.append(y[n])
    return newx, newy    

def get_startpoint(xlist = [], ylist = []):
    ymax = max(ylist)
    rt = xlist[np.argmax(ylist)]
    dig = len(str(int(ymax)))
    x = np.array(xlist)
    y = np.array(ylist)
    slopes = {}
    starts = {}
    ends = {}
    for n in range(0,len(x)-4):
        xs = x[n:n+5]
        ys = y[n:n+5]
        slope = ( ((np.mean(xs)*np.mean(ys))-np.mean(xs*ys)) /
           ((np.mean(xs)*np.mean(xs))-np.mean(xs*xs)) )
        slopes[x[n]] = round(slope/math.pow(10,int(dig)-1),4)
    for k,v in slopes.items():
	if k < rt and v < 1 and ylist[xlist.index(k)] < ymax/10:
	    starts[k] = v
	    if starts[k] < 0:
		starts[k] == 0
	elif k > rt and v > -1 and ylist[xlist.index(k)] < ymax/10:
	    ends[k] = v
	    if ends[k] > 0:
		ends[k] == 0
	else: 
	    continue
    if starts and ends:
        startpoint = max(starts)    
	endpoint = min(ends)
    else:
	return [],[]
    xl = []
    yl = []
    for n, i in enumerate(xlist):
        if i >= startpoint and i <= endpoint:
            xl.append(xlist[n])
            yl.append(ylist[n])
    return xl, yl
