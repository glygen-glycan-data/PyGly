import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')


from collections import defaultdict
from collections import OrderedDict
import csv,sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from regression import IteratedRobustLinearRegression as IRLR

if len(sys.argv) < 4:
    print 'pleae enter the table.tsv, library file, fdr file name!'
    exit(1)

dictrts = {}
dictscores = {}

tablefile = csv.DictReader(open(sys.argv[1]), delimiter='\t')

for row in tablefile:
    tg_id = row['peptide_group_label']
    rt = round(float(row['RT'])/60,3)
#    score = float(row['xx_lda_prelim_score'])
    score = float(row['main_var_xx_swath_prelim_score'])
    if tg_id not in dictrts:
	dictrts[tg_id] = rt
	dictscores[tg_id] = score
    else:
	if score > dictscores[tg_id]:
       	    dictscores[tg_id] = score
	    dictrts[tg_id] = rt
	else:
	    continue

threshold = float(open(sys.argv[3]).read().split()[0])

goodscores = {}

for tgid in dictscores:
    if dictscores[tgid] >= threshold:
        goodscores[tgid] = dictscores[tgid] 
    else:
	continue

x = []
y = []
points = []

dictassayrts = {}

libfile = csv.DictReader(open(sys.argv[2]), delimiter='\t')
for row in libfile:
    tg_id = row['transition_group_id']
    assay_rt = float(row['Tr_recalibrated'])
    if tg_id in goodscores and tg_id not in dictassayrts:
	    dictassayrts[tg_id] = assay_rt
	    x.append(assay_rt)
	    y.append(dictrts[tg_id])
	    points.append((assay_rt,dictrts[tg_id]))
    else:
	continue

assayrt_x = np.array(x)
rt_y = np.array(y)
slope, intercept, r_value, p_value, std_err = stats.linregress(assayrt_x,rt_y)

irlr = IRLR(minpoints=5,max_rvalue=0.99)

points1 = irlr.fit(points)

final_x = []
final_y = []
final_points = []
for tg in dictassayrts:
    nrt = dictassayrts[tg]
    rt = dictrts[tg]
    if (nrt, rt) in points1:
	final_x.append(nrt)
	final_y.append(rt)
finalnrt_x = np.array(final_x)
finalrt_y = np.array(final_y)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(finalnrt_x,finalrt_y)

plt.plot(assayrt_x, rt_y, '.')
plt.plot(assayrt_x, intercept + slope*assayrt_x, 'r')
plt.plot(finalnrt_x, intercept1 + slope1*finalnrt_x)
plt.xlabel('assay_rt')
plt.ylabel('rt(min)')
plt.savefig(sys.argv[1].split('.')[0]+'.png') 
plt.show()
plt.close()

print '<?xml version="1.0" encoding="UTF-8"?>'
print '<TrafoXML version="1.0" xsi:noNamespaceSchemaLocation="https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/TrafoXML_1_0.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
print '  <Transformation name="linear">'
print '    <Pairs count="2">'
print '      <Pair from="0" to="'+str(0-intercept1/slope1)+'"/>'
print '      <Pair from="1" to="'+str(1/(slope1*60)-intercept1/slope1)+'"/>'
print '    </Pairs>'
print '  </Transformation>'
print '</TrafoXML>'
