import sys
from html.parser import HTMLParser
import numpy as np
import cv2
import os

import boundingboxes

'''
This uses the boundingboxes class file used in the image extractor.
Parses the image map from the HTML to create a bounding box.
We pad the borders of the box in this case because svg coordinates are very tight to the object,
and YOLO consistently returns smaller detections than the boxes we train on (for some reason).
This is optional, depends on your use case.
'''

# little parser, to get the coordinates and name of the monosaccharide, and make coordinates into a python list
# includes a function to directly create a list of the mapped monosaccharides; we call that later

image_file = sys.argv[1]    
image = cv2.imread(image_file)
image_map_file = sys.argv[2]
image_map = open(image_map_file)

method = sys.argv[3]
outfile = sys.argv[4]

# this is set up for monosaccharides, as coded in the image extractor
# change this if you're svg parsing a different thing.
class_dict = {
    "-" : 0,
    "|" : 1,
    "/" : 2,
    "\\" : 3
    }

monos={}
links=[]
for i in image_map:
    seq=i.split('\t')
    if seq[0] == 'm':
        monos[int(seq[1])] = seq[2:]
    elif seq[0] == 'l':
        links.append(seq[1:])
# creating the bounding box specs (x min, width, y min, height, class) from the html map and the corresponding image
f = open(outfile,"a")
tmp = ''
for link in links:
    lenth = int(monos[int(link[0])][-1])
    x_list = []
    y_list = []
    coordC = []
    for l in link:
        for coord in monos[int(l)][1:-1]:
            [x,y] = coord.split(',')
            x_list.append(int(x))
            y_list.append(int(y))
        coordC.append((x_list[-1],y_list[-1]))

    ori = ''
    ax = int(coordC[0][0])
    ay = int(coordC[0][1])
    bx = int(coordC[1][0])
    by = int(coordC[1][1])
    if ax == bx:
        ori = class_dict['|']
    elif ay == by:
        ori = class_dict['-']
    elif (ax - bx)*(ay - by) < 0:
        ori = class_dict['/']
    else:
        ori = class_dict["\\"]

    if method == 'o':
        x_min = np.min(x_list)
        y_min = np.min(y_list)
        x_max = np.max(x_list)
        y_max = np.max(y_list)
        w = x_max - x_min
        h = y_max - y_min
    elif method == 'c':
        x_min = np.min([ax,bx])
        y_min = np.min([ay,by])
        x_max = np.max([ax,bx])
        y_max = np.max([ay,by])
        w = x_max - x_min
        h = y_max - y_min
        if h <= 10:
            y_min -= lenth/2
            h = lenth
        if w <= 10:
            x_min -= lenth/2
            w = lenth
    else:
        print('wrong method')
        sys.exit[1]
    box = boundingboxes.Training(image, x=x_min, width=w, y=y_min, height=h,class_=0)
    box.corner_to_center()
    if method != 'c':
        box.fix_pad_borders(10)
    box.abs_to_rel()
    boxlist = box.to_list()
    for x in boxlist:
        tmp += str(x) + " "
    tmp = tmp[:-1]
    tmp += "\n"
f.write(tmp)
f.close()
