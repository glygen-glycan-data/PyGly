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
class GetMapCoords(HTMLParser):
    def handle_starttag(self, tag, attrs):
        if tag == 'area':
            for attr in attrs:
                if attr[0] == 'title':
                    mono = attr[1]
                if attr[0] == 'coords':
                    coord_list = self.make_coordinate_list(attr[1])
            self.mono_list.append([mono, coord_list])
    def make_coordinate_list(self,coords):
        coord_list = coords.split()
        coord_list = [ coord.strip(",") for coord in coord_list ]
        coord_list = [ coord.split(",") for coord in coord_list ]
        return coord_list
    def make_mono_list(self,text):
        self.mono_list = []
        self.feed(text)
        return self.mono_list
                
image_file = sys.argv[1]    
image = cv2.imread(image_file)
image_map_file = sys.argv[2]
image_map_html = open(image_map_file)

outfile = sys.argv[3]

# this is set up for monosaccharides, as coded in the image extractor
# change this if you're svg parsing a different thing.
class_dict = {
    "GlcNAc" : 0,
    "NeuAc" : 1,
    "Fuc" : 2,
    "Man" : 3,
    "GalNAc" : 4,
    "Gal" : 5,
    "Glc" : 6,
    "NeuGc" : 7,
    }

parser = GetMapCoords()

monos = parser.make_mono_list(image_map_html.read())

parser.close()
image_map_html.close()

# creating the bounding box specs (x min, width, y min, height, class) from the html map and the corresponding image
f = open(outfile,"a")
for mono in monos:
    x_list = []
    y_list = []
    name = mono[0]
    try:
        class_ = class_dict[name]
    except KeyError:
        os.remove(image_file)
        f.close()
        os.remove(outfile)
        break
    for coord in mono[1]:
        [x,y] = coord
        x_list.append(int(x))
        y_list.append(int(y))
    x_min = np.min(x_list)
    y_min = np.min(y_list)
    x_max = np.max(x_list)
    y_max = np.max(y_list)
    w = x_max - x_min
    h = y_max - y_min
    box = boundingboxes.Training(image, x=x_min, width=w, y=y_min, height=h, class_=class_)
    box.corner_to_center()
    box.pad_borders()
    box.abs_to_rel()
    boxlist = box.to_list()
    for x in boxlist:
        f.write(str(x) + " ")
    f.write("\n")
f.close()
