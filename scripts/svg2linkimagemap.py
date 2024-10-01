#!/bin/env python2

"""
This script converts a subset of SVG into an HTML imagemap

Note *subset*.  It only handles <path> elements, for which it only pays
attention to the M and L commands.  Futher, it only notices the "translate"
transform.

It was written to generate the examples in the documentation for maphilight,
and thus is very squarely aimed at handling several SVG maps from wikipedia.
It *assumes* that all the <path>s it will need are inside a <g>.  Any <path>
outside of a <g> will be ignored.

It takes several possible arguments, in the form:
$ svn2imagemap.py FILENAME [x y [group1 group2 ... groupN]]

FILENAME must be the name of an SVG file.  All other arguments are optional.

x and y, if present, are the dimensions of the image you'll be creating from
the SVG.  If not present, it assumes the values of the width and height
attributes in the SVG file.

group1 through groupN are group ids.  If only want particular groups used,
enter their ids here and all others will be ignored.
"""

"""
This is modified from https://davidlynch.org/blog/2008/03/creating-an-image-map-from-svg/
This code looks for monosaccharide names (Glc, Gal, etc) and paths in our produced svgs and creates an html image map from them.
We use the image map to create YOLO-format bounding boxes later.
Why? Because this is the first svg path -> png coordinate parser I found.
"""

import os
import re
import sys
import xml.dom.minidom

import parse_path

if len(sys.argv) == 1:
    sys.exit("svn2imagemap.py FILENAME [x y [group1 group2 ... groupN]]")
if not os.path.exists(sys.argv[1]):
    sys.exit("Input file does not exist")
x, y, groups = None, None, None
if len(sys.argv) >= 3:
    x = float(sys.argv[2])
    y = float(sys.argv[3])
    if len(sys.argv) > 3:
        groups = sys.argv[4:]


svg_file = xml.dom.minidom.parse(sys.argv[1])
svg = svg_file.getElementsByTagName('svg')[0]
#print svg
#### viewBox="0 0 342 160"
#print ("hello", svg.getAttribute('viewBox'))
svg_viewbox = svg.getAttribute('viewBox').split()
#print(svg_viewbox)
svg_width = svg_viewbox[2]
svg_height = svg_viewbox[3] 
raw_width = float(svg_width)
raw_height = float(svg_height)
width_ratio = x and (x / raw_width) or 1
height_ratio = y and (y / raw_height) or 1


svg_file = xml.dom.minidom.parse(sys.argv[1])
svg = svg_file.getElementsByTagName('svg')[0]
#print svg

#### viewBox="0 0 342 160"

#print ("hello", svg.getAttribute('viewBox'))
svg_viewbox = svg.getAttribute('viewBox').split()
#print(svg_viewbox)
svg_width = svg_viewbox[2]
svg_height = svg_viewbox[3] 
raw_width = float(svg_width)
raw_height = float(svg_height)
width_ratio = x and (x / raw_width) or 1
height_ratio = y and (y / raw_height) or 1

if groups:
    #print "groups"
    elements = [g for g in svg.getElementsByTagName('g') if (g.hasAttribute('ID') and g.getAttribute('ID') in groups)]
    elements.extend([p for p in svg.getElementsByTagName('path') if (p.hasAttribute('ID') and p.getAttribute('ID') in groups)])
    #print elements
else:
    #print "not groups"
    elements = svg.getElementsByTagName('g')

parsed_groups = {}
for e in elements:
    pointset_count = 0
    if e.nodeName == 'g':
        #print "nodeName g"
        for node in e.childNodes:
            #print node.nodeName
            if node.nodeName == 'defs':
                clipPaths = node.childNodes
                for clipPath in clipPaths:
                    if clipPath.nodeName == 'clipPath':
                        #print "clipPath found"
                        clipPathID = clipPath.getAttribute('id')
                        #print clipPathID
                        svgpaths = clipPath.childNodes
                        paths = []
                        for path in svgpaths:
                            #print path
                            if path.nodeName == 'path':
                                #print "path" 
                                points = parse_path.get_points(path.getAttribute('d'))
                                #print points
                                for pointset in points:
                                    paths.append([clipPathID, pointset])
                                    #print paths, pointset_count
                                    pointset_count += 1
                        parsed_groups[clipPathID] = paths
    else:
        #print "nodeName not g"
        points = parse_path.get_points(e.getAttribute('d'))
        for pointset in points:
            paths.append([e.getAttribute('ID'), pointset])
    if e.hasAttribute('transform'):
        print e.getAttribute('ID'), e.getAttribute('transform')
        for transform in re.findall(r'(\w+)\((-?\d+.?\d*),(-?\d+.?\d*)\)', e.getAttribute('transform')):
            if transform[0] == 'translate':
                x_shift = float(transform[1])
                y_shift = float(transform[2])
                for path in paths:
                    path[1] = [(p[0] + x_shift, p[1] + y_shift) for p in path[1]]

groups = {}
for e in elements:
    if e.hasAttribute('ID'):
        data_type = e.getAttribute("data.type")
        if data_type == 'Monosaccharide':
            gid = e.getAttribute("ID")
            name = e.getAttribute("data.residueName") 
            stylestring = e.firstChild.getAttribute("style")
            stylestring = stylestring.split("(")[1]
            stylestring = stylestring.split(")")[0]
            pathname = stylestring.strip("#")
            groups[gid] = []           
            groups[gid].append(str(name))
            for i in parsed_groups[pathname][0][1]:
                groups[gid].append(i)
            lenth = int(e.childNodes[1].getAttribute("height"))
            cx = int(e.childNodes[1].getAttribute("x")) + lenth/2
            cy = int(e.childNodes[1].getAttribute("y")) + lenth/2
            groups[gid].append((cx,cy))
            groups[gid].append(lenth)
        cx = 0
        cy = 0
        lenth = 0
        if data_type == 'Linkage':
            gid = e.getAttribute("ID")     
            groups[gid] = []

out = []
for g in groups:
    if g[0] == 'l':
        t = g.split(':')[1].split(',')
        out.append('l\t'+t[0]+'\t'+t[1])
    if g[0] == 'r':
        tmp = 'm\t'+g.split(':')[-1]+'\t'
        for p in groups[g]:
            if type(p) == tuple:
                tmp += str(int(p[0]*width_ratio)) +',' + str(int(p[1]*height_ratio))+'\t'
            else:
                tmp += str(p) + '\t'
        out.append(tmp[:-1])
out.sort()        
outfile = open(sys.argv[1].replace('.svg', '_map.txt'), 'w')
outfile.write('\n'.join(out))
