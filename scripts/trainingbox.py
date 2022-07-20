#!/bin/env python3

import cv2, sys

for fn in sys.argv[1:]:

    print(fn)
    # Load image, grayscale, Gaussian blur, Otsu's threshold
    image = cv2.imread(fn)
    original = image.copy()
    origshape = original.shape
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    thresh = cv2.threshold(gray, 254, 255, cv2.THRESH_BINARY_INV)[1]
    
    # Find enclosing boundingbox and crop ROI
    coords = cv2.findNonZero(thresh)
    x,y,w,h = cv2.boundingRect(coords)
    origx,origy,origw,origh = x,y,w,h
    
    padding = 0.03

    ymax = image.shape[0]
    xmax = image.shape[1]

    xmid = x+w/2
    ymid = y+h/2
    
    h += min(20,origh*0.1)
    w += min(20,origw*0.1)

    # h *= (1+padding)
    # w *= (1+padding)

    w = min(w,2*(xmax-xmid),2*xmid)
    h = min(h,2*(ymax-ymid),2*ymid)

    # cv2.rectangle(image, (origx, origy), (origx+origw, origy+origh), (36,255,12), 2)
    # cv2.rectangle(image, (int(xmid-w/2), int(ymid-h/2)), (int(xmid+w/2), int(ymid+h/2)), (36,255,12), 2)

    ymid = ymid/ymax
    h = h/ymax
    xmid = xmid/xmax
    w = w/xmax
    
    wh = open(fn.rsplit('.',1)[0]+'.txt','w')
    wh.write("0 %.6f %.6f %.6f %.6f\n"%(xmid,ymid,w,h))
    wh.close()
    
    # cv2.imshow('thresh', thresh)
    # cv2.imshow('close', close)
    # cv2.imshow('image', image)
    # cv2.imshow('crop', crop)
    # cv2.waitKey()
