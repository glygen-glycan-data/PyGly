# -*- coding: utf-8 -*-
"""
use heuristic mono finding colour ranges to make ranges of blue/green/red/etc
"""
import sys
import cv2
import numpy as np
import random

# old_color_range_dict = {
#     "black_lower" : np.array([0,0,0]),
#     "black_upper" : np.array([180,171,235]),
#     "red_lower_l" : np.array([0,18,120]),
#     "red_upper_l" : np.array([19,255,255]),
#     "yellow_lower" : np.array([20,62,104]),
#     "yellow_upper" : np.array([36,255,255]),
#     "green_lower" : np.array([37,40,70]),
#     "green_upper" : np.array([77,255,255]),
#     "blue_lower" : np.array([87,40,60]),
#     "blue_upper" : np.array([134,255,255]),
#     "purple_lower" : np.array([135,20,40]),
#     "purple_upper" : np.array([169,255,225]),
#     "red_lower_h" : np.array([156,40,120]),
#     "red_upper_h" : np.array([180,255,255]),
#     }

color_range_dict = {
    "black_lower" : np.array([0,0,0]),
    "black_upper" : np.array([180,171,235]),
    "red_lower_l" : np.array([0,98,120]),
    "red_upper_l" : np.array([9,255,255]),
    "yellow_lower" : np.array([18,33,222]),
    "yellow_upper" : np.array([30,255,255]),
    "green_lower" : np.array([37,30,70]),
    "green_upper" : np.array([73,255,255]),
    "blue_lower" : np.array([99,123,77]),
    "blue_upper" : np.array([120,255,255]),
    "purple_lower" : np.array([128,58,74]),
    "purple_upper" : np.array([163,255,225]),
    "red_lower_h" : np.array([170,111,135]),
    "red_upper_h" : np.array([180,255,255]),
    "light_blue_lower": np.array([85,41,201]),
    "light_blue_upper": np.array([108,121,255]),
    }

img_file = sys.argv[1]


image = cv2.imread(img_file)
hsv=cv2.cvtColor(image,cv2.COLOR_BGR2HSV)

yellow_mask = cv2.inRange(hsv, color_range_dict['yellow_lower'], color_range_dict['yellow_upper'])
purple_mask = cv2.inRange(hsv, color_range_dict['purple_lower'], color_range_dict['purple_upper'])
red_mask_l = cv2.inRange(hsv, color_range_dict['red_lower_l'], color_range_dict['red_upper_l'])
red_mask_h = cv2.inRange(hsv, color_range_dict['red_lower_h'], color_range_dict['red_upper_h'])
red_mask = red_mask_l + red_mask_h
green_mask = cv2.inRange(hsv, color_range_dict['green_lower'], color_range_dict['green_upper'])
blue_mask = cv2.inRange(hsv, color_range_dict['blue_lower'], color_range_dict['blue_upper'])
black_mask = cv2.inRange(hsv, color_range_dict['black_lower'], color_range_dict['black_upper'])

# store these mask into array
mask_array = (red_mask, yellow_mask, green_mask, blue_mask, purple_mask, black_mask)
mask_array_name = ("red_mask", "yellow_mask", "green_mask", "blue_mask", "purple_mask", "black_mask")
mask_dict = dict(zip(mask_array_name, mask_array))


replacement_color_dict = {}

for mask_name in mask_dict.keys():
    color = mask_name.split("_")[0]
    if color == "black":
        continue
    hsv_range_dict = {"0": [],
                      "1": [],
                      "2": []}
    
    for name in color_range_dict.keys():
        if name.startswith(color):
            if "lower" in name:
                lowername = name
                uppername = lowername.replace("lower","upper")
                for idx in range(0,3):
                    low = color_range_dict[lowername][idx]
                    high = color_range_dict[uppername][idx]
                    valuerange = range(low, high)
                    [ hsv_range_dict[str(idx)].append(r) for r in valuerange ]
    
    h = random.choice(hsv_range_dict["0"])
    s = random.choice(hsv_range_dict["1"])
    v = random.choice(hsv_range_dict["2"])
    replacement_color_dict[color] = (h,s,v)
#print(replacement_color_dict)   

for color, value in replacement_color_dict.items():
    for name, mask in mask_dict.items():
        if name.startswith(color):
            hsv[mask>0] = value

img = cv2.cvtColor(hsv,cv2.COLOR_HSV2BGR)
blur = random.choice([True,False])
if blur:
    img = cv2.blur(img, (4,4))
# cv2.imshow("test",img)
# cv2.waitKey(0)
# cv2.destroyAllWindows
cv2.imwrite(img_file, img)
