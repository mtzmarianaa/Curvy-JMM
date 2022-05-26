import numpy as np
import cv2 as cv
import os

os.chdir(os.path.dirname(os.path.abspath(__file__))) 

im = cv.imread('ommatidium_diagram2.jpg')
cv.imshow('image',im)
imgray = cv.cvtColor(im, cv.COLOR_BGR2GRAY)
ret, thresh = cv.threshold(imgray, 127, 255, 0)
contours, hierarchy = cv.findContours(thresh, cv.RETR_TREE, cv.CHAIN_APPROX_SIMPLE)

list_points = np.array( contours[1] )
for k in range(2, 10):
    np.append(list_points, contours[k] )
    
print( list_points.shape  )
#print(list_points)

# Drawing the contours on top of the original image
cnt_vec = np.arange(-9 ,0)
for cnt_i in cnt_vec: 
    cnt = contours[cnt_i]
    cv.drawContours(im, [cnt], 0, (255, 255, 0), 3 )
cv.imshow('img', im)
cv.imwrite('contours.jpg', im)  
cv.waitKey()