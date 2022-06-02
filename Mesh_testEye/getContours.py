import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import cv2 as cv
import os

os.chdir(os.path.dirname(os.path.abspath(__file__))) 

im = cv.imread('ommatidium_diagram2.jpg')
cv.imshow('image',im)
imgray = cv.cvtColor(im, cv.COLOR_BGR2GRAY)
ret, thresh = cv.threshold(imgray, 127, 255, 0)
contours, hierarchy = cv.findContours(thresh, cv.RETR_CCOMP, cv.CHAIN_APPROX_NONE)

# Drawing the contours on top of the original image from the OPENCV method
# cnt_vec = np.arange(-9 ,0)
# for cnt_i in cnt_vec: 
#     cnt = contours[cnt_i]
#     cv.drawContours(im, [cnt], 0, (255, 255, 0), 3 )
# cv.imshow('img', im)
# cv.imwrite('contours.jpg', im)  
# cv.waitKey()


# Get the list of points
x_coord = []
y_coord = []
#levels_to_use = [0,1,2,3,4,5,6,7,9]
levels_to_use = [4]
for c in levels_to_use:
    n_contour = contours[c]
    for d in range(len(n_contour)):
        XY_Coordinates = n_contour[d]
        x_coord.append(XY_Coordinates[0][0])
        y_coord.append(XY_Coordinates[0][1])

y_coord = [-1*y for y in y_coord]
points_df = pd.DataFrame({'x': x_coord, 'y': y_coord})

points_df.plot(x ='x', y='y', style="o", ms=1, colormap = 'viridis', mark_right = False, legend = False)	
plt.xlim([1000, 1500])
plt.ylim([-imgray.shape[0], 0])
plt.show()

# Saving such list

points_df.to_csv('EdgesData_4.csv', index=False, header = False)


