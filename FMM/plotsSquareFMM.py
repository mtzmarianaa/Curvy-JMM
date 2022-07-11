# Script to generate plots from the simple square domain and fmm

# SCRIPT TO VISUALIZE ERRORS (can I say this?)
from dataclasses import replace
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from numpy.linalg import norm 
from math import sqrt


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('magma')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

######################################################
## 1. Plot of the output

eik_vals = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/ComputedValues.bin")
eik_coords = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/MeshPoints.txt", delimiter=",")
exact_values = []
errors = []
for i in range(len(eik_coords)):
    sol = sqrt( eik_coords[i, 0]**2 + eik_coords[i, 1]**2 )
    exact_values += [sol]
    errors += [ abs( sol - eik_vals[i] ) ]

colors = colormap1(eik_vals)

plt.figure(1)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals, cmap = colormap2)
plt.title("Computed eikonal value, square")
plt.show(block=False)

plt.figure(2)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords[:, 0], eik_coords[:, 1], eik_vals, c= eik_vals, cmap=colormap2)
plt.title("Computed eikonal values, square")
plt.show(block = False)

plt.figure(3)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = exact_values, cmap = colormap2)
plt.title("Exact solution, square")
plt.show(block=False)

plt.figure(4)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords[:, 0], eik_coords[:, 1], exact_values, c= exact_values, cmap=colormap2)
plt.title("Exact solution, square")
plt.show(block = False)


plt.figure(5)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = errors, cmap = colormap2)
plt.title("Computed errors per point, square")
plt.show(block = False)

plt.figure(6)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords[:, 0], eik_coords[:, 1], errors, c = errors, cmap=colormap2)
plt.title("Computed errors per point, square")
plt.show(block = False)


plt.show()
