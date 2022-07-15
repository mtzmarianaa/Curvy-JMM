# Script to generate plots from the square with inverted triangle and fmm

# SCRIPT TO VISUALIZE ERRORS (can I say this?)
from dataclasses import replace
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from numpy.linalg import norm 
from math import sqrt
import matplotlib.tri as tri


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('magma')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

def averageSolutionFace(eikvals, faces):
    avSol = []
    for i in range(len(faces)):
        avSol += [( eikvals[int(faces[i, 0])]+ eikvals[int(faces[i, 1])]+eikvals[int(faces[i, 2])]   )/3]
    return avSol

######################################################
## 1. Plot of the output

eik_vals = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/ComputedValues.bin")
eik_coords = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/MeshPoints.txt", delimiter=",")
triangles = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/Faces.txt", delimiter=",")

averageSolution = averageSolutionFace(eik_vals, triangles)

# exact_values = []
# errors = []
# for i in range(len(eik_coords)):
#     sol = sqrt( eik_coords[i, 0]**2 + eik_coords[i, 1]**2 )
#     exact_values += [sol]
#     errors += [ abs( sol - eik_vals[i] ) ]

colors = colormap1(eik_vals)

plt.figure(1)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals, cmap = colormap2)
plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles, '-.', lw=0.5, c='#6800ff')
plt.title("Computed eikonal value, square with inverted triangle")
plt.show(block=False)

plt.figure(2)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords[:, 0], eik_coords[:, 1], eik_vals, c= eik_vals, cmap=colormap2)
plt.title("Computed eikonal values, square with inverted triangle")
plt.show(block = False)

plt.figure(3)
plt.tripcolor(eik_coords[:, 0], eik_coords[:, 1], triangles, averageSolution, cmap = colormap2)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals, cmap = colormap2, marker=".")
#plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles, '-.', lw=0.5, c='#6800ff')
plt.title("Computed eikonal value, square with inverted triangle")
plt.show(block=False)


# Now we use the linear interpolation thing to plot the contours

xi, yi = np.meshgrid(np.linspace(-10, 10, 100), np.linspace(-10, 10, 100))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords[:, 0], eik_coords[:, 1], triangles)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
zi_lin = interp_lin(xi, yi)
#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(4)
plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.plot(xi, yi, 'k-', lw=0.5, alpha=0.3)
plt.plot(xi.T, yi.T, 'k', lw=0.5, alpha=0.3)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals, cmap = colormap2, marker='.')
plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square with inverted triangle")


#####################################################################################
#####################################################################################
#####################################################################################
# Same but with different starting point

eik_vals2 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/ComputedValues2.bin")

averageSolution2 = averageSolutionFace(eik_vals2, triangles)

colors = colormap1(eik_vals2)

plt.figure(5)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals2, cmap = colormap2)
plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles, '-.', lw=0.5, c='#6800ff')
plt.title("Computed eikonal value, square with inverted triangle")
plt.show(block=False)

plt.figure(6)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords[:, 0], eik_coords[:, 1], eik_vals2, c= eik_vals, cmap=colormap2)
plt.title("Computed eikonal values, square with inverted triangle")
plt.show(block = False)

plt.figure(7)
plt.tripcolor(eik_coords[:, 0], eik_coords[:, 1], triangles, averageSolution2, cmap = colormap2)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals2, cmap = colormap2, marker=".")
#plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles, '-.', lw=0.5, c='#6800ff')
plt.title("Computed eikonal value, square with inverted triangle")
plt.show(block=False)

# Now we use the linear interpolation thing to plot the contours

xi, yi = np.meshgrid(np.linspace(-10, 10, 100), np.linspace(-10, 10, 100))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords[:, 0], eik_coords[:, 1], triangles)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals2)
zi_lin = interp_lin(xi, yi)
#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(8)
plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.plot(xi, yi, 'k-', lw=0.5, alpha=0.3)
plt.plot(xi.T, yi.T, 'k', lw=0.5, alpha=0.3)
plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals2, cmap = colormap2, marker='.')
plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square with inverted triangle")


plt.show()