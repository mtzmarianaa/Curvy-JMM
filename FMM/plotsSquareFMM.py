# Script to generate plots from the simple square domain and fmm

# SCRIPT TO VISUALIZE ERRORS (can I say this?)
from dataclasses import replace
from filecmp import cmp
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


def average_edge_length(eik_coords, faces):
    #for each edge we calculate its length and then we calculate the average edge length for a triangulation
    sum = 0
    nEdges = 0
    n_points = len(eik_coords)
    counted = np.zeros((n_points, n_points))
    for i in range(len(faces)):
        p1 = int(faces[i, 0])
        p2 = int(faces[i, 1])
        p3 = int(faces[i, 2])
        if counted[p1, p2] == 0:
            counted[p1, p2] = 1
            nEdges += 1
            sum += sqrt(  (eik_coords[p1, 0] - eik_coords[p2, 0])**2 +  (eik_coords[p1, 1] - eik_coords[p2, 1])**2 )
        if counted[p1, p3] == 0:
            counted[p1, p3] = 1
            nEdges += 1
            sum += sqrt(  (eik_coords[p1, 0] - eik_coords[p3, 0])**2 +  (eik_coords[p1, 1] - eik_coords[p3, 1])**2 )
        if counted[p2, p3] == 0:
            counted[p2, p3] = 1
            nEdges += 1
            sum += sqrt(  (eik_coords[p2, 0] - eik_coords[p3, 0])**2 +  (eik_coords[p2, 1] - eik_coords[p3, 1])**2 )
    return sum/nEdges
        
        
def exact_solution(xi, yi):
    lenxi = len(xi)
    lenyi = len(yi)
    zi = np.zeros((lenxi, lenyi))
    for i in range(lenxi):
        for j in range(lenyi):
            zi[i,j] = sqrt( xi[i, j]**2 + yi[i, j]**2   )
    return zi   
        
n = 0
averageH = []
errorNorm = []

######################################################
######################################################
######################################################
#### H1
## 1. Plot of the output

eik_vals_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H4_ComputedValues.bin")
eik_coords_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H4_MeshPoints.txt", delimiter=",")
triangles_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H4_Faces.txt", delimiter=",")


exact_values = []
errors = []
for i in range(len(eik_coords_H1)):
    sol = sqrt( eik_coords_H1[i, 0]**2 + eik_coords_H1[i, 1]**2 )
    exact_values += [sol]
    errors += [ abs( sol - eik_vals_H1[i] ) ]

colors = colormap1(eik_vals_H1)


plt.figure(1)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_vals_H1, c= eik_vals_H1, cmap=colormap2)
plt.title("Computed eikonal values, square")
plt.show(block = False)


plt.figure(2)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], exact_values, c= exact_values, cmap=colormap2)
plt.title("Exact solution, square")
plt.show(block = False)


plt.figure(3)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], errors, c = errors, cmap=colormap2)
plt.title("Computed errors per point, square")
plt.show(block = False)




# We interpolate the solution on the triangles_H1 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H1)
zi_lin = interp_lin(xi, yi)

# Contours of the errors in 3D and 2D
solution_interpolated = exact_solution(xi, yi)
errors_abs = abs(zi_lin - solution_interpolated)

plt.figure(6)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Errors');
plt.title("3D point wise errors")
plt.show(block = False)


# Plot the absolute errors in 2D

plt.figure(7)
plt.axis('equal')
im7 = plt.imshow( errors_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors, square")
plt.show(block = False)
plt.colorbar(im7)


# The absolute errors in 2D with the triangulation

plt.figure(8)
plt.axis('equal')
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#04007e')
im8 = plt.imshow( errors_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors and triangulation, square")
plt.show(block = False)
plt.colorbar(im8)


#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(9)
plt.axis('equal')
im9 = plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], c = eik_vals_H1, cmap = colormap2)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square")
plt.show(block = False)
plt.colorbar(im9)


plt.figure(10)
plt.axis('equal')
im10 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Linear interpolation, square")
plt.show(block = False)
plt.colorbar(im10)


averageH += [average_edge_length(eik_coords_H1, triangles_H1)]
errorNorm += [norm( np.subtract(eik_vals_H1, exact_values)  )/norm( exact_values )]

n += 1


######################################################
######################################################
######################################################
#### H2
## 1. Plot of the output



# plt.figure(n+1)
# plt.plot(averageH, errorNorm, marker = '-.', c = '#6800ff')
# plt.title("Convergence")
# plt.xlabel("Average edge length")
# plt.ylabel("Error")


plt.show()
