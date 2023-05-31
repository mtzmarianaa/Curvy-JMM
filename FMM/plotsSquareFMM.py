# Script to generate plots from the simple square domain and fmm

# SCRIPT TO VISUALIZE errors_H1 (can I say this?)
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
nPointsH = []

######################################################
######################################################
######################################################
#### H1
## 1. Plot of the output

eik_vals_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_ComputedValues.bin")
eik_coords_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_MeshPoints.txt", delimiter=",")
triangles_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_Faces.txt", delimiter=",")
eik_grads_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_ComputedGradients.bin");
eik_grads_H1 = eik_grads_H1.reshape(len(eik_coords_H1), 2)

exact_values_H1 = []
errors_H1 = []
for i in range(len(eik_coords_H1)):
    sol = sqrt( eik_coords_H1[i, 0]**2 + eik_coords_H1[i, 1]**2 )
    exact_values_H1 += [sol]
    errors_H1 += [ abs( sol - eik_vals_H1[i] ) ]



plt.figure(1)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_vals_H1, c= eik_vals_H1, cmap=colormap2)
plt.title("Computed eikonal values, square h1")
plt.show(block = False)


plt.figure(2)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], exact_values_H1, c= exact_values_H1, cmap=colormap2)
plt.title("Exact solution, square h1")
plt.show(block = False)


plt.figure(3)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], errors_H1, c = errors_H1, cmap=colormap2)
plt.title("Computed errors_H1 per point, square h1")
plt.show(block = False)




# We interpolate the solution on the triangles_H1 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H1)
zi_lin = interp_lin(xi, yi)

# Contours of the errors_H1 in 3D and 2D
solution_interpolated = exact_solution(xi, yi)
errors_H1_abs = abs(zi_lin - solution_interpolated)

plt.figure(4)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_H1_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('errors_H1');
plt.title("3D point wise errors, square h1")
plt.show(block = False)


# Plot the absolute errors_H1 in 2D

plt.figure(5)
plt.axis('equal')
im5 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors, square h1")
plt.show(block = False)
plt.colorbar(im5)


# The absolute errors_H1 in 2D with the triangulation

plt.figure(6)
plt.axis('equal')
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#04007e')
im6 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors and triangulation, square h1")
plt.show(block = False)
plt.colorbar(im6)


#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(7)
plt.axis('equal')
im7 = plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], c = eik_vals_H1, cmap = colormap2)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square h1")
plt.show(block = False)
plt.colorbar(im7)


plt.figure(8)
plt.axis('equal')
im8 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Linear interpolation, square h1")
plt.show(block = False)
plt.colorbar(im8)


plt.figure(9)
plt.axis('equal')
im9 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.quiver(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_grads_H1[:, 0], eik_grads_H1[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, square h1")
plt.show(block = False)
plt.colorbar(im9)


averageH += [average_edge_length(eik_coords_H1, triangles_H1)]
errorNorm += [norm( np.subtract(eik_vals_H1, exact_values_H1)  )/norm( exact_values_H1 )]
nPointsH += [len(eik_coords_H1)]

######################################################
######################################################
######################################################
#### H2
## 1. Plot of the output

eik_vals_H2 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H2_ComputedValues.bin")
eik_coords_H2 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H2_MeshPoints.txt", delimiter=",")
triangles_H2 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H2_Faces.txt", delimiter=",")
eik_grads_H2 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H2_ComputedGradients.bin");
eik_grads_H2 = eik_grads_H2.reshape(len(eik_coords_H2), 2)

exact_values_H2 = []
errors_H2 = []
for i in range(len(eik_coords_H2)):
    sol = sqrt( eik_coords_H2[i, 0]**2 + eik_coords_H2[i, 1]**2 )
    exact_values_H2 += [sol]
    errors_H2 += [ abs( sol - eik_vals_H2[i] ) ]



plt.figure(10)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], eik_vals_H2, c= eik_vals_H2, cmap=colormap2)
plt.title("Computed eikonal values, square h2")
plt.show(block = False)


plt.figure(11)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], exact_values_H2, c= exact_values_H2, cmap=colormap2)
plt.title("Exact solution, square h2")
plt.show(block = False)


plt.figure(12)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], errors_H2, c = errors_H2, cmap=colormap2)
plt.title("Computed errors per point, square h2")
plt.show(block = False)




# We interpolate the solution on the triangles_H2 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H2)
zi_lin = interp_lin(xi, yi)

# Contours of the errors_H2 in 3D and 2D
solution_interpolated = exact_solution(xi, yi)
errors_H2_abs = abs(zi_lin - solution_interpolated)

plt.figure(13)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_H2_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('errors_H2');
plt.title("3D point wise errors, square h2")
plt.show(block = False)


# Plot the absolute errors_H2 in 2D

plt.figure(14)
plt.axis('equal')
im14 = plt.imshow( errors_H2_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors, square h2")
plt.show(block = False)
plt.colorbar(im14)


# The absolute errors_H2 in 2D with the triangulation

plt.figure(15)
plt.axis('equal')
plt.triplot(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2, '-.', lw=0.5, c='#04007e')
im15 = plt.imshow( errors_H2_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors and triangulation, square h2")
plt.show(block = False)
plt.colorbar(im15)


#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(16)
plt.axis('equal')
im16 = plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], c = eik_vals_H2, cmap = colormap2)
plt.triplot(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square h2")
plt.show(block = False)
plt.colorbar(im16)


plt.figure(17)
plt.axis('equal')
im17 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Linear interpolation, square h2")
plt.show(block = False)
plt.colorbar(im17)


plt.figure(18)
plt.axis('equal')
im18 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.quiver(eik_coords_H2[:, 0], eik_coords_H2[:, 1], eik_grads_H2[:, 0], eik_grads_H2[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, square h2")
plt.show(block = False)
plt.colorbar(im18)


averageH += [average_edge_length(eik_coords_H2, triangles_H2)]
errorNorm += [norm( np.subtract(eik_vals_H2, exact_values_H2)  )/norm( exact_values_H2 )]
nPointsH += [len(eik_coords_H2)]

######################################################
######################################################
######################################################
#### H3
## 1. Plot of the output

eik_vals_H3 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H3_ComputedValues.bin")
eik_coords_H3 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H3_MeshPoints.txt", delimiter=",")
triangles_H3 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H3_Faces.txt", delimiter=",")
eik_grads_H3 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H3_ComputedGradients.bin");
eik_grads_H3 = eik_grads_H3.reshape(len(eik_coords_H3), 2)

exact_values_H3 = []
errors_H3 = []
for i in range(len(eik_coords_H3)):
    sol = sqrt( eik_coords_H3[i, 0]**2 + eik_coords_H3[i, 1]**2 )
    exact_values_H3 += [sol]
    errors_H3 += [ abs( sol - eik_vals_H3[i] ) ]



plt.figure(19)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], eik_vals_H3, c= eik_vals_H3, cmap=colormap2)
plt.title("Computed eikonal values, square h3")
plt.show(block = False)


plt.figure(20)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], exact_values_H3, c= exact_values_H3, cmap=colormap2)
plt.title("Exact solution, square h3")
plt.show(block = False)


plt.figure(21)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], errors_H3, c = errors_H3, cmap=colormap2)
plt.title("Computed errors per point, square h3")
plt.show(block = False)




# We interpolate the solution on the triangles_H3 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H3)
zi_lin = interp_lin(xi, yi)

# Contours of the errors_H3 in 3D and 2D
solution_interpolated = exact_solution(xi, yi)
errors_H3_abs = abs(zi_lin - solution_interpolated)

plt.figure(22)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_H3_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('errors_H3');
plt.title("3D point wise errors, square h3")
plt.show(block = False)


# Plot the absolute errors_H3 in 2D

plt.figure(23)
plt.axis('equal')
im23 = plt.imshow( errors_H3_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors, square h3")
plt.show(block = False)
plt.colorbar(im23)


# The absolute errors_H3 in 2D with the triangulation

plt.figure(24)
plt.axis('equal')
plt.triplot(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3, '-.', lw=0.5, c='#04007e')
im24 = plt.imshow( errors_H3_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors and triangulation, square h3")
plt.show(block = False)
plt.colorbar(im24)


#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(25)
plt.axis('equal')
im25 = plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], c = eik_vals_H3, cmap = colormap2)
plt.triplot(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square h3")
plt.show(block = False)
plt.colorbar(im25)


plt.figure(26)
plt.axis('equal')
im26 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Linear interpolation, square h3")
plt.show(block = False)
plt.colorbar(im26)


plt.figure(27)
plt.axis('equal')
im27 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.quiver(eik_coords_H3[:, 0], eik_coords_H3[:, 1], eik_grads_H3[:, 0], eik_grads_H3[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, square h3")
plt.show(block = False)
plt.colorbar(im27)


averageH += [average_edge_length(eik_coords_H3, triangles_H3)]
errorNorm += [norm( np.subtract(eik_vals_H3, exact_values_H3)  )/norm( exact_values_H3 )]
nPointsH += [len(eik_coords_H3)]


######################################################
######################################################
######################################################
#### H4
## 1. Plot of the output

eik_vals_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H4_ComputedValues.bin")
eik_coords_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H4_MeshPoints.txt", delimiter=",")
triangles_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H4_Faces.txt", delimiter=",")
eik_grads_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H4_ComputedGradients.bin");
eik_grads_H4 = eik_grads_H4.reshape(len(eik_coords_H4), 2)

exact_values_H4 = []
errors_H4 = []
for i in range(len(eik_coords_H4)):
    sol = sqrt( eik_coords_H4[i, 0]**2 + eik_coords_H4[i, 1]**2 )
    exact_values_H4 += [sol]
    errors_H4 += [ abs( sol - eik_vals_H4[i] ) ]



plt.figure(28)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], eik_vals_H4, c= eik_vals_H4, cmap=colormap2)
plt.title("Computed eikonal values, square h4")
plt.show(block = False)


plt.figure(29)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], exact_values_H4, c= exact_values_H4, cmap=colormap2)
plt.title("Exact solution, square h4")
plt.show(block = False)


plt.figure(30)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], errors_H4, c = errors_H4, cmap=colormap2)
plt.title("Computed errors per point, square h4")
plt.show(block = False)




# We interpolate the solution on the triangles_H4 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H4)
zi_lin = interp_lin(xi, yi)

# Contours of the errors_H4 in 3D and 2D
solution_interpolated = exact_solution(xi, yi)
errors_H4_abs = abs(zi_lin - solution_interpolated)

plt.figure(31)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_H4_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('errors_H4');
plt.title("3D point wise errors, square h4")
plt.show(block = False)


# Plot the absolute errors_H4 in 2D

plt.figure(32)
plt.axis('equal')
im14 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors, square h4")
plt.show(block = False)
plt.colorbar(im14)


# The absolute errors_H4 in 2D with the triangulation

plt.figure(33)
plt.axis('equal')
plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#04007e')
im33 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors and triangulation, square h4")
plt.show(block = False)
plt.colorbar(im33)


#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(34)
plt.axis('equal')
im34 = plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], c = eik_vals_H4, cmap = colormap2)
plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square h4")
plt.show(block = False)
plt.colorbar(im34)


plt.figure(35)
plt.axis('equal')
im35 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Linear interpolation, square h4")
plt.show(block = False)
plt.colorbar(im35)


plt.figure(36)
plt.axis('equal')
im36 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.quiver(eik_coords_H4[:, 0], eik_coords_H4[:, 1], eik_grads_H4[:, 0], eik_grads_H4[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, square h4")
plt.show(block = False)
plt.colorbar(im36)


averageH += [average_edge_length(eik_coords_H4, triangles_H4)]
errorNorm += [norm( np.subtract(eik_vals_H4, exact_values_H4)  )/norm( exact_values_H4 )]
nPointsH += [len(eik_coords_H4)]

######################################################
######################################################
######################################################
#### H5
## 1. Plot of the output

eik_vals_H5 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H5_ComputedValues.bin")
eik_coords_H5 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H5_MeshPoints.txt", delimiter=",")
triangles_H5 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H5_Faces.txt", delimiter=",")
eik_grads_H5 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H5_ComputedGradients.bin");
eik_grads_H5 = eik_grads_H5.reshape(len(eik_coords_H5), 2)

exact_values_H5 = []
errors_H5 = []
for i in range(len(eik_coords_H5)):
    sol = sqrt( eik_coords_H5[i, 0]**2 + eik_coords_H5[i, 1]**2 )
    exact_values_H5 += [sol]
    errors_H5 += [ abs( sol - eik_vals_H5[i] ) ]



plt.figure(37)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], eik_vals_H5, c= eik_vals_H5, cmap=colormap2)
plt.title("Computed eikonal values, square h5")
plt.show(block = False)


plt.figure(38)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], exact_values_H5, c= exact_values_H5, cmap=colormap2)
plt.title("Exact solution, square h5")
plt.show(block = False)


plt.figure(39)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], errors_H5, c = errors_H5, cmap=colormap2)
plt.title("Computed errors per point, square h5")
plt.show(block = False)




# We interpolate the solution on the triangles_H5 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H5)
zi_lin = interp_lin(xi, yi)

# Contours of the errors_H5 in 3D and 2D
solution_interpolated = exact_solution(xi, yi)
errors_H5_abs = abs(zi_lin - solution_interpolated)

plt.figure(40)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_H5_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('errors_H5');
plt.title("3D point wise errors, square h5")
plt.show(block = False)


# Plot the absolute errors_H5 in 2D

plt.figure(41)
plt.axis('equal')
im41 = plt.imshow( errors_H5_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors, square h5")
plt.show(block = False)
plt.colorbar(im41)


# The absolute errors_H5 in 2D with the triangulation

plt.figure(42)
plt.axis('equal')
plt.triplot(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5, '-.', lw=0.5, c='#04007e')
im42 = plt.imshow( errors_H5_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors and triangulation, square h5")
plt.show(block = False)
plt.colorbar(im42)


#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(43)
plt.axis('equal')
im43 = plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], c = eik_vals_H5, cmap = colormap2)
plt.triplot(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square h5")
plt.show(block = False)
plt.colorbar(im43)


plt.figure(44)
plt.axis('equal')
im44 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Linear interpolation, square h5")
plt.show(block = False)
plt.colorbar(im44)


plt.figure(45)
plt.axis('equal')
im45 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.quiver(eik_coords_H5[:, 0], eik_coords_H5[:, 1], eik_grads_H5[:, 0], eik_grads_H5[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, square h5")
plt.show(block = False)
plt.colorbar(im45)


averageH += [average_edge_length(eik_coords_H5, triangles_H5)]
errorNorm += [norm( np.subtract(eik_vals_H5, exact_values_H5)  )/norm( exact_values_H5 )]
nPointsH += [len(eik_coords_H5)]

######################################################
######################################################
######################################################
#### H6
## 1. Plot of the output

eik_vals_H6 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H6_ComputedValues.bin")
eik_coords_H6 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H6_MeshPoints.txt", delimiter=",")
triangles_H6 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H6_Faces.txt", delimiter=",")
eik_grads_H6 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H6_ComputedGradients.bin");
eik_grads_H6 = eik_grads_H6.reshape(len(eik_coords_H6), 2)

exact_values_H6 = []
errors_H6 = []
for i in range(len(eik_coords_H6)):
    sol = sqrt( eik_coords_H6[i, 0]**2 + eik_coords_H6[i, 1]**2 )
    exact_values_H6 += [sol]
    errors_H6 += [ abs( sol - eik_vals_H6[i] ) ]



plt.figure(46)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], eik_vals_H6, c= eik_vals_H6, cmap=colormap2)
plt.title("Computed eikonal values, square h6")
plt.show(block = False)


plt.figure(47)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], exact_values_H6, c= exact_values_H6, cmap=colormap2)
plt.title("Exact solution, square h6")
plt.show(block = False)


plt.figure(48)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], errors_H6, c = errors_H6, cmap=colormap2)
plt.title("Computed errors per point, square h6")
plt.show(block = False)




# We interpolate the solution on the triangles_H6 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-5, 5, 1000), np.linspace(-5, 5, 1000))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H6)
zi_lin = interp_lin(xi, yi)

# Contours of the errors_H6 in 3D and 2D
solution_interpolated = exact_solution(xi, yi)
errors_H6_abs = abs(zi_lin - solution_interpolated)

plt.figure(49)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_H6_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('errors_H6');
plt.title("3D point wise errors, square h6")
plt.show(block = False)


# Plot the absolute errors_H6 in 2D

plt.figure(50)
plt.axis('equal')
im50 = plt.imshow( errors_H6_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors, square h6")
plt.show(block = False)
plt.colorbar(im50)


# The absolute errors_H6 in 2D with the triangulation

plt.figure(51)
plt.axis('equal')
plt.triplot(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6, '-.', lw=0.5, c='#04007e')
im51 = plt.imshow( errors_H6_abs, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Point wise absolute errors and triangulation, square h6")
plt.show(block = False)
plt.colorbar(im51)


#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
plt.figure(52)
plt.axis('equal')
im52 = plt.contourf(xi, yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], c = eik_vals_H6, cmap = colormap2)
plt.triplot(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, square h6")
plt.show(block = False)
plt.colorbar(im52)


plt.figure(53)
plt.axis('equal')
im53 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.title("Linear interpolation, square h6")
plt.show(block = False)
plt.colorbar(im53)


plt.figure(54)
plt.axis('equal')
im54 = plt.imshow( zi_lin, cmap = colormap2, extent=[-5,5,-5,5]  )
plt.quiver(eik_coords_H6[:, 0], eik_coords_H6[:, 1], eik_grads_H6[:, 0], eik_grads_H6[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, square h6")
plt.show(block = False)
plt.colorbar(im54)


averageH += [average_edge_length(eik_coords_H6, triangles_H6)]
errorNorm += [norm( np.subtract(eik_vals_H6, exact_values_H6)  )/norm( exact_values_H6 )]
nPointsH += [len(eik_coords_H6)]

######################################################
######################################################
######################################################
######################################################
################## ERRORS ############################
################### EACH #############################
####################  H  #############################
######################################################
######################################################


plt.figure(55)
plt.loglog(averageH, errorNorm, c = '#6800ff')
plt.title("l2 errors and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")

plt.figure(56)
plt.loglog(averageH, nPointsH, c = '#6800ff')
plt.title("l2 errors and number of points in triangulation")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")


plt.show()
