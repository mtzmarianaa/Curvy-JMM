# Script to generate plots from the square with inverted triangle and fmm

# SCRIPT TO VISUALIZE ERRORS (can I say this?)
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm 
from math import sqrt
import matplotlib.tri as tri
import matplotlib.animation as animation
import tabulate
import AnalyticSolutionTestGeometry as AnSol


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('magma')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

def rotate(angle):
    ax.view_init(azim=angle)

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

def constrainOnBoundary(x):
    '''
    function for the constrain |x[0]| = x[1] (i.e. the specification that xlambda is ON the boundary between the two regions)
    '''
    return x[1] - abs(x[0])


def IndexRefractionRegions(x):
    '''
    Slowness function according to the two sections present in this particular domain
    '''
    if constrainOnBoundary(x)> 0:
        s = 2
    else:
        s = 1
    return s
        
       
        
n = 0
averageH = []
errorNorm = []
nPointsH = []

times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Times.bin")

######################################################
######################################################
######################################################
#### H1
## 1. Plot of the output

eik_vals_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1_ComputedValues.bin")
eik_coords_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1_MeshPoints.txt", delimiter=",")
triangles_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1_Faces.txt", delimiter=",")
eik_grads_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1_ComputedGradients.bin");
eik_grads_H1 = eik_grads_H1.reshape(len(eik_coords_H1), 2)

exact_values_H1 = []
errors_H1 = []
for i in range(len(eik_coords_H1)):
    xi = eik_coords_H1[i, 0]
    yi = eik_coords_H1[i, 1]
    sol = exact_solution1(xi, yi)
    exact_values_H1 += [sol]
    errors_H1 += [ abs( sol - eik_vals_H1[i] ) ]



my_dpi=96
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_vals_H1, c= eik_vals_H1, cmap=colormap2)
plt.title("Computed eikonal values, test geometry h1")
plt.show(block = False)
rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_ComputedValues.gif', dpi=80, writer='imagemagick')


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], exact_values_H1, c= exact_values_H1, cmap=colormap2)
plt.title("Exact solution, test geometry h1")
plt.show(block = False)
rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_ExactSolution.gif', dpi=80, writer='imagemagick')


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], errors_H1, c = errors_H1, cmap=colormap2)
plt.title("Computed errors per point, test geometry h1")
plt.show(block = False)
rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_PointsPointErrors.gif', dpi=80, writer='imagemagick')



# We interpolate the solution on the triangles_H1 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-10, 10, 500), np.linspace(-10, 10, 500))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H1)
zi_lin = interp_lin(xi, -yi)

# Contours of the errors_H1 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(len(xi)):
    for j in range(len(yi)):
        solution_interpolated[i, j] = exact_solution1(  xi[i, j], -yi[i,j]  )
errors_H1_abs = abs(zi_lin - solution_interpolated)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im4 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-10,10,-10,10]  )
plt.title("Exact solution, test geometry h1")
plt.show(block = False)
plt.colorbar(im4)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_ExactSolution.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
ax = plt.axes(projection='3d')
ax.contour3D(xi, yi, errors_H1_abs , 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('errors_H1');
plt.title("3D point wise errors, test geometry h1")
plt.show(block = False)
rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_PointErrors.gif', dpi=80, writer='imagemagick')


# Plot the absolute errors_H1 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im6 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-10,10,-10,10]  )
plt.title("Point wise absolute errors, test geometry h1")
plt.show(block = False)
plt.colorbar(im6)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H1 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#04007e')
im7 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-10,10,-10,10]  )
plt.title("Point wise absolute errors and triangulation, test geometry h1")
plt.show(block = False)
plt.colorbar(im7)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im8 = plt.contourf(xi, -yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], c = eik_vals_H1, cmap = colormap2)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry h1")
plt.show(block = False)
plt.colorbar(im8)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im9 = plt.imshow( zi_lin, cmap = colormap2, extent=[-10,10,-10,10]  )
plt.title("Linear interpolation, test geometry h1")
plt.show(block = False)
plt.colorbar(im9)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im10 = plt.imshow( zi_lin, cmap = colormap2, extent=[-10,10,-10,10]  )
plt.quiver(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_grads_H1[:, 0], eik_grads_H1[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry h1")
plt.show(block = False)
plt.colorbar(im10)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H1/H1_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H1, triangles_H1)]
errorNorm += [norm( np.subtract(eik_vals_H1, exact_values_H1)  )/norm( exact_values_H1 )]
nPointsH += [len(eik_coords_H1)]