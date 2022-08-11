# Script to generate plots from the square with just a circle

# SCRIPT TO VISUALIZE ERRORS 
from cmath import asin
from string import whitespace
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm 
from math import sqrt, cos, sin, pi, atan
from scipy.optimize import minimize, minimize_scalar
import matplotlib.animation as animation
import tabulate
from matplotlib.patches import Arc
from matplotlib.colors import ListedColormap
from matplotlib.ticker import FuncFormatter
from analyticSol_circle import trueSolution
import matplotlib.tri as tri


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)


nx = 36*10
ny = 42*10
my_dpi=96
eta1 = 1.0
eta2 = 1.452
x0 = np.array([-15, -10])
center = np.array([0,0])
R = 10
eps = np.finfo(np.float64).resolution

#########################################################
####### USEFUL FUNCTIONS

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

def getPathFromIndex(eik_coords, index_to_get_path_from, parents_path, lambdas):
    '''
    Given the coordinates, list of parents and list of optimal lambdas we get the path of index_to_get_path_from
    '''
    path = [(eik_coords[index_to_get_path_from, 0], eik_coords[index_to_get_path_from, 1])] # the path ends/starts in the index of interest
    queue_indices = [index_to_get_path_from]
    queue = [ (parents_path[index_to_get_path_from, 0],  parents_path[index_to_get_path_from, 1]) ] # queue, the indices in the current level considered
    while( queue ): #while there are elements on the list
        # for the pairs of parents in the queue we calculate the mean of the paths (i.e. means of (1-lambda)parent0 + lambda*parent1)
        sum_x = 0
        sum_y = 0
        count = 0
        len_queue = len(queue) #number of pairs of parents in the queue
        for i in range(len_queue):
            lamb = lambdas[ queue_indices[i]]
            sum_x += (1 - lamb)*eik_coords[ queue[i][0] , 0] + lamb*eik_coords[ queue[i][1], 0 ] #add the coordinates of the xlambdas
            sum_y += (1 - lamb)*eik_coords[ queue[i][0] , 1] + lamb*eik_coords[ queue[i][1], 1 ]
            count += 1
        path.extend( [  (sum_x/count, sum_y/count)  ] ) #add the mean of the xlambdas to que path
        queue_indices = [ n[0] for n in queue ] + [ n[1] for n in queue ]
        if len(set(queue_indices)) == 0:
            queue = [] # if all the new indices are the same then it means that we've reached a source/starting point
            path.extend(  (eik_coords[queue_indices[0], 0], eik_coords[queue_indices[0], 1]  )   ) # we add the source to the path
        else:
            queue = [  (parents_path[ind, 0], parents_path[ind, 0] ) for ind in queue_indices   ] # we have a new queue
    # we need to reverse this list because it starts at the end point and ends in the source (just for aesthetics)
    path.reverse()
    return path
        
        
        
    


n = 0
averageH = []
errorNorm = []
nPointsH = []

times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/Times.bin")
print(len(times))
# Compute the analytic solution in a grid

xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
true_solGrid = np.zeros(xi.shape)
type_solution = np.zeros(xi.shape)



for i in range(ny):
    for j in range(nx):
        sol, typeSol = trueSolution(  xi[i, j], yi[i,j], x0, center, R, eta1, eta2  )
        true_solGrid[i, j] = sol
        type_solution[i, j] = typeSol

# We plot the true solution

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
im1 = plt.imshow( true_solGrid, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
fig.subplots_adjust(wspace=None)
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.title("Exact solution, test geometry just base")
plt.show(block = False)
plt.colorbar(im1)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ExactSolution8.png', dpi=my_dpi * 10)


# Plot the type of solution
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2 = plt.imshow( type_solution, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Type of solution solution, test geometry just base")
plt.show(block = False)
plt.colorbar(im2)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/TypeSolution8.png', dpi=my_dpi * 10)

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# # Plot the surface.
# surf = ax.plot_surface(xi, yi, true_solGrid, cmap=colormap2, linewidth=0, antialiased=False)
# plt.show(block = False)

# Plot the contours in 2D
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im_bar1 = plt.contourf(xi, yi, true_solGrid, cmap = colormap2, levels = 25)
plt.title("Exact solution, test geometry just base")
plt.show(block = False)
figName_Contour = '/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ExactSolution_Contour8.png'
plt.savefig(figName_Contour, dpi=my_dpi * 10)

# # Plot in 3D and save the gif
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(xi, yi, true_solGrid, c= true_solGrid, cmap=colormap2)
# plt.title("Exact solution, test geometry just base")
# plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ExactSolution_Contour8.gif', dpi=80, writer='Pillow')


# # Save the computed values (in case they are useful)

# np.savetxt('/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/TrueSolutionGrid.txt', true_solGrid, delimiter =', ', fmt = '%.8f' )


######################################################
######################################################
######################################################
#### H1
## 1. Plot of the output

eik_vals_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedValues.bin")
eik_coords_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_MeshPoints.txt", delimiter=",")
triangles_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Faces.txt", delimiter=",")
eik_grads_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedGradients.bin");
eik_grads_H1 = eik_grads_H1.reshape(len(eik_coords_H1), 2)

exact_values_H1 = []
errors_H1 = []
for i in range(len(eik_coords_H1)):
    xi_coords = eik_coords_H1[i, 0]
    yi_coords = eik_coords_H1[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H1 += [sol]
    errors_H1 += [ abs( sol - eik_vals_H1[i] ) ]
    


# We interpolate the solution on the triangles_H1 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H1)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

# Contours of the errors_H1 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(ny):
    for j in range(nx):
        tau, type_path = trueSolution(xi[i,j], yi[i,j], x0, center, R, eta1, eta2)
        solution_interpolated[i, j] = tau
errors_H1_abs = abs(zi_linP - solution_interpolated)



#Plot the absolute errors_H1 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im6 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base h1")
plt.show(block = False)
plt.colorbar(im6)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H1 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#04007e')
im7 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base h1")
plt.show(block = False)
plt.colorbar(im7)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im8 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], c = eik_vals_H1, cmap = colormap2)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry just base h1")
plt.show(block = False)
plt.colorbar(im8)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im9 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base h1")
plt.show(block = False)
plt.colorbar(im9)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im10 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_grads_H1[:, 0], eik_grads_H1[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base h1")
plt.show(block = False)
plt.colorbar(im10)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H1, triangles_H1)]
errorNorm += [norm( errors_H1  )/norm( exact_values_H1 )]
nPointsH += [len(eik_coords_H1)]

######################################################
######################################################
######################################################
#### H2
## 1. Plot of the output

eik_vals_H2 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedValues.bin")
eik_coords_H2 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_MeshPoints.txt", delimiter=",")
triangles_H2 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Faces.txt", delimiter=",")
eik_grads_H2 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedGradients.bin");
eik_grads_H2 = eik_grads_H2.reshape(len(eik_coords_H2), 2)

exact_values_H2 = []
errors_H2 = []
for i in range(len(eik_coords_H2)):
    xi_coords = eik_coords_H2[i, 0]
    yi_coords = eik_coords_H2[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H2 += [sol]
    errors_H2 += [ abs( sol - eik_vals_H2[i] ) ]
    


# We interpolate the solution on the triangles_H2 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H2)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

# Contours of the errors_H2 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(ny):
    for j in range(nx):
        tau, type_path = trueSolution(xi[i,j], yi[i,j], x0, center, R, eta1, eta2)
        solution_interpolated[i, j] = tau
errors_H2_abs = abs(zi_linP - solution_interpolated)



#Plot the absolute errors_H2 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im26 = plt.imshow( errors_H2_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H2")
plt.show(block = False)
plt.colorbar(im26)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H2 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2, '-.', lw=0.5, c='#04007e')
im27 = plt.imshow( errors_H2_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H2")
plt.show(block = False)
plt.colorbar(im27)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im28 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], c = eik_vals_H2, cmap = colormap2)
plt.triplot(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H2")
plt.show(block = False)
plt.colorbar(im28)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im29 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H2")
plt.show(block = False)
plt.colorbar(im29)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im210 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H2[:, 0], eik_coords_H2[:, 1], eik_grads_H2[:, 0], eik_grads_H2[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H2")
plt.show(block = False)
plt.colorbar(im210)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H2, triangles_H2)]
errorNorm += [norm( errors_H2  )/norm( exact_values_H2 )]
nPointsH += [len(eik_coords_H2)]


######################################################
######################################################
######################################################
#### H3
## 1. Plot of the output

eik_vals_H3 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_ComputedValues.bin")
eik_coords_H3 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_MeshPoints.txt", delimiter=",")
triangles_H3 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_Faces.txt", delimiter=",")
eik_grads_H3 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_ComputedGradients.bin");
eik_grads_H3 = eik_grads_H3.reshape(len(eik_coords_H3), 2)

exact_values_H3 = []
errors_H3 = []
for i in range(len(eik_coords_H3)):
    xi_coords = eik_coords_H3[i, 0]
    yi_coords = eik_coords_H3[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H3 += [sol]
    errors_H3 += [ abs( sol - eik_vals_H3[i] ) ]
    


# We interpolate the solution on the triangles_H3 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H3)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

# Contours of the errors_H3 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(ny):
    for j in range(nx):
        tau, type_path = trueSolution(xi[i,j], yi[i,j], x0, center, R, eta1, eta2)
        solution_interpolated[i, j] = tau
errors_H3_abs = abs(zi_linP - solution_interpolated)



#Plot the absolute errors_H3 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im36 = plt.imshow( errors_H3_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H3")
plt.show(block = False)
plt.colorbar(im36)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H3 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3, '-.', lw=0.5, c='#04007e')
im37 = plt.imshow( errors_H3_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H3")
plt.show(block = False)
plt.colorbar(im37)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im38 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], c = eik_vals_H3, cmap = colormap2)
plt.triplot(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H3")
plt.show(block = False)
plt.colorbar(im38)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im39 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H3")
plt.show(block = False)
plt.colorbar(im39)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im310 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H3[:, 0], eik_coords_H3[:, 1], eik_grads_H3[:, 0], eik_grads_H3[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H3")
plt.show(block = False)
plt.colorbar(im310)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H3, triangles_H3)]
errorNorm += [norm( errors_H3  )/norm( exact_values_H3 )]
nPointsH += [len(eik_coords_H3)]



######################################################
######################################################
######################################################
#### H4
## 1. Plot of the output

eik_vals_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_ComputedValues.bin")
eik_coords_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_MeshPoints.txt", delimiter=",")
triangles_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_Faces.txt", delimiter=",")
eik_grads_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_ComputedGradients.bin");
eik_grads_H4 = eik_grads_H4.reshape(len(eik_coords_H4), 2)

exact_values_H4 = []
errors_H4 = []
for i in range(len(eik_coords_H4)):
    xi_coords = eik_coords_H4[i, 0]
    yi_coords = eik_coords_H4[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H4 += [sol]
    errors_H4 += [ abs( sol - eik_vals_H4[i] ) ]
    


# We interpolate the solution on the triangles_H4 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H4)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

# Contours of the errors_H4 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(ny):
    for j in range(nx):
        tau, type_path = trueSolution(xi[i,j], yi[i,j], x0, center, R, eta1, eta2)
        solution_interpolated[i, j] = tau
errors_H4_abs = abs(zi_linP - solution_interpolated)



#Plot the absolute errors_H4 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im46 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H4")
plt.show(block = False)
plt.colorbar(im46)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H4 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#04007e')
im47 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H4")
plt.show(block = False)
plt.colorbar(im47)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im48 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], c = eik_vals_H4, cmap = colormap2)
plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H4")
plt.show(block = False)
plt.colorbar(im48)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im49 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H4")
plt.show(block = False)
plt.colorbar(im49)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im410 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H4[:, 0], eik_coords_H4[:, 1], eik_grads_H4[:, 0], eik_grads_H4[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H4")
plt.show(block = False)
plt.colorbar(im410)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H4, triangles_H4)]
errorNorm += [norm( errors_H4  )/norm( exact_values_H4 )]
nPointsH += [len(eik_coords_H4)]


######################################################
######################################################
######################################################
#### H5
## 1. Plot of the output

eik_vals_H5 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_ComputedValues.bin")
eik_coords_H5 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_MeshPoints.txt", delimiter=",")
triangles_H5 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_Faces.txt", delimiter=",")
eik_grads_H5 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_ComputedGradients.bin");
eik_grads_H5 = eik_grads_H5.reshape(len(eik_coords_H5), 2)

exact_values_H5 = []
errors_H5 = []
for i in range(len(eik_coords_H5)):
    xi_coords = eik_coords_H5[i, 0]
    yi_coords = eik_coords_H5[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H5 += [sol]
    errors_H5 += [ abs( sol - eik_vals_H5[i] ) ]
    


# We interpolate the solution on the triangles_H5 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H5)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

# Contours of the errors_H5 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(ny):
    for j in range(nx):
        tau, type_path = trueSolution(xi[i,j], yi[i,j], x0, center, R, eta1, eta2)
        solution_interpolated[i, j] = tau
errors_H5_abs = abs(zi_linP - solution_interpolated)



#Plot the absolute errors_H5 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im56 = plt.imshow( errors_H5_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H5")
plt.show(block = False)
plt.colorbar(im56)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H5 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5, '-.', lw=0.5, c='#04007e')
im57 = plt.imshow( errors_H5_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H5")
plt.show(block = False)
plt.colorbar(im57)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im58 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], c = eik_vals_H5, cmap = colormap2)
plt.triplot(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H5")
plt.show(block = False)
plt.colorbar(im58)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im59 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H5")
plt.show(block = False)
plt.colorbar(im59)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im510 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H5[:, 0], eik_coords_H5[:, 1], eik_grads_H5[:, 0], eik_grads_H5[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H5")
plt.show(block = False)
plt.colorbar(im510)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H5, triangles_H5)]
errorNorm += [norm( errors_H5  )/norm( exact_values_H5 )]
nPointsH += [len(eik_coords_H5)]



######################################################
######################################################
######################################################
#### H6
## 1. Plot of the output

eik_vals_H6 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_ComputedValues.bin")
eik_coords_H6 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_MeshPoints.txt", delimiter=",")
triangles_H6 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_Faces.txt", delimiter=",")
eik_grads_H6 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_ComputedGradients.bin");
eik_grads_H6 = eik_grads_H6.reshape(len(eik_coords_H6), 2)

exact_values_H6 = []
errors_H6 = []
for i in range(len(eik_coords_H6)):
    xi_coords = eik_coords_H6[i, 0]
    yi_coords = eik_coords_H6[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H6 += [sol]
    errors_H6 += [ abs( sol - eik_vals_H6[i] ) ]
    


# We interpolate the solution on the triangles_H6 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H6)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

# Contours of the errors_H6 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(ny):
    for j in range(nx):
        tau, type_path = trueSolution(xi[i,j], yi[i,j], x0, center, R, eta1, eta2)
        solution_interpolated[i, j] = tau
errors_H6_abs = abs(zi_linP - solution_interpolated)



#Plot the absolute errors_H6 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im66 = plt.imshow( errors_H6_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H6")
plt.show(block = False)
plt.colorbar(im66)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H6 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6, '-.', lw=0.5, c='#04007e')
im67 = plt.imshow( errors_H6_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H6")
plt.show(block = False)
plt.colorbar(im67)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im68 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], c = eik_vals_H6, cmap = colormap2)
plt.triplot(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H6")
plt.show(block = False)
plt.colorbar(im68)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im69 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H6")
plt.show(block = False)
plt.colorbar(im69)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im610 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H6[:, 0], eik_coords_H6[:, 1], eik_grads_H6[:, 0], eik_grads_H6[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H6")
plt.show(block = False)
plt.colorbar(im610)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H6, triangles_H6)]
errorNorm += [norm( errors_H6  )/norm( exact_values_H6 )]
nPointsH += [len(eik_coords_H6)]


######################################################
######################################################
######################################################
#### H7
## 1. Plot of the output

eik_vals_H7 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_ComputedValues.bin")
eik_coords_H7 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_MeshPoints.txt", delimiter=",")
triangles_H7 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_Faces.txt", delimiter=",")
eik_grads_H7 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_ComputedGradients.bin");
eik_grads_H7 = eik_grads_H7.reshape(len(eik_coords_H7), 2)

exact_values_H7 = []
errors_H7 = []
for i in range(len(eik_coords_H7)):
    xi_coords = eik_coords_H7[i, 0]
    yi_coords = eik_coords_H7[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H7 += [sol]
    errors_H7 += [ abs( sol - eik_vals_H7[i] ) ]
    


# We interpolate the solution on the triangles_H7 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H7[:, 0], eik_coords_H7[:, 1], triangles_H7)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H7)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

# Contours of the errors_H7 in 3D and 2D
solution_interpolated = np.zeros(zi_lin.shape)
for i in range(ny):
    for j in range(nx):
        tau, type_path = trueSolution(xi[i,j], yi[i,j], x0, center, R, eta1, eta2)
        solution_interpolated[i, j] = tau
errors_H7_abs = abs(zi_linP - solution_interpolated)



#Plot the absolute errors_H7 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im76 = plt.imshow( errors_H7_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H7")
plt.show(block = False)
plt.colorbar(im76)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_PointErrors.png', dpi=my_dpi * 10)



# The absolute errors_H7 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H7[:, 0], eik_coords_H7[:, 1], triangles_H7, '-.', lw=0.5, c='#04007e')
im77 = plt.imshow( errors_H7_abs, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H7")
plt.show(block = False)
plt.colorbar(im77)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im78 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H7[:, 0], eik_coords_H7[:, 1], c = eik_vals_H7, cmap = colormap2)
plt.triplot(eik_coords_H7[:, 0], eik_coords_H7[:, 1], triangles_H7, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H7")
plt.show(block = False)
plt.colorbar(im78)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im79 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H7")
plt.show(block = False)
plt.colorbar(im79)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im710 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H7[:, 0], eik_coords_H7[:, 1], eik_grads_H7[:, 0], eik_grads_H7[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H7")
plt.show(block = False)
plt.colorbar(im710)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H7, triangles_H7)]
errorNorm += [norm( errors_H7  )/norm( exact_values_H7 )]
nPointsH += [len(eik_coords_H7)]


######################################################
######################################################
######################################################
######################################################
################## ERRORS ############################
################### EACH #############################
####################  H  #############################
######################################################
######################################################


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH, errorNorm, c = '#6800ff')
plt.title("l2 errors and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_EdgeLength.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH, nPointsH, c = '#6800ff')
plt.title("l2 errors and number of points in triangulation")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_nPoints.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH, times, c = '#6800ff')
plt.title("Average edge length and time taken to solve")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/EdgeLength_Times.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(times, errorNorm, c = '#6800ff')
plt.title("Time taken to solve and l2 errors")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errors.png', dpi=my_dpi * 10)


table_sqTr = {"Average h": averageH, "Time taken": times, "l2 errors": errorNorm, "Points in triangulation": nPointsH}

print(tabulate(table_sqTr, headers="keys", tablefmt="latex"))

plt.show()