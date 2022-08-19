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
import pandas as pd


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('Spectral_r')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)


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
    j = 0
    path = [(eik_coords[index_to_get_path_from, 0], eik_coords[index_to_get_path_from, 1])] # the path ends/starts in the index of interest
    queue_indices = [index_to_get_path_from]
    queue = [ (parents_path[index_to_get_path_from, 0],  parents_path[index_to_get_path_from, 1]) ] # queue, the indices in the current level considered
    while( queue ): #while there are elements on the list
        # for the pairs of parents in the queue we calculate the mean of the paths (i.e. means of (1-lambda)parent0 + lambda*parent1)
        # print(j)
        # print(queue)
        j += 1
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
        #queue_indices = list(set(queue_indices))
        if len(set(queue_indices)) == 1:
            if queue_indices[0] ==  queue[0][0] and queue_indices[0] == queue[0][1]:
                queue = []
            else:
                queue = list(set(queue))
            #path.extend(  [(eik_coords[queue_indices[0], 0], eik_coords[queue_indices[0], 1]  )]   ) # we add the source to the path
        elif len(set(queue_indices)) == 2:
            queue_indices = pd.Series(queue_indices).drop_duplicates().tolist()
            queue = [  (parents_path[ind, 0], parents_path[ind, 0] ) for ind in queue_indices   ]
        else:
            queue = [  (parents_path[ind, 0], parents_path[ind, 0] ) for ind in queue_indices   ] # we have a new queue
            #queue = list(set(queue))
    # we need to reverse this list because it starts at the end point and ends in the source (just for aesthetics)
    path.reverse()
    print(path)
    return path
        
        
        
    


n = 0
averageH = []
errorNorm = []
nPointsH = []

# times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/Times_ARTIFICIAL_ARTIFICIAL.bin")
# Compute the analytic solution in a grid

xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
true_solGrid = np.zeros(xi.shape)
type_solution = np.zeros(xi.shape)



for i in range(ny):
    for j in range(nx):
        sol, typeSol = trueSolution(  xi[i, j], yi[i,j], x0, center, R, eta1, eta2  )
        true_solGrid[i, j] = sol
        type_solution[i, j] = typeSol

# # We plot the true solution

# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# ax = plt.gca()
# im1 = plt.imshow( true_solGrid, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
# fig.subplots_adjust(wspace=None)
# ax.set_xlim(-18,18)
# ax.set_ylim(-18, 24)
# plt.title("Exact solution, test geometry just base")
# plt.show(block = False)
# plt.colorbar(im1)
# #plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ExactSolution8_ARTIFICIAL.png', dpi=my_dpi * 10)


# # Plot the type of solution
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# ax = plt.gca()
# ax.set_xlim(-18,18)
# ax.set_ylim(-18, 24)
# im1 = plt.imshow( type_solution, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
# plt.title("Type of solution solution, test geometry just base")
# plt.show(block = False)
# plt.colorbar(im1)
# #plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/TypeSolution8_ARTIFICIAL.png', dpi=my_dpi * 10)

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# # # Plot the surface.
# # surf = ax.plot_surface(xi, yi, true_solGrid, cmap=colormap2, linewidth=0, antialiased=False)
# # plt.show(block = False)

# # Plot the contours in 2D
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# ax = plt.gca()
# ax.set_xlim(-18,18)
# ax.set_ylim(-18, 24)
# im_bar1 = plt.contourf(xi, yi, true_solGrid, cmap = colormap2, levels = 25)
# plt.title("Exact solution, test geometry just base")
# plt.show(block = False)
# figName_Contour = '/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ExactSolution_Contour8_ARTIFICIAL.png'
# #plt.savefig(figName_Contour, dpi=my_dpi * 10)

# # # Plot in 3D and save the gif
# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(xi, yi, true_solGrid, c= true_solGrid, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base")
# # plt.show(block = False)
# # # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ExactSolution_Contour8.gif', dpi=80, writer='Pillow')


# # # Save the computed values (in case they are useful)

# # np.savetxt('/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/TestIndex/TrueSolutionGrid.txt', true_solGrid, delimiter =', ', fmt = '%.8f' )



######################################################
######################################################
######################################################
#### H1
## 1. Plot of the output

eik_vals_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedValues_ARTIFICIAL.bin")
eik_coords_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_MeshPoints.txt", delimiter=",")
triangles_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Faces.txt", delimiter=",")
eik_grads_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedGradients_ARTIFICIAL.bin");
eik_grads_H1 = eik_grads_H1.reshape(len(eik_coords_H1), 2)
eik_parents_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Parents_ARTIFICIAL.bin", dtype=np.int32)
eik_parents_H1 = eik_parents_H1.reshape(len(eik_coords_H1), 2)
eik_lambdas_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_LambdasOpt_ARTIFICIAL.bin")

exact_values_H1 = []
errorsAbs_H1 = []
errors_H1 = []
for i in range(len(eik_coords_H1)):
    xi_coords = eik_coords_H1[i, 0]
    yi_coords = eik_coords_H1[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H1 += [sol]
    errorsAbs_H1 += [ abs( sol - eik_vals_H1[i] ) ]
    errors_H1 += [ sol - eik_vals_H1[i] ]


# We interpolate the solution on the triangles_H1 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H1)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

#Contours of the errorsAbs_H1 in 3D and 2D
errors_inter_H1 = true_solGrid - zi_linP
errorsAbs_inter_H1 = abs(true_solGrid - zi_linP )



#Plot the absolute errorsAbs_H1 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_2 = plt.imshow( errorsAbs_inter_H1, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_2)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors_ARTIFICIAL.png', dpi=my_dpi * 10)

# Signed point wise errors
vmax = np.max( errorsAbs_inter_H1 )
vmin = -1*vmax
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18) 
ax.set_ylim(-18, 24)
im2_3 = plt.imshow( errors_inter_H1, cmap = colormap3, extent=[-18,18,-18,24], origin='lower', vmin = vmin, vmax = vmax  )
plt.title("Signed point wise absolute errors, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_3)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_SignPointErrors_ARTIFICIAL.png', dpi=my_dpi * 10)


# The absolute errorsAbs_H1 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#ffffff')
im2_4 = plt.imshow( errorsAbs_inter_H1, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_4)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors_Mesh_ARTIFICIAL.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_5 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], c = eik_vals_H1, cmap = colormap2)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_5)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Mesh_ARTIFICIAL.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_ARTIFICIAL.png', dpi=my_dpi * 10)


# Plotting the paths to certain indices

###### Path reg A1
print("First trial node")
path1_H1 = getPathFromIndex(eik_coords_H1, 1, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in path1_H1], [p[1] for p in path1_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[1], 0], eik_coords_H1[[1], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg A1 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path1_ARTIFICIAL.png', dpi=my_dpi * 10)


###### Path on circle
print("Second trial node")
path2_H1 = getPathFromIndex(eik_coords_H1, 6, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in path2_H1], [p[1] for p in path2_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[6], 0], eik_coords_H1[[6], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on circle , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path2_ARTIFICIAL.png', dpi=my_dpi * 10)


###### Path on reg1
print("Third trial node")
path3_H1 = getPathFromIndex(eik_coords_H1, 2, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in path3_H1], [p[1] for p in path3_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[2], 0], eik_coords_H1[[2], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg 1 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path3_ARTIFICIAL.png', dpi=my_dpi * 10)


###### Path on reg3
print("Fourth trial node")
Path4_H1 = getPathFromIndex(eik_coords_H1, 3, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in Path4_H1], [p[1] for p in Path4_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[3], 0], eik_coords_H1[[3], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg3 type1 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path4_ARTIFICIAL.png', dpi=my_dpi * 10)


###### Path on reg3 boundary of type 1 and type 2
print("Fifth trial node")
Path5_H1 = getPathFromIndex(eik_coords_H1, 4, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in Path5_H1], [p[1] for p in Path5_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[4], 0], eik_coords_H1[[4], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg3 boundary type 1/2 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path5_ARTIFICIAL.png', dpi=my_dpi * 10)


###### Path on reg3 type 2
print("Fifth trial node")
patH1_H1 = getPathFromIndex(eik_coords_H1, 5, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in patH1_H1], [p[1] for p in patH1_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[5], 0], eik_coords_H1[[5], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg3 type 2 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_PatH1_ARTIFICIAL.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_grads_H1[:, 0], eik_grads_H1[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_13)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Grad_ARTIFICIAL.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H1, triangles_H1)]
errorNorm += [norm( errorsAbs_H1  )/norm( exact_values_H1 )]
nPointsH += [len(eik_coords_H1)]

print("Error with artificial triangles", errorNorm[0])




####### Same but with non artificial triangles



eik_vals_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedValues.bin")
eik_coords_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_MeshPoints.txt", delimiter=",")
triangles_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Faces.txt", delimiter=",")
eik_grads_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedGradients.bin");
eik_grads_H1 = eik_grads_H1.reshape(len(eik_coords_H1), 2)
eik_parents_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Parents.bin", dtype=np.int32)
eik_parents_H1 = eik_parents_H1.reshape(len(eik_coords_H1), 2)
eik_lambdas_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_LambdasOpt.bin")

exact_values_H1 = []
errorsAbs_H1 = []
errors_H1 = []
for i in range(len(eik_coords_H1)):
    xi_coords = eik_coords_H1[i, 0]
    yi_coords = eik_coords_H1[i, 1]
    sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
    exact_values_H1 += [sol]
    errorsAbs_H1 += [ abs( sol - eik_vals_H1[i] ) ]
    errors_H1 += [ sol - eik_vals_H1[i] ]


# We interpolate the solution on the triangles_H1 (so that we get a smooth plot + Sam´s idea)


# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H1)
zi_lin = interp_lin(xi, -yi+6)
zi_linP = interp_lin(xi, yi)

#Contours of the errorsAbs_H1 in 3D and 2D
errors_inter_H1 = true_solGrid - zi_linP
errorsAbs_inter_H1 = abs(true_solGrid - zi_linP )



#Plot the absolute errorsAbs_H1 in 2D

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_2 = plt.imshow( errorsAbs_inter_H1, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_2)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors.png', dpi=my_dpi * 10)

# Signed point wise errors
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18) 
ax.set_ylim(-18, 24)
im2_3 = plt.imshow( errors_inter_H1, cmap = colormap3, extent=[-18,18,-18,24], origin='lower', vmin = vmin, vmax = vmax  )
plt.title("Signed point wise absolute errors, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_3)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_SignPointErrors.png', dpi=my_dpi * 10)


# The absolute errorsAbs_H1 in 2D with the triangulation

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#ffffff')
im2_4 = plt.imshow( errorsAbs_inter_H1, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Point wise absolute errors and triangulation, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_4)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_5 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], c = eik_vals_H1, cmap = colormap2)
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.title("Linear interpolation, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_5)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.title("Linear interpolation, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt.png', dpi=my_dpi * 10)


# Plotting the paths to certain indices

###### Path reg A1
print("First trial node")
path1_H1 = getPathFromIndex(eik_coords_H1, 1, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in path1_H1], [p[1] for p in path1_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[1], 0], eik_coords_H1[[1], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg A1 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path1.png', dpi=my_dpi * 10)


###### Path on circle
print("Second trial node")
path2_H1 = getPathFromIndex(eik_coords_H1, 6, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in path2_H1], [p[1] for p in path2_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[6], 0], eik_coords_H1[[6], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on circle , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path2.png', dpi=my_dpi * 10)


###### Path on reg1
print("Third trial node")
path3_H1 = getPathFromIndex(eik_coords_H1, 2, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in path3_H1], [p[1] for p in path3_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[2], 0], eik_coords_H1[[2], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg 1 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path3.png', dpi=my_dpi * 10)


###### Path on reg3
print("Fourth trial node")
Path4_H1 = getPathFromIndex(eik_coords_H1, 3, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in Path4_H1], [p[1] for p in Path4_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[3], 0], eik_coords_H1[[3], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg3 type1 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path4.png', dpi=my_dpi * 10)


###### Path on reg3 boundary of type 1 and type 2
print("Fifth trial node")
Path5_H1 = getPathFromIndex(eik_coords_H1, 4, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in Path5_H1], [p[1] for p in Path5_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[4], 0], eik_coords_H1[[4], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg3 boundary type 1/2 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path5.png', dpi=my_dpi * 10)


###### Path on reg3 type 2
print("Fifth trial node")
Path6_H1 = getPathFromIndex(eik_coords_H1, 5, eik_parents_H1, eik_lambdas_H1)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca( )
plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.2, c='#6800ff')
plt.plot( [p[0] for p in Path6_H1], [p[1] for p in Path6_H1], marker = ".", c = "#ffffff", linewidth=2 )
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.scatter( eik_coords_H1[[5], 0], eik_coords_H1[[5], 1], marker='o', c = "#c100ff" )
plt.title("Linear interpolation with path, on reg3 type 2 , test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_6)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Path6.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
ax = plt.gca()
ax.set_xlim(-18,18)
ax.set_ylim(-18, 24)
im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
plt.quiver(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_grads_H1[:, 0], eik_grads_H1[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H1")
plt.show(block = False)
plt.colorbar(im2_13)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Grad.png', dpi=my_dpi * 10)



averageH += [average_edge_length(eik_coords_H1, triangles_H1)]
errorNorm += [norm( errorsAbs_H1  )/norm( exact_values_H1 )]
nPointsH += [len(eik_coords_H1)]

print("Error without artificial triangles", norm( errorsAbs_H1  )/norm( exact_values_H1 ))


plt.show()