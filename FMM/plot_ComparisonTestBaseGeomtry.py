# COMPARISON BETWEEN CLASSICAL TRIANGLES AND ARTIFICIAL TRIANGLES


# SCRIPT TO VISUALIZE ERRORS 

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm 
from math import sqrt
import matplotlib.animation as animation
from tabulate import tabulate
from matplotlib.patches import Arc
from analyticSol_circle import trueSolution
import matplotlib.tri as tri
import pandas as pd


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('Spectral_r')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

saveFigures = True
nx = 36*10
ny = 42*10
my_dpi=96
eta1 = 1.0
eta2 = 1.452
x0 = np.array([-15, -10])
center = np.array([0,0])
R = 10
eps = np.finfo(np.float64).resolution

Hs = ["H0", "H-1", "H-2", "H-3", "H-4", "H-5", "H-6", "H1", "H2", "H3", "H4", "H5", "H6", "H7"]

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
        # print("current length of the queue of parents: ", len_queue, "\n")
        for i in range(len_queue):
            lamb = lambdas[ queue_indices[i]]
            sum_x += (1 - lamb)*eik_coords[ queue[i][0] , 0] + lamb*eik_coords[ queue[i][1], 0 ] #add the coordinates of the xlambdas
            sum_y += (1 - lamb)*eik_coords[ queue[i][0] , 1] + lamb*eik_coords[ queue[i][1], 1 ]
            count += 1
        path.extend( [  (sum_x/count, sum_y/count)  ] ) #add the mean of the xlambdas to que path
        # print("Current addition: ", (sum_x/count, sum_y/count))
        queue_indices = [ n[0] for n in queue ] + [ n[1] for n in queue ]
        #queue_indices = list(set(queue_indices))
        if len(set(queue_indices)) == 1:
            if queue_indices[0] ==  queue[0][0] or queue_indices[0] == queue[0][1]: #meaning that its parent is itself
                print(queue_indices)
                print(queue)
                queue = []
            else:
                queue = list(set(queue))
            #path.extend(  [(eik_coords[queue_indices[0], 0], eik_coords[queue_indices[0], 1]  )]   ) # we add the source to the path
        elif len(set(queue_indices)) == 2:
            queue_indices = pd.Series(queue_indices).drop_duplicates().tolist()
            queue = [  (parents_path[ind, 0], parents_path[ind, 1] ) for ind in queue_indices   ]
        else:
            queue_indices = pd.Series(queue_indices).drop_duplicates().tolist()
            queue = [  (parents_path[ind, 0], parents_path[ind, 1] ) for ind in queue_indices   ] # we have a new queue
            queue = pd.Series(queue).drop_duplicates().tolist()
            #queue = list(set(queue))
        if(path[-1] == path[-2]):
            # print("\n\n\nWe started just to add things that were already there")
            # print(queue_indices)
            # print(queue)
            # print(path)
            path.extend( [ (-15, -10) ] )
            return path
    # we need to reverse this list because it starts at the end point and ends in the source (just for aesthetics)
    path.reverse()
    return path
        

# times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/Times_ARTIFICIAL.bin")
# Compute the analytic solution in a grid

xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
true_solGrid = np.zeros(xi.shape)
type_solution = np.zeros(xi.shape)



for i in range(ny):
    for j in range(nx):
        sol, typeSol = trueSolution(  xi[i, j], yi[i,j], x0, center, R, eta1, eta2  )
        true_solGrid[i, j] = sol
        type_solution[i, j] = typeSol
        
averageH_orig = []
errorNorm_orig = []
nPointsH_orig = []
averageH_artf = []
errorNorm_artf = []
nPointsH_artf = []
times_orig_vec = []
times_artf_vec = []
        
for stringPart in Hs:
    # We want to plot for each of the H's we're considering
    ######
    ######       FOR THE ORIGINAL TRIANGLES (IN MESH) UPDATES
    times_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Times.bin")
    eik_vals_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_ComputedValues.bin")
    eik_coords_orig = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_MeshPoints.txt", delimiter=",")
    triangles_orig = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Faces.txt", delimiter=",")
    eik_grads_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_ComputedGradients.bin");
    eik_grads_orig = eik_grads_orig.reshape(len(eik_coords_orig), 2)
    eik_parents_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Parents.bin", dtype=np.int32)
    eik_parents_orig = eik_parents_orig.reshape(len(eik_coords_orig), 2)
    eik_lambdas_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_LambdasOpt.bin")

    exact_values_orig = []
    errorsAbs_orig = []
    errors_orig = []
    for i in range(len(eik_coords_orig)):
        xi_coords = eik_coords_orig[i, 0]
        yi_coords = eik_coords_orig[i, 1]
        sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
        exact_values_orig += [sol]
        errorsAbs_orig += [ abs( sol - eik_vals_orig[i] ) ]
        errors_orig += [ sol - eik_vals_orig[i] ]
    # We interpolate the solution on the triangles_orig (so that we get a smooth plot + Sam´s idea)
    # We need a triangulation object thing
    triang = tri.Triangulation(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig)
    # To be able to use LinearTriInterpolator
    interp_lin = tri.LinearTriInterpolator(triang, eik_vals_orig)
    zi_lin = interp_lin(xi, -yi+6)
    zi_linP = interp_lin(xi, yi)
    #Contours of the errorsAbs_orig in 3D and 2D
    errors_inter_orig = true_solGrid - zi_linP
    errorsAbs_inter_orig = abs(true_solGrid - zi_linP )
    #Plot the absolute errorsAbs_orig in 2D
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_2 = plt.imshow( errorsAbs_inter_orig, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.title("Point wise absolute errors, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_2)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_PointErrors.png', dpi=my_dpi * 10)

    # Signed point wise errors
    vmax = np.max( errorsAbs_inter_orig )
    vmin = -1*vmax
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18) 
    ax.set_ylim(-18, 24)
    im2_3 = plt.imshow( errors_inter_orig, cmap = colormap3, extent=[-18,18,-18,24], origin='lower', vmin = vmin, vmax = vmax  )
    plt.title("Signed point wise absolute errors, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_3)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_SignPointErrors.png', dpi=my_dpi * 10)
    # The absolute errorsAbs_orig in 2D with the triangulation

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#ffffff')
    im2_4 = plt.imshow( errorsAbs_inter_orig, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.title("Point wise absolute errors and triangulation, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_4)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_PointErrors_Mesh.png', dpi=my_dpi * 10)

    #Now we can plot + plot the triangulation + dots on top
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_5 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
    plt.scatter(eik_coords_orig[:, 0], eik_coords_orig[:, 1], c = eik_vals_orig, cmap = colormap2)
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.title("Linear interpolation, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_5)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Mesh.png', dpi=my_dpi * 10)
    
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.title("Linear interpolation, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt.png', dpi=my_dpi * 10)

    # Plotting the paths to certain indices

    ###### Path reg A1
    print("First trial node")
    path1_orig = getPathFromIndex(eik_coords_orig, 1, eik_parents_orig, eik_lambdas_orig)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path1_orig], [p[1] for p in path1_orig], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_orig[[1], 0], eik_coords_orig[[1], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg A1 , test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path1.png', dpi=my_dpi * 10)
    
    ###### Path on circle
    print("Second trial node")
    path2_orig = getPathFromIndex(eik_coords_orig, 6, eik_parents_orig, eik_lambdas_orig)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path2_orig], [p[1] for p in path2_orig], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_orig[[6], 0], eik_coords_orig[[6], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on circle , test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path2.png', dpi=my_dpi * 10)

    ###### Path on reg1
    print("Third trial node")
    path3_orig = getPathFromIndex(eik_coords_orig, 2, eik_parents_orig, eik_lambdas_orig)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path3_orig], [p[1] for p in path3_orig], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_orig[[2], 0], eik_coords_orig[[2], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg 1 , test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path3.png', dpi=my_dpi * 10)

    ###### Path on reg3
    print("Fourth trial node")
    Path4_orig = getPathFromIndex(eik_coords_orig, 3, eik_parents_orig, eik_lambdas_orig)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in Path4_orig], [p[1] for p in Path4_orig], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_orig[[3], 0], eik_coords_orig[[3], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg3 type1 , test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path4.png', dpi=my_dpi * 10)

    ###### Path on reg3 boundary of type 1 and type 2
    print("Fifth trial node")
    Path5_orig = getPathFromIndex(eik_coords_orig, 4, eik_parents_orig, eik_lambdas_orig)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in Path5_orig], [p[1] for p in Path5_orig], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_orig[[4], 0], eik_coords_orig[[4], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg3 boundary type 1/2 , test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path5.png', dpi=my_dpi * 10)

    ###### Path on reg3 type 2
    print("Fifth trial node")
    path6_orig = getPathFromIndex(eik_coords_orig, 5, eik_parents_orig, eik_lambdas_orig)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path6_orig], [p[1] for p in path6_orig], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_orig[[5], 0], eik_coords_orig[[5], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg3 type 2 , test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path6.png', dpi=my_dpi * 10)

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.quiver(eik_coords_orig[:, 0], eik_coords_orig[:, 1], eik_grads_orig[:, 0], eik_grads_orig[:, 1])
    plt.title("Linear interpolation and computed eikonal gradient, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_13)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Grad.png', dpi=my_dpi * 10)

    averageH_orig += [average_edge_length(eik_coords_orig, triangles_orig)]
    errorNorm_orig += [norm( errorsAbs_orig  )/norm( exact_values_orig )]
    nPointsH_orig += [len(eik_coords_orig)]
    times_orig_vec += [times_orig[0]]

    
    ######
    ######       FOR THE UPDATES WITH ARTIFICIAL TRIANGLES
    times_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Times_ARTIFICIAL.bin")
    eik_vals_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_ComputedValues_ARTIFICIAL.bin")
    eik_coords_artf = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_MeshPoints.txt", delimiter=",")
    triangles_artf = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Faces.txt", delimiter=",")
    eik_grads_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_ComputedGradients_ARTIFICIAL.bin");
    eik_grads_artf = eik_grads_artf.reshape(len(eik_coords_artf), 2)
    eik_parents_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Parents_ARTIFICIAL.bin", dtype=np.int32)
    eik_parents_artf = eik_parents_artf.reshape(len(eik_coords_artf), 2)
    eik_lambdas_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_LambdasOpt_ARTIFICIAL.bin")
    
    exact_values_artf = []
    errorsAbs_artf = []
    errors_artf = []
    for i in range(len(eik_coords_artf)):
        xi_coords = eik_coords_artf[i, 0]
        yi_coords = eik_coords_artf[i, 1]
        sol = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
        exact_values_artf += [sol]
        errorsAbs_artf += [ abs( sol - eik_vals_artf[i] ) ]
        errors_artf += [ sol - eik_vals_artf[i] ]
    # We interpolate the solution on the triangles_artf (so that we get a smooth plot + Sam´s idea)
    # We need a triangulation object thing
    triang = tri.Triangulation(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf)
    # To be able to use LinearTriInterpolator
    interp_lin = tri.LinearTriInterpolator(triang, eik_vals_artf)
    zi_lin = interp_lin(xi, -yi+6)
    zi_linP = interp_lin(xi, yi)
    #Contours of the errorsAbs_artf in 3D and 2D
    errors_inter_artf = true_solGrid - zi_linP
    errorsAbs_inter_artf = abs(true_solGrid - zi_linP )
    #Plot the absolute errorsAbs_artf in 2D
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_2 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.title("Point wise absolute errors, test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_2)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_PointErrors_ARTIFICIAL.png', dpi=my_dpi * 10)

    # Signed point wise errors
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18) 
    ax.set_ylim(-18, 24)
    im2_3 = plt.imshow( errors_inter_artf, cmap = colormap3, extent=[-18,18,-18,24], origin='lower', vmin = vmin, vmax = vmax  )
    plt.title("Signed point wise absolute errors, test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_3)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_SignPointErrors_ARTIFICIAL.png', dpi=my_dpi * 10)
    # The absolute errorsAbs_artf in 2D with the triangulation

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#ffffff')
    im2_4 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.title("Point wise absolute errors and triangulation, test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_4)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_PointErrors_Mesh_ARTIFICIAL.png', dpi=my_dpi * 10)

    #Now we can plot + plot the triangulation + dots on top
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_5 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
    plt.scatter(eik_coords_artf[:, 0], eik_coords_artf[:, 1], c = eik_vals_artf, cmap = colormap2)
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.title("Linear interpolation, test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_5)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Mesh_ARTIFICIAL.png', dpi=my_dpi * 10)
    
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.title("Linear interpolation, test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_ARTIFICIAL.png', dpi=my_dpi * 10)

    # Plotting the paths to certain indices

    ###### Path reg A1
    print("First trial node")
    path1_artf = getPathFromIndex(eik_coords_artf, 1, eik_parents_artf, eik_lambdas_artf)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path1_artf], [p[1] for p in path1_artf], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_artf[[1], 0], eik_coords_artf[[1], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg A1 , test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path1_ARTIFICIAL.png', dpi=my_dpi * 10)
    
    ###### Path on circle
    print("Second trial node")
    path2_artf = getPathFromIndex(eik_coords_artf, 6, eik_parents_artf, eik_lambdas_artf)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path2_artf], [p[1] for p in path2_artf], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_artf[[6], 0], eik_coords_artf[[6], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on circle , test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path2_ARTIFICIAL.png', dpi=my_dpi * 10)

    ###### Path on reg1
    print("Third trial node")
    path3_artf = getPathFromIndex(eik_coords_artf, 2, eik_parents_artf, eik_lambdas_artf)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path3_artf], [p[1] for p in path3_artf], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_artf[[2], 0], eik_coords_artf[[2], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg 1 , test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path3_ARTIFICIAL.png', dpi=my_dpi * 10)

    ###### Path on reg3
    print("Fourth trial node")
    Path4_artf = getPathFromIndex(eik_coords_artf, 3, eik_parents_artf, eik_lambdas_artf)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in Path4_artf], [p[1] for p in Path4_artf], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_artf[[3], 0], eik_coords_artf[[3], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg3 type1 , test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path4_ARTIFICIAL.png', dpi=my_dpi * 10)

    ###### Path on reg3 boundary of type 1 and type 2
    print("Fifth trial node")
    Path5_artf = getPathFromIndex(eik_coords_artf, 4, eik_parents_artf, eik_lambdas_artf)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in Path5_artf], [p[1] for p in Path5_artf], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_artf[[4], 0], eik_coords_artf[[4], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg3 boundary type 1/2 , test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path5_ARTIFICIAL.png', dpi=my_dpi * 10)

    ###### Path on reg3 type 2
    print("Fifth trial node")
    path6_artf = getPathFromIndex(eik_coords_artf, 5, eik_parents_artf, eik_lambdas_artf)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca( )
    plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.plot( [p[0] for p in path6_artf], [p[1] for p in path6_artf], marker = ".", c = "#ffffff", linewidth=2 )
    circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
    ax.add_patch(circle_b)
    plt.scatter( eik_coords_artf[[5], 0], eik_coords_artf[[5], 1], marker='o', c = "#c100ff" )
    plt.title("Linear interpolation with path, on reg3 type 2 , test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Path6_ARTIFICIAL.png', dpi=my_dpi * 10)

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.quiver(eik_coords_artf[:, 0], eik_coords_artf[:, 1], eik_grads_artf[:, 0], eik_grads_artf[:, 1])
    plt.title("Linear interpolation and computed eikonal gradient, test geometry just base " + stringPart + "artificial triangles")
    plt.show(block = False)
    plt.colorbar(im2_13)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Grad_ARTIFICIAL.png', dpi=my_dpi * 10)

    averageH_artf += [average_edge_length(eik_coords_artf, triangles_artf)]
    errorNorm_artf += [norm( errorsAbs_artf  )/norm( exact_values_artf )]
    nPointsH_artf += [len(eik_coords_artf)]
    times_artf_vec += [times_artf[0]]
    
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
plt.loglog(averageH_orig, errorNorm_orig, c = '#6800ff', linestyle='--', marker='o')
plt.title("l2 errors and average edge length, triangles in mesht")
plt.xlabel("Average edge length")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_EdgeLength_orig.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_orig, nPointsH_orig, c = '#6800ff', linestyle='--', marker='o')
plt.title("l2 errors and number of points in triangulation, triangles in mesh")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_nPoints_orig.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_orig, times_orig_vec, c = '#6800ff', linestyle='--', marker='o')
plt.title("Average edge length and time taken to solve, triangles in mesh")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/EdgeLength_Times_orig.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(times_orig_vec, errorNorm_orig, c = '#6800ff', linestyle='--', marker='o')
plt.title("Time taken to solve and l2 errors, triangles in mesh")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errors_orig.png', dpi=my_dpi * 10)


table_orig = {"Average h": averageH_orig, "Time taken": times_orig_vec, "l2 errors": errorNorm_orig, "Points in triangulation": nPointsH_orig}

# with artificial triangles

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_artf, errorNorm_artf, c = '#5993b3', linestyle='--', marker='o')
plt.title("l2 errors and average edge length, artificial triangles")
plt.xlabel("Average edge length")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_EdgeLength_artf.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_artf, nPointsH_artf, c = '#5993b3', linestyle='--', marker='o')
plt.title("l2 errors and number of points in triangulation, artificial triangles")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_nPoints_artf.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_artf, times_artf_vec, c = '#5993b3', linestyle='--', marker='o')
plt.title("Average edge length and time taken to solve, artificial triangles")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/EdgeLength_Times_artf.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(times_artf_vec, errorNorm_artf, c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l2 errors, artificial triangles")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errors_artf.png', dpi=my_dpi * 10)


table_artf = {"Average h": averageH_artf, "Time taken": times_artf_vec, "l2 errors": errorNorm_artf, "Points in triangulation": nPointsH_artf}



#### We plot the comparison

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_artf, errorNorm_artf, c = '#5993b3', linestyle='--', marker='o', label = 'Artificial triangles')
plt.loglog(averageH_orig, errorNorm_orig, c = '#6800ff', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.title("l2 errors and average edge length, artificial triangles")
plt.xlabel("Average edge length")
plt.ylabel("Error")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_EdgeLength.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_artf, nPointsH_artf, c = '#5993b3', linestyle='--', marker='o', label = 'Artificial triangles')
plt.loglog(averageH_orig, nPointsH_orig, c = '#6800ff', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.title("l2 errors and number of points in triangulation, artificial triangles")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_nPoints.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(averageH_artf, times_artf_vec, c = '#5993b3', linestyle='--', marker='o', label = 'Artificial triangles')
plt.loglog(averageH_orig, times_orig_vec, c = '#6800ff', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.title("Average edge length and time taken to solve, artificial triangles")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/EdgeLength_Times.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(times_artf_vec, errorNorm_artf, c = '#5993b3', linestyle='--', marker='o', label = 'Artificial triangles')
plt.loglog(times_orig_vec, errorNorm_orig, c = '#6800ff', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.title("Time taken to solve and l2 errors, artificial triangles")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errors.png', dpi=my_dpi * 10)

print("\n\n\nTable with original method:\n")
print(tabulate(table_orig, headers="keys", tablefmt="latex"))

print("\n\n\nTable with artificial triangles:\n")
print(tabulate(table_artf, headers="keys", tablefmt="latex"))

# plt.show()