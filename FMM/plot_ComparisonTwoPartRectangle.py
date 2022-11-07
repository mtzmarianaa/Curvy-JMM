# Generates error plots for different triangulations - SQUARE WITH ONE TYPE OF INDEX OF REFRACTION


# SCRIPT TO VISUALIZE ERRORS 

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, log, exp, acos, atan, atan2
from tabulate import tabulate
import matplotlib.tri as tri
import pandas as pd
import colorcet as cc
from scipy.optimize import minimize_scalar
import matplotlib.colors as clr


#plt.ion()

colormap2  = clr.LinearSegmentedColormap.from_list('Retro',
                                                   [(0,    '#120c52'),
                                                    (0.25, '#0d0093'),
                                                    (0.60, '#7035c0'),
                                                    (1,    '#e800ff')], N=256)
colormap1 = colormap2

colormap3  = clr.LinearSegmentedColormap.from_list('Retro_div',
                                                   [(0,    '#120c52'),
                                                    (0.5, '#ffffff'),
                                                    (1,    '#e800ff')], N=256)

colormap4  = clr.LinearSegmentedColormap.from_list('Type_sol',
                                                   [(0,    '#00064c'),
                                                    (0.2, '#000eab'),
                                                    (0.4, '#0061ff'),
                                                    (0.6, '#ada8ff'),
                                                    (0.8, '#8700ff'),
                                                    (1,    '#c500ff')], N=256)


saveFigures = False
nx = 10*100
ny = 10*100
my_dpi=96
eta1 = 1.0
eta2 = 1.452
eps = np.finfo(np.float64).resolution

Hs = ["H2", "H3", "H4", "H5"]

#########################################################
####### USEFUL FUNCTIONS

def angle_error( trueGradient, numericalGradient  ):
    '''
    This function calculates the error (angle) between the true gradient and the numerical gradient
    '''
    if ( norm(trueGradient) == 0.0 or norm(numericalGradient) == 0.0  ):
        angle_between = 0.0
    else:
        dProd = np.dot( trueGradient, numericalGradient  )/(  norm(trueGradient)*norm(numericalGradient)  )
        if( dProd<-1.0 or dProd>1.0  ):
            dProd = max( -1.0, min( 1.0, dProd  )  )
        angle_between = acos( dProd  ) # the angle between the two vectors
    return angle_between
    

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
        #queue_indices = liSt(set(queue_indices))
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


def trueSolutionTwoPartSquare(xi, yi, eta1, eta2):
    '''
    Analytic solution for the square with edges on (-5, -5), (5, -5), (5, 5), (-5, 5)
    and two types of indices of refraction separated at y = 2.5
    '''
    if ( yi <= 2.5):
        # We're in the part where it is reachable directly
        sol = sqrt( xi**2 + yi**2 )
        grad = [ eta1*xi/sol, eta1*yi/sol  ]
        return sol, grad
    else:
        # Snells law, we're changing region here, we need an optimization subproblem
        def fToOptimizeUpperRegion(alph):
            p_x, p_y = (1-alph)*(-5) + (alph)*5, 2.5
            return eta1*sqrt( p_x**2 + p_y**2 ) + eta2*sqrt( (p_x - xi)**2 + (p_y - yi)**2  )
        # Solve with scipy
        bounds = (0, 1)
        opt_alpha = minimize_scalar( fToOptimizeUpperRegion, None, bounds, method = 'bounded', tol = eps).x
        sol = fToOptimizeUpperRegion(opt_alpha)
        p_x, p_y = (1-opt_alpha)*(-5) + (opt_alpha)*5, 2.5 # Optimum point of inflection, we need this for the gradient
        d_x, d_y = xi - p_x, yi - p_y
        norm_grad = eta2/sqrt( d_x**2 + d_y**2  )
        grad = [ norm_grad*d_x, norm_grad*d_y  ]
        return sol, grad


    

# Compute the analytic solution in a grid

xi, yi = np.meshgrid(np.linspace(-5, 5, nx), np.linspace(-5, 5, ny))
true_solGrid = np.zeros(xi.shape)



for i in range(ny):
    for j in range(nx):
        true_solGrid[i, j], grad = trueSolutionTwoPartSquare(xi[i, j], yi[i, j], eta1, eta2) 
        

averageH = []
errorl2_eik = []
errorl2_eikIndex1 = []
errorl2_eikIndex2 = []
nPoints = []
times_vec = []
errorl2_grad = []
errorl2_gradIndex1 = []
errorl2_gradIndex2 = []
errorl1_eik = []
errorl1_eikIndex1 = []
errorl1_eikIndex2 = []
errorl1_grad = []
errorl1_gradIndex1 = []
errorl1_gradIndex2 = []
angleError_grad = []
angleError_gradIndex1 = []
angleError_gradIndex2 = []
        
for stringPart in Hs:
    # We want to plot for each of the H's we're considering
    ######
    print(stringPart)

    
    ######
    ######      
    times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_Times.bin")
    eik_vals = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_ComputedValues.bin")
    eik_coords = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_MeshPoints.txt", delimiter=",")
    triangles_points = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_Faces.txt", delimiter=",")
    eik_grads = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_ComputedGradients.bin")
    eik_grads = eik_grads.reshape(len(eik_coords), 2)
    eik_parents = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_Parents.bin", dtype=np.int32)
    eik_parents = eik_parents.reshape(len(eik_coords), 2)
    eik_lambdas = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_LambdasOpt.bin")
    faces_label = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_FacesLabel.txt", dtype=np.int32)
    mesh_tris = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTwoPartSquare/" + stringPart + "/" + stringPart + "_Faces.txt", dtype=np.int32, delimiter = ",")
    mesh_tris = mesh_tris.reshape(len(faces_label), 3)
    
    exact_values_artf = []
    errorsAbs_artf = []
    errors_artf = []
    true_grads = np.zeros(eik_grads.shape)
    errorGradH = []
    normTrueGrads = []
    errorl1GradH = []
    norml1TrueGrads = []
    point_errors_grads = []
    trueAngleGrads = []
    indices_1 = []
    indices_2 = []
    
    for i in range(len(eik_coords)):
        xi_coords = eik_coords[i, 0]
        yi_coords = eik_coords[i, 1]
        solution, gradient = trueSolutionTwoPartSquare(xi_coords, yi_coords, eta1, eta2) 
        sol = solution
        true_grads[i, 0] = gradient[0]
        true_grads[i, 1] = gradient[1]
        normTrueGrads += [ true_grads[i,0]**2 + true_grads[i,1]**2  ]
        norml1TrueGrads += [ norm(true_grads[i, :], 1)   ]
        errorGradH += [ sqrt( (eik_grads[i, 0] - true_grads[i, 0])**2 + (eik_grads[i, 1] - true_grads[ i , 1 ])**2 )]
        errorl1GradH += [ norm( np.subtract(true_grads[i, :], eik_grads[i, :]), 1 )  ]
        exact_values_artf += [sol]
        errorsAbs_artf += [ abs( sol - eik_vals[i] ) ]
        errors_artf += [ sol - eik_vals[i] ]
        point_errors_grads += [ angle_error( true_grads[i, :], eik_grads[i, :]  ) ]
        trueAngleGrads += [ atan2( true_grads[i, 1], true_grads[i, 0]  )   ]
        if(yi_coords<= 2.5):
            indices_1 += [i]
        else:
            indices_2 += [i]
        
    
    #The first one belongs to the source, no gradient there
    errorGradH.pop(0)
    normTrueGrads.pop(0)
    errorl1GradH.pop(0)
    norml1TrueGrads.pop(0)
    point_errors_grads.pop(0)
    trueAngleGrads.pop(0)
    errorGradH = [0] + errorGradH
    normTrueGrads = [1] + normTrueGrads
    errorl1GradH = [0] + errorl1GradH
    norml1TrueGrads = [1] + norml1TrueGrads
    point_errors_grads = [0] + point_errors_grads
    trueAngleGrads = [0] + trueAngleGrads
    # We interpolate the solution on the triangles_points (so that we get a smooth plot + Sam´s idea)
    # We need a triangulation object thing
    triang = tri.Triangulation(eik_coords[:, 0], eik_coords[:, 1], triangles_points)
    # To be able to use LinearTriInterpolator
    interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
    zi_lin = interp_lin(xi, yi)
    zi_linP = interp_lin(xi, yi)
    #Contours of the errorsAbs_artf in 3D and 2D
    errors_inter_artf = true_solGrid - zi_linP
    errorsAbs_inter_artf = abs(true_solGrid - zi_linP )
    vMaxAbs = np.amax(errorsAbs_inter_artf)
    #Plot the absolute errorsAbs_artf in 2D
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    im2_2 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-5, 5, -5, 5], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors, " + stringPart)
    #plt.show(block = False)
    plt.colorbar(im2_2)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_PointErrors.png', dpi=my_dpi * 10)

    # Signed point wise errors
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    im2_3 = plt.imshow( errors_inter_artf, cmap = colormap3, extent=[-5, 5, -5, 5], origin='lower', vmin=-vMaxAbs, vmax=vMaxAbs  )
    plt.title("Signed point wise errors" )
    #plt.show(block = False)
    plt.colorbar(im2_3)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_SignPointErrors.png', dpi=my_dpi * 10)
    # The absolute errorsAbs_artf in 2D with the triangulation

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#ffffff')
    im2_4 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-5, 5, -5, 5], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors and triangulation, " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_4)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_PointErrors_Mesh.png', dpi=my_dpi * 10)

    #Now we can plot + plot the triangulation + dots on top
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    im2_5 = plt.contourf(xi, -yi, zi_lin, cmap = colormap2, levels = 30)
    # plt.scatter(eik_coords[:, 0], eik_coords[:, 1], c = eik_vals, cmap = colormap2)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#6800ff')
    plt.title("Level sets")
    #plt.show(block = False)
    plt.colorbar(im2_5)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_LinearInt_Mesh.png', dpi=my_dpi * 10)
        
    
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    im2_5 = plt.contourf(xi, -yi, zi_lin, cmap = colormap2, levels = 20)
    plt.title("Level sets" )
    #plt.show(block = False)
    plt.colorbar(im2_5)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_LevelSets.png', dpi=my_dpi * 10)
    
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-5, 5, -5, 5], origin='lower'  )
    plt.title("Linear interpolation, " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_LinearInt.png', dpi=my_dpi * 10)

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-5, 5, -5, 5], origin='lower'  )
    plt.quiver(eik_coords[:, 0], eik_coords[:, 1], eik_grads[:, 0], eik_grads[:, 1])
    plt.title("Linear interpolation and computed eikonal gradient, " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_13)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_LinearInt_Grad_A.png', dpi=my_dpi * 10)

    # Plot the errors for the gradients
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    plt.axis('equal')
    #plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
    imError = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = point_errors_grads, cmap = colormap1)
    plt.colorbar(imError)
    plt.title("Angle error in gradients, " + stringPart)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/' + stringPart + "/" + stringPart + '_GradAngleErrors.png', dpi=my_dpi * 10)

    r_iter = range( len( eik_coords  )  )
    averageH += [average_edge_length(eik_coords, triangles_points)]
    errorl2_eik += [norm( errorsAbs_artf  )/norm( exact_values_artf )]
    errorl2_eikIndex1 += [ norm( [ errorsAbs_artf[i] for i in r_iter if i in indices_1 ]  )/norm( [ exact_values_artf[i] for i in r_iter if i in indices_1  ]  )  ]
    errorl2_eikIndex2 += [ norm( [ errorsAbs_artf[i] for i in r_iter if i in indices_2 ]  )/norm( [ exact_values_artf[i] for i in r_iter if i in indices_2  ]  )  ]
    errorl1_eik += [ norm(errorsAbs_artf, 1)/norm( exact_values_artf, 1 )  ]
    errL1ind1 = norm( [ errorsAbs_artf[i] for i in r_iter if i in indices_1], 1)/norm( [ exact_values_artf[i] for i in r_iter if i in indices_1 ], 1  )
    errL1ind2 = norm( [ errorsAbs_artf[i] for i in r_iter if i in indices_2], 1)/norm( [ exact_values_artf[i] for i in r_iter if i in indices_2 ], 1  )
    errorl1_eikIndex1 += [ errL1ind1 ]
    errorl1_eikIndex2 += [ errL1ind2 ]
    nPoints += [len(eik_coords)]
    times_vec += [times[0]]
    errorl2_grad += [ norm(errorGradH) /norm(normTrueGrads) ]
    errGradL2ind1 = norm( [errorGradH[i] for i in r_iter if i in indices_1  ]  )/norm( [ normTrueGrads[i] for i in r_iter if i in indices_1  ]  )
    errGradL2ind2 = norm( [errorGradH[i] for i in r_iter if i in indices_2  ]  )/norm( [ normTrueGrads[i] for i in r_iter if i in indices_2  ]  )
    errorl2_gradIndex1 += [ errGradL2ind1 ]
    errorl2_gradIndex2 += [ errGradL2ind2 ]
    errorl1_grad += [ norm(errorl1GradH, 1)/norm(norml1TrueGrads, 1)    ]
    errGradL1ind1 = norm( [errorl1GradH[i] for i in r_iter if i in indices_1], 1  )/norm([ norml1TrueGrads[i] for i in r_iter if i in indices_1], 1 )
    errGradL1ind2 = norm( [errorl1GradH[i] for i in r_iter if i in indices_2], 1  )/norm([ norml1TrueGrads[i] for i in r_iter if i in indices_2], 1 )
    errorl1_gradIndex1 += [ errGradL1ind1 ]
    errorl1_gradIndex2 += [ errGradL1ind2 ]
    angleError_grad += [ norm(point_errors_grads)/norm(trueAngleGrads)   ]
    errGradAngle1 = norm( [ point_errors_grads[i] for i in r_iter if i in indices_1 ]  )/norm( [ trueAngleGrads[i] for i in r_iter if i in indices_1  ] )
    errGradAngle2 = norm( [ point_errors_grads[i] for i in r_iter if i in indices_2 ]  )/norm( [ trueAngleGrads[i] for i in r_iter if i in indices_2  ] )
    angleError_gradIndex1 += [ errGradAngle1 ]
    angleError_gradIndex2 += [ errGradAngle2 ]

    # Print so that we know its done
    print(stringPart, " done\n")

######################################################
######################################################
######################################################
######################################################
################## ERRORS ############################
################### EACH #############################
####################  H  #############################
######################################################
######################################################


# First we need to order these things so that the plots look nice

info_frameErrors = pd.DataFrame(  data = {'H': Hs,'nPoints': nPoints,'Edge Length': averageH,
                                          'l2 error Eikonal': errorl2_eik, 'l1 error Eikonal': errorl1_eik,
                                          'l2 error Eikonal index1': errorl2_eikIndex1, 'l2 error Eikonal index2': errorl2_eikIndex2,
                                          'l1 error Eikonal index1': errorl1_eikIndex1, 'l1 error Eikonal index2': errorl1_eikIndex2,
                                          'l2 error gradients': errorl2_grad, 'l1 error gradients': errorl1_grad,
                                          'l2 error gradients index1': errorl2_gradIndex1, 'l2 error gradients index2': errorl2_gradIndex2,
                                          'l1 error gradients index1': errorl1_gradIndex1, 'l1 error gradients index2': errorl1_gradIndex2,
                                          'angle error gradients': angleError_grad,
                                          'angle error gradients index1': angleError_gradIndex1, 'angle error gradients index2': angleError_gradIndex2,
                                          'Time to solve (s)': times_vec}  )

# Sort them according to the average edge length original

info_frameErrors = info_frameErrors.sort_values( by = ['Edge Length'], ignore_index = True )
print(info_frameErrors)

# LEAST SQUARES FIT

# LOG OF EVERYTHING
logL2ErrorEik = np.log(info_frameErrors['l2 error Eikonal'])                     # L2 Eikonal
logL2ErrorEik_1 = np.log( info_frameErrors['l2 error Eikonal index1'] )          # L2 Eikonal index 1
logL2ErrorEik_2 = np.log( info_frameErrors['l2 error Eikonal index2'] )          # L2 Eikonal index 2
logL1ErrorEik = np.log(info_frameErrors['l1 error Eikonal'])                     # L1 Eikonal
logL1ErrorEik_1 = np.log( info_frameErrors['l1 error Eikonal index1'] )          # L1 Eikonal index 1
logL1ErrorEik_2 = np.log( info_frameErrors['l1 error Eikonal index2'] )          # L1 Eikonal index 2
logL2ErrorGrad = np.log(info_frameErrors['l2 error gradients'])                  # GRADIENTS ->
logL2ErrorGrad_1 = np.log( info_frameErrors['l2 error gradients index1'] )
logL2ErrorGrad_2 = np.log( info_frameErrors['l2 error gradients index2'] )
logL1ErrorGrad = np.log(info_frameErrors['l1 error gradients'])
logL1ErrorGrad_1 = np.log( info_frameErrors['l1 error gradients index1'] )
logL1ErrorGrad_2 = np.log( info_frameErrors['l1 error gradients index2'] )
logAngleErrorGrad = np.log( info_frameErrors['angle error gradients']  )
logAngleErrorGrad_1 = np.log( info_frameErrors['angle error gradients index1'] )
logAngleErrorGrad_2 = np.log( info_frameErrors['angle error gradients index2'] )
logH = np.log(info_frameErrors['Edge Length'])
logN = np.log(info_frameErrors['nPoints'])

# GET THE POLYNOMIALS AND THEIR COEFFICIENTS
# EIKONAL AND H
logPoly_l2Eik_h = np.polyfit( logH, logL2ErrorEik, deg = 1 )
p_l2Eik_h = np.poly1d(logPoly_l2Eik_h)
logPoly_l2EikInd1_h = np.polyfit( logH, logL2ErrorEik_1, deg = 1 )
p_l2EikInd1_h = np.poly1d(logPoly_l2EikInd1_h)
logPoly_l2EikInd2_h = np.polyfit( logH, logL2ErrorEik_2, deg = 1 )
p_l2EikInd2_h = np.poly1d(logPoly_l2EikInd2_h)

logPoly_l1Eik_h = np.polyfit( logH, logL1ErrorEik, deg = 1 )
p_l1Eik_h = np.poly1d(logPoly_l1Eik_h)
logPoly_l1EikInd1_h = np.polyfit( logH, logL1ErrorEik_1, deg = 1 )
p_l1EikInd1_h = np.poly1d(logPoly_l1EikInd1_h)
logPoly_l1EikInd2_h = np.polyfit( logH, logL1ErrorEik_2, deg = 1 )
p_l1EikInd2_h = np.poly1d(logPoly_l1EikInd2_h)

# EIKONAL AND N
logPoly_l2Eik_n = np.polyfit( logN, logL2ErrorEik, deg = 1 )
p_l2Eik_n = np.poly1d(logPoly_l2Eik_n)
logPoly_l2EikInd1_n = np.polyfit( logN, logL2ErrorEik_1, deg = 1 )
p_l2EikInd1_n = np.poly1d(logPoly_l2EikInd1_n)
logPoly_l2EikInd2_n = np.polyfit( logN, logL2ErrorEik_2, deg = 1 )
p_l2EikInd2_n = np.poly1d(logPoly_l2EikInd2_n)

logPoly_l1Eik_n = np.polyfit( logN, logL1ErrorEik, deg = 1 )
p_l1Eik_n = np.poly1d(logPoly_l1Eik_n)
logPoly_l1EikInd1_n = np.polyfit( logN, logL1ErrorEik_1, deg = 1 )
p_l1EikInd1_n = np.poly1d(logPoly_l1EikInd1_n)
logPoly_l1EikInd2_n = np.polyfit( logN, logL1ErrorEik_2, deg = 1 )
p_l1EikInd2_n = np.poly1d(logPoly_l1EikInd2_n)

# GRADIENT AND H
logPoly_l2Grad_h = np.polyfit( logH, logL2ErrorGrad, deg = 1 )
p_l2Grad_h = np.poly1d(logPoly_l2Grad_h)
logPoly_l2GradInd1_h = np.polyfit( logH, logL2ErrorGrad_1, deg = 1 )
p_l2GradInd1_h = np.poly1d(logPoly_l2GradInd1_h)
logPoly_l2GradInd2_h = np.polyfit( logH, logL2ErrorGrad_2, deg = 1 )
p_l2GradInd2_h = np.poly1d(logPoly_l2GradInd2_h)

logPoly_l1Grad_h = np.polyfit( logH, logL1ErrorGrad, deg = 1 )
p_l1Grad_h = np.poly1d(logPoly_l1Grad_h)
logPoly_l1GradInd1_h = np.polyfit( logH, logL1ErrorGrad_1, deg = 1 )
p_l1GradInd1_h = np.poly1d(logPoly_l1GradInd1_h)
logPoly_l1GradInd2_h = np.polyfit( logH, logL1ErrorGrad_2, deg = 1 )
p_l1GradInd2_h = np.poly1d(logPoly_l1GradInd2_h)

logPoly_angleErr_h = np.polyfit( logH, logAngleErrorGrad, deg = 1 )
p_angleErr_h = np.poly1d(logPoly_angleErr_h)
logPoly_angleErrInd1_h = np.polyfit( logH, logAngleErrorGrad_1, deg = 1 )
p_angleErrInd1_h = np.poly1d(logPoly_angleErrInd1_h)
logPoly_angleErrInd2_h = np.polyfit( logH, logAngleErrorGrad_2, deg = 1 )
p_angleErrInd2_h = np.poly1d(logPoly_angleErrInd2_h)

# GRADIENT AND N
logPoly_l2Grad_n = np.polyfit( logN, logL2ErrorGrad, deg = 1 )
p_l2Grad_n = np.poly1d(logPoly_l2Grad_n)
logPoly_l2GradInd1_n = np.polyfit( logN, logL2ErrorGrad_1, deg = 1 )
p_l2GradInd1_n = np.poly1d(logPoly_l2GradInd1_n)
logPoly_l2GradInd2_n = np.polyfit( logN, logL2ErrorGrad_2, deg = 1 )
p_l2GradInd2_n = np.poly1d(logPoly_l2GradInd2_n)

logPoly_l1Grad_n = np.polyfit( logN, logL1ErrorGrad, deg = 1 )
p_l1Grad_n = np.poly1d(logPoly_l1Grad_n)
logPoly_l1GradInd1_n = np.polyfit( logN, logL1ErrorGrad_1, deg = 1 )
p_l1GradInd1_n = np.poly1d(logPoly_l1GradInd1_n)
logPoly_l1GradInd2_n = np.polyfit( logN, logL1ErrorGrad_2, deg = 1 )
p_l1GradInd2_n = np.poly1d(logPoly_l1GradInd2_n)

logPoly_angleErr_n = np.polyfit( logN, logAngleErrorGrad, deg = 1 )
p_angleErr_n = np.poly1d(logPoly_angleErr_n)
logPoly_angleErrInd1_n = np.polyfit( logN, logAngleErrorGrad_1, deg = 1 )
p_angleErrInd1_n = np.poly1d(logPoly_angleErrInd1_n)
logPoly_angleErrInd2_n = np.polyfit( logN, logAngleErrorGrad_2, deg = 1 )
p_angleErrInd2_n = np.poly1d(logPoly_angleErrInd2_n)

distanceHs = max(info_frameErrors['Edge Length']) - min(info_frameErrors['Edge Length'])
minHToPlot =  min(info_frameErrors['Edge Length']) 
maxHtoPlot =  max(info_frameErrors['Edge Length']) 

distanceNs = max(info_frameErrors['nPoints']) - min(info_frameErrors['nPoints'])
minNToPlot =  min(info_frameErrors['nPoints']) 
maxNtoPlot =  max(info_frameErrors['nPoints']) 

# We put the coefficients in a nice readable table 

table_polynomialCoefs = {'Fitting': ['l2 errors eikonal vs h', 'l1 errors eikonal vs h', 'l2 errors eikonal vs n', 'l1 errors eikonal vs n', 
                                     'l2 errors grad vs h', 'l1 errors grad vs h', 'l2 errors grad vs n', 'l1 errors grad vs n',
                                     'angle error grad vs h', 'angle error grad vs n'],
                        'All': [logPoly_l2Eik_h[0], logPoly_l1Eik_h[0], logPoly_l2Eik_n[0], logPoly_l1Eik_n[0],
                               logPoly_l2Grad_h[0], logPoly_l1Grad_h[0], logPoly_l2Grad_n[0], logPoly_l1Grad_n[0],
                                logPoly_angleErr_h[0], logPoly_angleErr_n[0] ],
                         'Index 1': [logPoly_l2EikInd1_h[0], logPoly_l1EikInd1_h[0], logPoly_l2EikInd1_n[0], logPoly_l1EikInd1_n[0],
                               logPoly_l2GradInd1_h[0], logPoly_l1GradInd1_h[0], logPoly_l2GradInd1_n[0], logPoly_l1GradInd1_n[0],
                                logPoly_angleErrInd1_h[0], logPoly_angleErrInd1_n[0] ],
                         'Index 2': [logPoly_l2EikInd2_h[0], logPoly_l1EikInd2_h[0], logPoly_l2EikInd2_n[0], logPoly_l1EikInd2_n[0],
                               logPoly_l2GradInd2_h[0], logPoly_l1GradInd2_h[0], logPoly_l2GradInd2_n[0], logPoly_l1GradInd2_n[0],
                                logPoly_angleErrInd2_h[0], logPoly_angleErrInd2_n[0] ]}  
info_polynomialCoefs = pd.DataFrame( data = table_polynomialCoefs )
print(info_polynomialCoefs)

# SIMPLE PLOTS

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l2 error Eikonal'], c = '#006ad5', marker='o')
leg = "c1 = " + str( round(logPoly_l2Eik_h[0], 3) ) + "  c0 = " + str( round(logPoly_l2Eik_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l2Eik_h(log(minHToPlot)) ), exp( p_l2Eik_h(log(maxHtoPlot)) )], c = "#0055aa", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors Eikonal and average edge length, square two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/Errors_EdgeLength_artf.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l1 error Eikonal'], c = '#008ade', marker='o')
leg = "c1 = " + str( round(logPoly_l1Eik_h[0], 3) ) + "  c0 = " + str( round(logPoly_l1Eik_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l1Eik_h(log(minHToPlot)) ), exp( p_l1Eik_h(log(maxHtoPlot)) )], c = "#0071b6", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors Eikonal and average edge length, square two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/Errorsl1_EdgeLength_artf.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.plot(info_frameErrors['Edge Length'] , info_frameErrors['l2 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Relative l2 errors Eikonal and average edge length, square one index of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/NOLOG_Errors_EdgeLength_artf.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l2 error Eikonal'], c = '#3e7393', marker='o')
leg = "c1 = " + str( round(logPoly_l2Eik_n[0], 3) ) + "  c0 = " + str( round(logPoly_l2Eik_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l2Eik_n(log(minNToPlot)) ), exp( p_l2Eik_n(log(maxNtoPlot)))], c = "#335e78", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors Eikonal and number of points in triangulation, square two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/Errors_nPoints.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l1 error Eikonal'], c = '#5badd6', marker='o')
leg = "c1 = " + str( round(logPoly_l1Eik_n[0], 3) ) + "  c0 = " + str( round(logPoly_l1Eik_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l1Eik_n(log(minNToPlot))), exp( p_l1Eik_n(log(maxNtoPlot)) )], c = "#4d92b5", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors Eikonal and number of points in triangulation, square two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/Errorsl1_nPoints.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l2 error gradients'], c = '#31008c', marker='o')
leg = "c1 = " + str( round(logPoly_l2Grad_h[0], 3) ) + "  c0 = " + str( round(logPoly_l2Grad_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l2Grad_h(log(minHToPlot)) ), exp( p_l2Grad_h(log(maxHtoPlot)) )], c = "#230065", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors gradients and average edge length, square two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/ErrorsGrad_EdgeLength_artf.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l1 error gradients'], c = '#4500c4', marker='o')
leg = "c1 = " + str( round(logPoly_l1Grad_h[0], 3) ) + "  c0 = " + str( round(logPoly_l1Grad_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l1Grad_h(log(minHToPlot)) ), exp( p_l1Grad_h(log(maxHtoPlot)) )], c = "#3800a1", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors gradients and average edge length, square two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/ErrorsGradl1_EdgeLength_artf.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.plot(info_frameErrors['Edge Length'] , info_frameErrors['l2 error gradients'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Relative l2 errors gradients and average edge length, square two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/NOLOG_ErrorsGrad_EdgeLength_artf.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l2 error gradients'], c = '#543988', marker='o')
leg = "c1 = " + str( round(logPoly_l2Grad_n[0], 3) ) + "  c0 = " + str( round(logPoly_l2Grad_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l2Grad_n(log(minNToPlot)) ), exp( p_l2Grad_n(log(maxNtoPlot)) ) ] , c = "#452f6d", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors gradients and number of points in triangulation, square two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/ErrorsGrad_nPoints.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l1 error gradients'], c = '#754ebd', marker='o')
leg = "c1 = " + str( round(logPoly_l1Grad_n[0], 3) ) + "  c0 = " + str( round(logPoly_l1Grad_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l1Grad_n(log(minNToPlot)) ), exp( p_l1Grad_n(log(maxNtoPlot)) )], c = "#6442a2", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors gradients and number of points in triangulation, square two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/ErrorsGradl1_nPoints.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'], info_frameErrors['angle error gradients'], c = '#754ebd', marker='o')
leg = "c1 = " + str( round(logPoly_angleErr_h[0], 3) ) + "  c0 = " + str( round(logPoly_angleErr_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_angleErr_h(log(minHToPlot)) ), exp( p_angleErr_h(log(maxHtoPlot)) )], c = "#6442a2", linestyle='--', label = leg)
plt.legend()
plt.title("Relative angle errors gradients and average edge length, square two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/AngleErrorsGrad_EdgeLength_artf.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['angle error gradients'], c = '#754ebd', marker='o')
leg = "c1 = " + str( round(logPoly_angleErr_n[0], 3) ) + "  c0 = " + str( round(logPoly_angleErr_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_angleErr_n(log(minNToPlot)) ), exp( p_angleErr_n(log(maxNtoPlot)) )], c = "#6442a2", linestyle='--', label = leg)
plt.legend()
plt.title("Relative angle errors gradients and number of points in triangulation, square two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/AngleErrorsGrad_nPoints.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'], info_frameErrors['Time to solve (s)'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Average edge length and time taken to solve, square two indices of refraction")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/EdgeLength_times.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Time to solve (s)'], info_frameErrors['l2 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l2 errors Eikonal, square two indices of refraction")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/Times_Errors_artf.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Time to solve (s)'], info_frameErrors['l1 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l1 errors Eikonal, square two indices of refraction")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestTwoPartSquare/Times_Errorsl1_artf.png', dpi=my_dpi * 10)


table_artf = {"Average h": info_frameErrors["Edge Length"], "Time taken (s)": info_frameErrors["Time to solve (s)"], 
              "l2 errors Eikonal": info_frameErrors["l2 error Eikonal"], "l1 errors Eikonal": info_frameErrors["l1 error Eikonal"],
              "l2 errors gradients": info_frameErrors["l2 error gradients"], "l1 errors gradients": info_frameErrors["l1 error gradients"],
              "angle error gradients": info_frameErrors['angle error gradients'],
              "Points in triangulation": info_frameErrors["nPoints"]}


print("\n\n\n\nNice table errors\n\n")
print(tabulate(table_artf, headers="keys", tablefmt="latex"))

print("\n\n\n\nNice table coefficients from fitting\n\n")
print(tabulate(table_polynomialCoefs, headers = "keys", tablefmt = "latex"))



#plt.show()
