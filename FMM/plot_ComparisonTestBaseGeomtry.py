# COMPARISON BETWEEN CLASSICAL TRIANGLES AND ARTIFICIAL TRIANGLES


# SCRIPT TO VISUALIZE ERRORS 

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm 
from math import sqrt, log, exp
import matplotlib.animation as animation
from tabulate import tabulate
from matplotlib.patches import Arc
from analyticSol_circle import trueSolution
import matplotlib.tri as tri
import pandas as pd
import colorcet as cc

#plt.ion()


colormap1 = plt.cm.get_cmap('plasma')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('cet_diverging_bwg_20_95_c41')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

saveFigures = False
nx = 36*10
ny = 42*10
my_dpi=96
eta1 = 1.0
eta2 = 1.452
x0 = np.array([-15, -10])
center = np.array([0,0])
R = 10
eps = np.finfo(np.float64).resolution

#Hs = ["H1", "H8_5", "H2", "H3", "H4_5", "H3_5", "H4", "H13_5", "H2_5", "H5", "H1_5" , "H12_5" ,"H11_5" ,    "H10_5",  "H14_5", "H15_5", "H16_5", "H17_5"  ,"H6", "H5_5", "H6_5", "H7_5", "H20_5", "H8", "H19_5", "H18_5", "H9_5", "H11",  "H12", "H13", "H14",  "H15" , "H9", "H16", "H17", "H18", "H20", "H21"]
# #     1.75,  1.72,   1.6,  1.4,  1.35,   1.27,   1.2,  1.1,  1.05,   1.0,  0.95,   0.975,   0.9,    0.85,    0.8,    0.75,    0.7,      0.65   ,  0.6,    0.55,   0.5 ,  0.45,   0.4,  0.375,  0.35,  0.325,  0.3,   0.29,  0.28 ,0.275,  0.25, ,0.225,  0.2,  0.19,  0.18,  0.15,  0.05,  0.04
Hs = ["H12","H13","H14","H15","H16", "H9", "H17", "H18", "H20", "H21"]
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
        

# times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/Times.bin")
# Compute the analytic solution in a grid

xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
true_solGrid = np.zeros(xi.shape)
type_solution = np.zeros(xi.shape)



for i in range(ny):
    for j in range(nx):
        sol, typeSol, trueGrad  = trueSolution(  xi[i, j], yi[i,j], x0, center, R, eta1, eta2  )
        true_solGrid[i, j] = sol
        type_solution[i, j] = typeSol
        

averageH_artf = []
errorNorm_artf = []
nPointsH_artf = []
times_artf_vec = []
errorGrad = []
errorl1_eik = []
errorl1_grad = []
        
for stringPart in Hs:
    # We want to plot for each of the H's we're considering
    ######
    print(stringPart)

    
    ######
    ######       FOR THE UPDATES WITH ARTIFICIAL TRIANGLES
    times_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Times.bin")
    eik_vals_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_ComputedValues.bin")
    eik_coords_artf = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_MeshPoints.txt", delimiter=",")
    triangles_artf = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Faces.txt", delimiter=",")
    eik_grads_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_ComputedGradients.bin")
    eik_grads_artf = eik_grads_artf.reshape(len(eik_coords_artf), 2)
    eik_parents_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_Parents.bin", dtype=np.int32)
    eik_parents_artf = eik_parents_artf.reshape(len(eik_coords_artf), 2)
    eik_lambdas_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + stringPart + "/" + stringPart + "_LambdasOpt.bin")
    
    exact_values_artf = []
    errorsAbs_artf = []
    errors_artf = []
    true_grads = np.zeros(eik_grads_artf.shape)
    errorGradH = []
    normTrueGrads = []
    errorl1GradH = []
    norml1TrueGrads = []
    
    for i in range(len(eik_coords_artf)):
        xi_coords = eik_coords_artf[i, 0]
        yi_coords = eik_coords_artf[i, 1]
        sol, typeSolN, trueGradN = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
        true_grads[i, 0] = trueGradN[0]
        true_grads[i, 1] = trueGradN[1]
        normTrueGrads += [ trueGradN[0]**2 + trueGradN[1]**2  ]
        norml1TrueGrads += [ norm(trueGradN, 1)   ]
        errorGradH += [ sqrt( (eik_grads_artf[i, 0] - trueGradN[0])**2 +  (eik_grads_artf[i, 1] - trueGradN[1])**2 )  ]
        errorl1GradH += [ norm( np.subtract(trueGradN, eik_grads_artf[i, :]), 1 )  ]
        exact_values_artf += [sol]
        errorsAbs_artf += [ abs( sol - eik_vals_artf[i] ) ]
        errors_artf += [ sol - eik_vals_artf[i] ]
    
    #The first one belongs to the source, no gradient there
    errorGradH.pop(0)
    normTrueGrads.pop(0)
    errorl1GradH.pop(0)
    norml1TrueGrads.pop(0)
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
    vMaxAbs = np.amax(errorsAbs_inter_artf)
    #Plot the absolute errorsAbs_artf in 2D
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_2 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-18,18,-18,24], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors, test geometry just base " + stringPart)
    #plt.show(block = False)
    plt.colorbar(im2_2)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_PointErrors_D.png', dpi=my_dpi * 10)

    # Signed point wise errors
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18) 
    ax.set_ylim(-18, 24)
    im2_3 = plt.imshow( errors_inter_artf, cmap = colormap3, extent=[-18,18,-18,24], origin='lower', vmin=-vMaxAbs, vmax=vMaxAbs  )
    plt.title("Signed point wise absolute errors, test geometry just base " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_3)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_SignPointErrors_D.png', dpi=my_dpi * 10)
    # The absolute errorsAbs_artf in 2D with the triangulation

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#ffffff')
    im2_4 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-18,18,-18,24], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors and triangulation, test geometry just base " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_4)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_PointErrors_Mesh_D.png', dpi=my_dpi * 10)

    #Now we can plot + plot the triangulation + dots on top
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_5 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2, levels = 30)
    # plt.scatter(eik_coords_artf[:, 0], eik_coords_artf[:, 1], c = eik_vals_artf, cmap = colormap2)
    plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    plt.title("Linear interpolation, test geometry just base " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_5)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Mesh_D.png', dpi=my_dpi * 10)
        
    
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_5 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2, levels = 30)
    plt.title("Linear interpolation, test geometry just base " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_5)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LevelSets_D.png', dpi=my_dpi * 10)
    
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.title("Linear interpolation, test geometry just base " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_D.png', dpi=my_dpi * 10)

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-18,18)
    ax.set_ylim(-18, 24)
    im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    plt.quiver(eik_coords_artf[:, 0], eik_coords_artf[:, 1], eik_grads_artf[:, 0], eik_grads_artf[:, 1])
    plt.title("Linear interpolation and computed eikonal gradient, test geometry just base " + stringPart )
    #plt.show(block = False)
    plt.colorbar(im2_13)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + stringPart + "/" + stringPart + '_LinearInt_Grad_A_D.png', dpi=my_dpi * 10)


    averageH_artf += [average_edge_length(eik_coords_artf, triangles_artf)]
    errorNorm_artf += [norm( errorsAbs_artf  )/norm( exact_values_artf )]
    errorl1_eik += [ norm(errorsAbs_artf, 1)/norm( exact_values_artf, 1 )  ]
    nPointsH_artf += [len(eik_coords_artf)]
    times_artf_vec += [times_artf[0]]
    errorGrad += [ norm(errorGradH) /norm(normTrueGrads) ]
    errorl1_grad += [ norm(errorl1GradH)/norm(norml1TrueGrads)    ]
    plt.close('all')
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

info_frameErrors = pd.DataFrame(  data = {'H': Hs,'nPoints': nPointsH_artf,'Edge Length': averageH_artf, 'l2 error Eikonal': errorNorm_artf, 'l1 error Eikonal': errorl1_eik, 'l2 error gradients': errorGrad, 'l1 error gradients': errorl1_grad,  'Time to solve (s)': times_artf_vec}  )

# Sort them according to the average edge length original

info_frameErrors = info_frameErrors.sort_values( by = ['Edge Length'], ignore_index = True )
print(info_frameErrors)

# LEAST SQUARES FIT

logL2ErrorEik = np.log(info_frameErrors['l2 error Eikonal'])
logL1ErrorEik = np.log(info_frameErrors['l1 error Eikonal'])
logL2ErrorGrad = np.log(info_frameErrors['l2 error gradients'])
logL1ErrorGrad = np.log(info_frameErrors['l1 error gradients'])
logH = np.log(info_frameErrors['Edge Length'])
logN = np.log(info_frameErrors['nPoints'])

logPoly_l2Eik_h = np.polyfit( logH, logL2ErrorEik, deg = 1 )
p_l2Eik_h = np.poly1d(logPoly_l2Eik_h)
logPoly_l1Eik_h = np.polyfit( logH, logL1ErrorEik, deg = 1 )
p_l1Eik_h = np.poly1d(logPoly_l1Eik_h)
logPoly_l2Eik_n = np.polyfit( logN, logL2ErrorEik, deg = 1 )
p_l2Eik_n = np.poly1d(logPoly_l2Eik_n)
logPoly_l1Eik_n = np.polyfit( logN, logL1ErrorEik, deg = 1 )
p_l1Eik_n = np.poly1d(logPoly_l1Eik_n)
logPoly_l2Grad_h = np.polyfit( logH, logL2ErrorGrad, deg = 1 )
p_l2Grad_h = np.poly1d(logPoly_l2Grad_h)
logPoly_l1Grad_h = np.polyfit( logH, logL1ErrorGrad, deg = 1 )
p_l1Grad_h = np.poly1d(logPoly_l1Grad_h)
logPoly_l2Grad_n = np.polyfit( logN, logL2ErrorGrad, deg = 1 )
p_l2Grad_n = np.poly1d(logPoly_l2Grad_n)
logPoly_l1Grad_n = np.polyfit( logN, logL1ErrorGrad, deg = 1 )
p_l1Grad_n = np.poly1d(logPoly_l1Grad_n)

distanceHs = max(info_frameErrors['Edge Length']) - min(info_frameErrors['Edge Length'])
minHToPlot =  min(info_frameErrors['Edge Length']) 
maxHtoPlot =  max(info_frameErrors['Edge Length']) 

distanceNs = max(info_frameErrors['nPoints']) - min(info_frameErrors['nPoints'])
minNToPlot =  min(info_frameErrors['nPoints']) 
maxNtoPlot =  max(info_frameErrors['nPoints']) 

# We put the coefficients in a nice readable table

table_polynomialCoefs = {'Fitting': ['l2 errors eikonal vs h', 'l1 errors eikonal vs h', 'l2 errors eikonal vs n', 'l1 errors eikonal vs n', 
                                    'l2 errors grad vs h', 'l1 errors grad vs h', 'l2 errors grad vs n', 'l1 errors grad vs n'],
                        'c1': [logPoly_l2Eik_h[0], logPoly_l1Eik_h[0], logPoly_l2Eik_n[0], logPoly_l1Eik_n[0],logPoly_l2Grad_h[0], logPoly_l1Grad_h[0], logPoly_l2Grad_n[0], logPoly_l1Grad_n[0]  ],
                        'c0': [logPoly_l2Eik_h[1], logPoly_l1Eik_h[1], logPoly_l2Eik_n[1], logPoly_l1Eik_n[1],logPoly_l2Grad_h[1], logPoly_l1Grad_h[1], logPoly_l2Grad_n[1], logPoly_l1Grad_n[1]  ]}  
info_polynomialCoefs = pd.DataFrame( data = table_polynomialCoefs )
print(info_polynomialCoefs)

# SIMPLE PLOTS

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l2 error Eikonal'], c = '#006ad5', marker='o')
leg = "c1 = " + str( round(logPoly_l2Eik_h[0], 3) ) + "  c0 = " + str( round(logPoly_l2Eik_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l2Eik_h(log(minHToPlot)) ), exp( p_l2Eik_h(log(maxHtoPlot)) )], c = "#0055aa", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors Eikonal and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_EdgeLength_artf_D.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l1 error Eikonal'], c = '#008ade', marker='o')
leg = "c1 = " + str( round(logPoly_l1Eik_h[0], 3) ) + "  c0 = " + str( round(logPoly_l1Eik_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l1Eik_h(log(minHToPlot)) ), exp( p_l1Eik_h(log(maxHtoPlot)) )], c = "#0071b6", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors Eikonal and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errorsl1_EdgeLength_artf_D.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.plot(info_frameErrors['Edge Length'] , info_frameErrors['l2 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Relative l2 errors Eikonal and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/NOLOG_Errors_EdgeLength_artf_D.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l2 error Eikonal'], c = '#3e7393', marker='o')
leg = "c1 = " + str( round(logPoly_l2Eik_n[0], 3) ) + "  c0 = " + str( round(logPoly_l2Eik_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l2Eik_n(log(minNToPlot)) ), exp( p_l2Eik_n(log(maxNtoPlot)))], c = "#335e78", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors Eikonal and number of points in triangulation")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_nPoints_artf_D.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l1 error Eikonal'], c = '#5badd6', marker='o')
leg = "c1 = " + str( round(logPoly_l1Eik_n[0], 3) ) + "  c0 = " + str( round(logPoly_l1Eik_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l1Eik_n(log(minNToPlot))), exp( p_l1Eik_n(log(maxNtoPlot)) )], c = "#4d92b5", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors Eikonal and number of points in triangulation")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errorsl1_nPoints_artf_D.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l2 error gradients'], c = '#31008c', marker='o')
leg = "c1 = " + str( round(logPoly_l2Grad_h[0], 3) ) + "  c0 = " + str( round(logPoly_l2Grad_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l2Grad_h(log(minHToPlot)) ), exp( p_l2Grad_h(log(maxHtoPlot)) )], c = "#230065", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors gradients and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGrad_EdgeLength_artf_D.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l1 error gradients'], c = '#4500c4', marker='o')
leg = "c1 = " + str( round(logPoly_l1Grad_h[0], 3) ) + "  c0 = " + str( round(logPoly_l1Grad_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l1Grad_h(log(minHToPlot)) ), exp( p_l1Grad_h(log(maxHtoPlot)) )], c = "#3800a1", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors gradients and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGradl1_EdgeLength_artf_D.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.plot(info_frameErrors['Edge Length'] , info_frameErrors['l2 error gradients'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Relative l2 errors gradients and average edge length")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/NOLOG_ErrorsGrad_EdgeLength_artf_D.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l2 error gradients'], c = '#543988', marker='o')
leg = "c1 = " + str( round(logPoly_l2Grad_n[0], 3) ) + "  c0 = " + str( round(logPoly_l2Grad_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l2Grad_n(log(minNToPlot)) ), exp( p_l2Grad_n(log(maxNtoPlot)) ) ] , c = "#452f6d", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors gradients and number of points in triangulation")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGrad_nPoints_artf_D.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l1 error gradients'], c = '#754ebd', marker='o')
leg = "c1 = " + str( round(logPoly_l1Grad_n[0], 3) ) + "  c0 = " + str( round(logPoly_l1Grad_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l1Grad_n(log(minNToPlot)) ), exp( p_l1Grad_n(log(maxNtoPlot)) )], c = "#6442a2", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors gradients and number of points in triangulation")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGradl1_nPoints_artf_D.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'], info_frameErrors['Time to solve (s)'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Average edge length and time taken to solve")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/EdgeLength_Times_artf_D.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Time to solve (s)'], info_frameErrors['l2 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l2 errors Eikonal")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errors_artf_D.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Time to solve (s)'], info_frameErrors['l1 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l1 errors Eikonal")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
#plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errorsl1_artf_D.png', dpi=my_dpi * 10)


table_artf = {"Average h": info_frameErrors["Edge Length"], "Time taken (s)": info_frameErrors["Time to solve (s)"], 
              "l2 errors Eikonal": info_frameErrors["l2 error Eikonal"], "l1 errors Eikonal": info_frameErrors["l1 error Eikonal"],
              "l2 errors gradients": info_frameErrors["l2 error gradients"], "l1 errors gradients": info_frameErrors["l1 error gradients"],
              "Points in triangulation": info_frameErrors["nPoints"]}


print("\n\n\n\nNice table errors\n\n")
print(tabulate(table_artf, headers="keys", tablefmt="latex"))

print("\n\n\n\nNice table coefficients from fitting\n\n")
print(tabulate(table_polynomialCoefs, headers = "keys", tablefmt = "latex"))


plt.close('all')
#plt.show()
