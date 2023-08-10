
# COMPARISON BETWEEN CLASSICAL TRIANGLES AND ARTIFICIAL TRIANGLES


# SCRIPT TO VISUALIZE ERRORS 

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, log, exp, acos, atan, atan2, pi
import matplotlib.animation as animation
from tabulate import tabulate
from matplotlib.patches import Arc
from analyticSol_circle import trueSolution
import matplotlib.tri as tri
import pandas as pd
import colorcet as cc
import matplotlib.colors as clr
from pathlib import Path

# USEFUL THINGS TO CHANGE AND KEEP IN HAND EVERYTIME WE RUN THIS SCRIPT

saveFigures = True
nx = 36*10
ny = 42*10
my_dpi=96
eta1 = 1.0
eta2 = 1.452
x0 = np.array([-15, -10])
center = np.array([0,0])
R = 10.0
eps = np.finfo(np.float64).resolution
# The paths
path_read = '/Users/marianamartinez/Documents/Curvy-JMM/JMM/'
path_trueSol = '/Users/marianamartinez/Documents/Curvy-JMM/JMM/'
#path_trueSol = None

# The different triangulations to consider

Hs = ["H0", "H2", "H4", "H5"]


# Colors
colormap2  = clr.LinearSegmentedColormap.from_list('Retro',
                                                   [(0,    '#120c52'),
                                                    (0.25, '#0d0093'),
                                                    (0.60, '#7035c0'),
                                                    (1,    '#e800ff')], N=256)


colormap3  = clr.LinearSegmentedColormap.from_list('Retro_div',
                                                   [(0,    '#120c52'),
                                                    (0.5, '#ffffff'),
                                                    (1,    '#e800ff')], N=256)



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
    if angle_between > pi:
        return angle_between - pi
    else:
        return angle_between
    

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

        
# Compute the analytic solution in a grid if we don't have it already computed
if path_trueSol is not None:
    path_solution = path_read + 'true_solGrid_' + str(nx) + '_' + str(ny) + '.txt'
    path_type = path_read + 'type_Sol_' + str(nx) + '_' + str(ny) + '.txt'
    path_test1 = Path(path_solution)
    path_test2 = Path(path_type)
else:
    path_test1 = ""
    path_test2 = ""
print(path_test1)
if path_test1.is_file() and path_test2.is_file():
    # Meaning that we gave a path AND the file exists
    true_solGrid = np.genfromtxt(path_solution, delimiter = ',')
    type_Sol =  np.genfromtxt(path_type, dtype=np.int32)
    print("We read files")
    print("Size of true_solGrid", true_solGrid.shape)
    print("Size of type_Sol", type_Sol.shape)
else:
    # If this happens then it means that we haven't computed the true solution in this grid
    # it is useful to save this because it takes a lot of time
    xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
    true_solGrid = np.zeros(xi.shape)
    type_Sol = np.zeros(xi.shape)

    for i in range(ny):
        for j in range(nx):
            sol, typeSol, trueGrad  = trueSolution(  xi[i, j], yi[i,j], x0, center, R, eta1, eta2  )
            true_solGrid[i, j] = sol
            type_Sol[i, j] = typeSol
    # Save them
    np.savetxt(path_read + 'true_solGrid_' + str(nx) + '_' + str(ny) + '.txt', true_solGrid, delimiter = ', ', fmt = '%.8f')
    np.savetxt(path_read + 'type_Sol_' + str(nx) + '_' + str(ny) + '.txt', type_Sol.astype(int), delimiter = ', ', fmt = '%.0f')
        

averageH = []
errorl2_eik = []
errorl2_eikIndex1 = []
errorl2_eikIndex2 = []
errorl2_eikIndex3 = []
nPoints = []
times_vec = []
errorl2_grad = []
errorl2_gradIndex1 = []
errorl2_gradIndex2 = []
errorl2_gradIndex3 = []
errorl1_eik = []
errorl1_eikIndex1 = []
errorl1_eikIndex2 = []
errorl1_eikIndex3 = []
errorl1_grad = []
errorl1_gradIndex1 = []
errorl1_gradIndex2 = []
errorl1_gradIndex3 = []
angleError_grad = []
angleError_gradIndex1 = []
angleError_gradIndex2 = []
angleError_gradIndex3 = []
        
for stringPart in Hs:
    # We want to plot for each of the H's we're considering
    ######
    print(stringPart)
    ######
    ######      READ EVERYTHING
    times = np.fromfile(path_read + stringPart + "/" + stringPart + "_TimesFast.bin")
    eik_vals = np.fromfile(path_read + stringPart + "/" + stringPart + "_ComputedValuesFast.bin")
    eik_coords = np.genfromtxt(path_read + stringPart + "/" + stringPart + "_MeshPoints.txt", delimiter=",")
    triangles_points = np.genfromtxt(path_read + stringPart + "/" + stringPart + "_Faces.txt", delimiter=",")
    eik_grads = np.fromfile(path_read + stringPart + "/" + stringPart + "_ComputedGradientsFast.bin")
    eik_grads = eik_grads.reshape(len(eik_coords), 2)
    faces_label = np.genfromtxt(path_read + stringPart + "/" + stringPart + "_Indices.txt")
    mesh_tris = np.genfromtxt(path_read + stringPart + "/" + stringPart + "_Faces.txt", dtype=np.int32, delimiter = ",")
    mesh_tris = mesh_tris.reshape(len(faces_label), 3)
    
    # We might have already computed the true solution + gradients for this specific mesh, we try to find
    # those files to save up computation time
    path_TrueSolutionMesh = path_read + stringPart + "/" + stringPart + "_true_values.txt"
    path_TrueGradientsMesh = path_read + stringPart + "/" + stringPart + "_true_grads.txt"
    path_TrueTypeSol = path_read + stringPart + "/" + stringPart + "_true_type.txt"
    path_test1 = Path(path_TrueSolutionMesh)
    path_test2 = Path(path_TrueGradientsMesh)
    path_test3 = Path(path_TrueTypeSol)
    if path_test1.is_file() and path_test2.is_file() and path_test3.is_file():
        true_values = np.genfromtxt(path_TrueSolutionMesh, delimiter = ",")
        true_grads = np.genfromtxt(path_TrueGradientsMesh, delimiter = ",")
        true_type = np.genfromtxt(path_TrueTypeSol, dtype = np.int32)
        need_to_compute = False
        print("True solutions for this mesh known, just reading them")
    else:
        true_values = []
        true_grads = np.zeros(eik_grads.shape)
        true_type = []
        need_to_compute = True
    
    errorsAbs = []
    errors = []
    errorGradH = []
    normTrueGrads = []
    errorl1GradH = []
    norml1TrueGrads = []
    point_errors_grads = []
    trueAngleGrads = []
    indices_1 = []
    indices_2 = []
    indices_3 = []
    
    for i in range(len(eik_coords)):
        xi_coords = eik_coords[i, 0]
        yi_coords = eik_coords[i, 1]
        if need_to_compute:
            solution, typeSol, gradient = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
            true_grads[i, 0] = gradient[0]
            true_grads[i, 1] = gradient[1]
            true_values += [solution]
            true_type += [typeSol]
        normTrueGrads += [ true_grads[i,0]**2 + true_grads[i,1]**2  ]
        norml1TrueGrads += [ norm(true_grads[i, :], 1)   ]
        errorGradH += [ sqrt( (eik_grads[i, 0] - true_grads[i, 0])**2 + (eik_grads[i, 1] - true_grads[ i , 1 ])**2 )]
        errorl1GradH += [ norm( np.subtract(true_grads[i, :], eik_grads[i, :]), 1 )  ]
        errorsAbs += [ abs( true_values[i] - eik_vals[i] ) ]
        errors += [ true_values[i] - eik_vals[i] ]
        point_errors_grads += [ angle_error( true_grads[i, :], eik_grads[i, :]  ) ]
        trueAngleGrads += [ atan2( true_grads[i, 1], true_grads[i, 0]  )   ]
        if(true_type[i] == 3):
            indices_1 += [i]
        elif(true_type[i] == 1 or true_type[i] == 2):
            indices_3 += [i]
        else:
            indices_2 += [i]

    # We save the true values, gradients and type of solution if we computed them in this iteration
    if need_to_compute:
        np.savetxt(path_TrueSolutionMesh, true_values, delimiter = ", ", fmt = "%.8f")
        np.savetxt(path_TrueGradientsMesh, true_grads, delimiter = ", ", fmt = "%.8f")
        np.savetxt(path_TrueTypeSol, true_type, delimiter = ", ", fmt = "%.0f")
        print("Saved", path_TrueSolutionMesh)
    
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
    # If we are saving the true values and gradients we do so now

    # To improve performance
    # REGION 1
    errorsAbs_ind1 = []
    errorGradH_ind1 = []
    errorl1GradH_ind1  = []
    point_errors_grads_ind1 = []
    true_values_ind1 = []
    normTrueGrads_ind1 = []
    norml1TrueGrads_ind1 = []
    trueAngleGrads_ind1 = []
    
    for i in indices_1:
        errorsAbs_ind1.append( errorsAbs[i] )
        errorGradH_ind1.append( errorGradH[i] )
        errorl1GradH_ind1.append( errorl1GradH[i] )
        point_errors_grads_ind1.append( point_errors_grads[i] )
        true_values_ind1.append( true_values[i] )
        normTrueGrads_ind1.append( normTrueGrads[i] )
        norml1TrueGrads_ind1.append(  norml1TrueGrads[i] )
        trueAngleGrads_ind1.append(  trueAngleGrads[i] )

    # REGION 2
    errorsAbs_ind2 = []
    errorGradH_ind2 = []
    errorl1GradH_ind2  = []
    point_errors_grads_ind2 = []
    true_values_ind2 = []
    normTrueGrads_ind2 = []
    norml1TrueGrads_ind2 = []
    trueAngleGrads_ind2 = []
    
    for i in indices_2:
        errorsAbs_ind2.append( errorsAbs[i] )
        errorGradH_ind2.append( errorGradH[i] )
        errorl1GradH_ind2.append( errorl1GradH[i] )
        point_errors_grads_ind2.append( point_errors_grads[i] )
        true_values_ind2.append( true_values[i] )
        normTrueGrads_ind2.append( normTrueGrads[i] )
        norml1TrueGrads_ind2.append(  norml1TrueGrads[i] )
        trueAngleGrads_ind2.append(  trueAngleGrads[i] )

    # REGION 3
    errorsAbs_ind3 = []
    errorGradH_ind3 = []
    errorl1GradH_ind3  = []
    point_errors_grads_ind3 = []
    true_values_ind3 = []
    normTrueGrads_ind3 = []
    norml1TrueGrads_ind3 = []
    trueAngleGrads_ind3 = []
    
    for i in indices_3:
        errorsAbs_ind3.append( errorsAbs[i] )
        errorGradH_ind3.append( errorGradH[i] )
        errorl1GradH_ind3.append( errorl1GradH[i] )
        point_errors_grads_ind3.append( point_errors_grads[i] )
        true_values_ind3.append( true_values[i] )
        normTrueGrads_ind3.append( normTrueGrads[i] )
        norml1TrueGrads_ind3.append(  norml1TrueGrads[i] )
        trueAngleGrads_ind3.append(  trueAngleGrads[i] )

    
    # COMPUTE ALL THE ERRORS AND THE VARIABLE NEEDED TO DO SO
    r_iter = range( len( eik_coords  )  )
    averageH += [average_edge_length(eik_coords, triangles_points)]
    errorl2_eik += [norm( errorsAbs  )/norm( true_values )]
    errorl2_eikIndex1 += [ norm( errorsAbs_ind1  )/norm( true_values_ind1  )  ]
    errorl2_eikIndex2 += [ norm( errorsAbs_ind2  )/norm( true_values_ind2  )  ]
    errorl2_eikIndex3 += [ norm( errorsAbs_ind3  )/norm( true_values_ind3  )  ]
    errorl1_eik += [ norm(errorsAbs, 1)/norm( true_values, 1 )  ]
    errL1ind1 = norm(  errorsAbs_ind1, 1)/norm( true_values_ind1, 1  )
    errL1ind2 = norm(  errorsAbs_ind2, 1)/norm( true_values_ind2, 1  )
    errL1ind3 = norm(  errorsAbs_ind3, 1)/norm( true_values_ind3, 1  )
    errorl1_eikIndex1 += [ errL1ind1 ]
    errorl1_eikIndex2 += [ errL1ind2 ]
    errorl1_eikIndex3 += [ errL1ind3 ]
    nPoints += [len(eik_coords)]
    times_vec += [times[0]]
    errorl2_grad += [ norm(errorGradH) /norm(normTrueGrads) ]
    errGradL2ind1 = norm( errorGradH_ind1  )/norm( normTrueGrads_ind1  )
    errGradL2ind2 = norm( errorGradH_ind2  )/norm( normTrueGrads_ind2  )
    errGradL2ind3 = norm( errorGradH_ind3  )/norm( normTrueGrads_ind3  )
    errorl2_gradIndex1 += [ errGradL2ind1 ]
    errorl2_gradIndex2 += [ errGradL2ind2 ]
    errorl2_gradIndex3 += [ errGradL2ind3 ]
    errorl1_grad += [ norm(errorl1GradH, 1)/norm(norml1TrueGrads, 1)    ]
    errGradL1ind1 = norm( errorl1GradH_ind1, 1  )/norm( norml1TrueGrads_ind1, 1 )
    errGradL1ind2 = norm( errorl1GradH_ind2, 1  )/norm( norml1TrueGrads_ind2, 1 )
    errGradL1ind3 = norm( errorl1GradH_ind3, 1  )/norm( norml1TrueGrads_ind3, 1 )
    errorl1_gradIndex1 += [ errGradL1ind1 ]
    errorl1_gradIndex2 += [ errGradL1ind2 ]
    errorl1_gradIndex3 += [ errGradL1ind3 ]
    angleError_grad += [ norm(point_errors_grads)/norm(trueAngleGrads)   ]
    errGradAngle1 = norm( point_errors_grads_ind1  )/norm( trueAngleGrads_ind1 )
    errGradAngle2 = norm( point_errors_grads_ind2  )/norm( trueAngleGrads_ind2 )
    errGradAngle3 = norm( point_errors_grads_ind3  )/norm( trueAngleGrads_ind3 )
    angleError_gradIndex1 += [ errGradAngle1 ]
    angleError_gradIndex2 += [ errGradAngle2 ]
    angleError_gradIndex3 += [ errGradAngle3 ]

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
                                          'l2 error Eikonal index3': errorl2_eikIndex3,
                                          'l1 error Eikonal index1': errorl1_eikIndex1, 'l1 error Eikonal index2': errorl1_eikIndex2,
                                          'l1 error Eikonal index3': errorl1_eikIndex3,
                                          'l2 error gradients': errorl2_grad, 'l1 error gradients': errorl1_grad,
                                          'l2 error gradients index1': errorl2_gradIndex1, 'l2 error gradients index2': errorl2_gradIndex2,
                                          'l2 error gradients index3': errorl2_gradIndex3,
                                          'l1 error gradients index1': errorl1_gradIndex1, 'l1 error gradients index2': errorl1_gradIndex2,
                                          'l1 error gradients index3': errorl1_gradIndex3,
                                          'angle error gradients': angleError_grad,
                                          'angle error gradients index1': angleError_gradIndex1, 'angle error gradients index2': angleError_gradIndex2,
                                          'angle error gradients index3': angleError_gradIndex3,
                                          'Time to solve (s)': times_vec}  )

# Sort them according to the average edge length original

info_frameErrors = info_frameErrors.sort_values( by = ['Edge Length'], ignore_index = True )
print(info_frameErrors)

# LEAST SQUARES FIT

# LOG OF EVERYTHING
logL2ErrorEik = np.log(info_frameErrors['l2 error Eikonal'])                     # L2 Eikonal
logL2ErrorEik_1 = np.log( info_frameErrors['l2 error Eikonal index1'] )          # L2 Eikonal index 1
logL2ErrorEik_2 = np.log( info_frameErrors['l2 error Eikonal index2'] )          # L2 Eikonal index 2
logL2ErrorEik_3 = np.log( info_frameErrors['l2 error Eikonal index3'] )          # L2 Eikonal index 2
logL1ErrorEik = np.log(info_frameErrors['l1 error Eikonal'])                     # L1 Eikonal
logL1ErrorEik_1 = np.log( info_frameErrors['l1 error Eikonal index1'] )          # L1 Eikonal index 1
logL1ErrorEik_2 = np.log( info_frameErrors['l1 error Eikonal index2'] )          # L1 Eikonal index 2
logL1ErrorEik_3 = np.log( info_frameErrors['l1 error Eikonal index3'] )          # L1 Eikonal index 3
logL2ErrorGrad = np.log(info_frameErrors['l2 error gradients'])                  # GRADIENTS ->
logL2ErrorGrad_1 = np.log( info_frameErrors['l2 error gradients index1'] )
logL2ErrorGrad_2 = np.log( info_frameErrors['l2 error gradients index2'] )
logL2ErrorGrad_3 = np.log( info_frameErrors['l2 error gradients index3'] )
logL1ErrorGrad = np.log(info_frameErrors['l1 error gradients'])
logL1ErrorGrad_1 = np.log( info_frameErrors['l1 error gradients index1'] )
logL1ErrorGrad_2 = np.log( info_frameErrors['l1 error gradients index2'] )
logL1ErrorGrad_3 = np.log( info_frameErrors['l1 error gradients index3'] )
logAngleErrorGrad = np.log( info_frameErrors['angle error gradients']  )
logAngleErrorGrad_1 = np.log( info_frameErrors['angle error gradients index1'] )
logAngleErrorGrad_2 = np.log( info_frameErrors['angle error gradients index2'] )
logAngleErrorGrad_3 = np.log( info_frameErrors['angle error gradients index3'] )
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
logPoly_l2EikInd3_h = np.polyfit( logH, logL2ErrorEik_3, deg = 1 )
p_l2EikInd3_h = np.poly1d(logPoly_l2EikInd3_h)

logPoly_l1Eik_h = np.polyfit( logH, logL1ErrorEik, deg = 1 )
p_l1Eik_h = np.poly1d(logPoly_l1Eik_h)
logPoly_l1EikInd1_h = np.polyfit( logH, logL1ErrorEik_1, deg = 1 )
p_l1EikInd1_h = np.poly1d(logPoly_l1EikInd1_h)
logPoly_l1EikInd2_h = np.polyfit( logH, logL1ErrorEik_2, deg = 1 )
p_l1EikInd2_h = np.poly1d(logPoly_l1EikInd2_h)
logPoly_l1EikInd3_h = np.polyfit( logH, logL1ErrorEik_3, deg = 1 )
p_l1EikInd3_h = np.poly1d(logPoly_l1EikInd3_h)

# EIKONAL AND N
logPoly_l2Eik_n = np.polyfit( logN, logL2ErrorEik, deg = 1 )
p_l2Eik_n = np.poly1d(logPoly_l2Eik_n)
logPoly_l2EikInd1_n = np.polyfit( logN, logL2ErrorEik_1, deg = 1 )
p_l2EikInd1_n = np.poly1d(logPoly_l2EikInd1_n)
logPoly_l2EikInd2_n = np.polyfit( logN, logL2ErrorEik_2, deg = 1 )
p_l2EikInd2_n = np.poly1d(logPoly_l2EikInd2_n)
logPoly_l2EikInd3_n = np.polyfit( logN, logL2ErrorEik_3, deg = 1 )
p_l2EikInd3_n = np.poly1d(logPoly_l2EikInd3_n)

logPoly_l1Eik_n = np.polyfit( logN, logL1ErrorEik, deg = 1 )
p_l1Eik_n = np.poly1d(logPoly_l1Eik_n)
logPoly_l1EikInd1_n = np.polyfit( logN, logL1ErrorEik_1, deg = 1 )
p_l1EikInd1_n = np.poly1d(logPoly_l1EikInd1_n)
logPoly_l1EikInd2_n = np.polyfit( logN, logL1ErrorEik_2, deg = 1 )
p_l1EikInd2_n = np.poly1d(logPoly_l1EikInd2_n)
logPoly_l1EikInd3_n = np.polyfit( logN, logL1ErrorEik_3, deg = 1 )
p_l1EikInd3_n = np.poly1d(logPoly_l1EikInd3_n)

# GRADIENT AND H
logPoly_l2Grad_h = np.polyfit( logH, logL2ErrorGrad, deg = 1 )
p_l2Grad_h = np.poly1d(logPoly_l2Grad_h)
logPoly_l2GradInd1_h = np.polyfit( logH, logL2ErrorGrad_1, deg = 1 )
p_l2GradInd1_h = np.poly1d(logPoly_l2GradInd1_h)
logPoly_l2GradInd2_h = np.polyfit( logH, logL2ErrorGrad_2, deg = 1 )
p_l2GradInd2_h = np.poly1d(logPoly_l2GradInd2_h)
logPoly_l2GradInd3_h = np.polyfit( logH, logL2ErrorGrad_3, deg = 1 )
p_l2GradInd3_h = np.poly1d(logPoly_l2GradInd3_h)

logPoly_l1Grad_h = np.polyfit( logH, logL1ErrorGrad, deg = 1 )
p_l1Grad_h = np.poly1d(logPoly_l1Grad_h)
logPoly_l1GradInd1_h = np.polyfit( logH, logL1ErrorGrad_1, deg = 1 )
p_l1GradInd1_h = np.poly1d(logPoly_l1GradInd1_h)
logPoly_l1GradInd2_h = np.polyfit( logH, logL1ErrorGrad_2, deg = 1 )
p_l1GradInd2_h = np.poly1d(logPoly_l1GradInd2_h)
logPoly_l1GradInd3_h = np.polyfit( logH, logL1ErrorGrad_3, deg = 1 )
p_l1GradInd3_h = np.poly1d(logPoly_l1GradInd3_h)

logPoly_angleErr_h = np.polyfit( logH, logAngleErrorGrad, deg = 1 )
p_angleErr_h = np.poly1d(logPoly_angleErr_h)
logPoly_angleErrInd1_h = np.polyfit( logH, logAngleErrorGrad_1, deg = 1 )
p_angleErrInd1_h = np.poly1d(logPoly_angleErrInd1_h)
logPoly_angleErrInd2_h = np.polyfit( logH, logAngleErrorGrad_2, deg = 1 )
p_angleErrInd2_h = np.poly1d(logPoly_angleErrInd2_h)
logPoly_angleErrInd3_h = np.polyfit( logH, logAngleErrorGrad_3, deg = 1 )
p_angleErrInd3_h = np.poly1d(logPoly_angleErrInd3_h)

# GRADIENT AND N
logPoly_l2Grad_n = np.polyfit( logN, logL2ErrorGrad, deg = 1 )
p_l2Grad_n = np.poly1d(logPoly_l2Grad_n)
logPoly_l2GradInd1_n = np.polyfit( logN, logL2ErrorGrad_1, deg = 1 )
p_l2GradInd1_n = np.poly1d(logPoly_l2GradInd1_n)
logPoly_l2GradInd2_n = np.polyfit( logN, logL2ErrorGrad_2, deg = 1 )
p_l2GradInd2_n = np.poly1d(logPoly_l2GradInd2_n)
logPoly_l2GradInd3_n = np.polyfit( logN, logL2ErrorGrad_3, deg = 1 )
p_l2GradInd3_n = np.poly1d(logPoly_l2GradInd3_n)

logPoly_l1Grad_n = np.polyfit( logN, logL1ErrorGrad, deg = 1 )
p_l1Grad_n = np.poly1d(logPoly_l1Grad_n)
logPoly_l1GradInd1_n = np.polyfit( logN, logL1ErrorGrad_1, deg = 1 )
p_l1GradInd1_n = np.poly1d(logPoly_l1GradInd1_n)
logPoly_l1GradInd2_n = np.polyfit( logN, logL1ErrorGrad_2, deg = 1 )
p_l1GradInd2_n = np.poly1d(logPoly_l1GradInd2_n)
logPoly_l1GradInd3_n = np.polyfit( logN, logL1ErrorGrad_3, deg = 1 )
p_l1GradInd3_n = np.poly1d(logPoly_l1GradInd3_n)

logPoly_angleErr_n = np.polyfit( logN, logAngleErrorGrad, deg = 1 )
p_angleErr_n = np.poly1d(logPoly_angleErr_n)
logPoly_angleErrInd1_n = np.polyfit( logN, logAngleErrorGrad_1, deg = 1 )
p_angleErrInd1_n = np.poly1d(logPoly_angleErrInd1_n)
logPoly_angleErrInd2_n = np.polyfit( logN, logAngleErrorGrad_2, deg = 1 )
p_angleErrInd2_n = np.poly1d(logPoly_angleErrInd2_n)
logPoly_angleErrInd3_n = np.polyfit( logN, logAngleErrorGrad_3, deg = 1 )
p_angleErrInd3_n = np.poly1d(logPoly_angleErrInd3_n)

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
                                     logPoly_angleErrInd2_h[0], logPoly_angleErrInd2_n[0] ],
                         'Index 3': [logPoly_l2EikInd3_h[0], logPoly_l1EikInd3_h[0], logPoly_l2EikInd3_n[0], logPoly_l1EikInd3_n[0],
                               logPoly_l2GradInd3_h[0], logPoly_l1GradInd3_h[0], logPoly_l2GradInd3_n[0], logPoly_l1GradInd3_n[0],
                                     logPoly_angleErrInd3_h[0], logPoly_angleErrInd3_n[0] ]}  
info_polynomialCoefs = pd.DataFrame( data = table_polynomialCoefs )
print(info_polynomialCoefs)

# SIMPLE PLOTS

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l2 error Eikonal'], c = '#006ad5', marker='o')
leg = "c1 = " + str( round(logPoly_l2Eik_h[0], 3) ) + "  c0 = " + str( round(logPoly_l2Eik_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l2Eik_h(log(minHToPlot)) ), exp( p_l2Eik_h(log(maxHtoPlot)) )], c = "#0055aa", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors Eikonal and average edge length, circle two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_EdgeLengthCubic.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l1 error Eikonal'], c = '#008ade', marker='o')
leg = "c1 = " + str( round(logPoly_l1Eik_h[0], 3) ) + "  c0 = " + str( round(logPoly_l1Eik_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l1Eik_h(log(minHToPlot)) ), exp( p_l1Eik_h(log(maxHtoPlot)) )], c = "#0071b6", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors Eikonal and average edge length, circle two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errorsl1_EdgeLength.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.plot(info_frameErrors['Edge Length'] , info_frameErrors['l2 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Relative l2 errors Eikonal and average edge length, circle one index of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/NOLOG_Errors_EdgeLengthCubic.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l2 error Eikonal'], c = '#3e7393', marker='o')
leg = "c1 = " + str( round(logPoly_l2Eik_n[0], 3) ) + "  c0 = " + str( round(logPoly_l2Eik_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l2Eik_n(log(minNToPlot)) ), exp( p_l2Eik_n(log(maxNtoPlot)))], c = "#335e78", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors Eikonal and number of points in triangulation, circle two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_nPointsCubic.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l1 error Eikonal'], c = '#5badd6', marker='o')
leg = "c1 = " + str( round(logPoly_l1Eik_n[0], 3) ) + "  c0 = " + str( round(logPoly_l1Eik_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l1Eik_n(log(minNToPlot))), exp( p_l1Eik_n(log(maxNtoPlot)) )], c = "#4d92b5", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors Eikonal and number of points in triangulation, circle two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errorsl1_nPointsCubic.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l2 error gradients'], c = '#31008c', marker='o')
leg = "c1 = " + str( round(logPoly_l2Grad_h[0], 3) ) + "  c0 = " + str( round(logPoly_l2Grad_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l2Grad_h(log(minHToPlot)) ), exp( p_l2Grad_h(log(maxHtoPlot)) )], c = "#230065", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors gradients and average edge length, circle two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGrad_EdgeLengthCubic.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'] , info_frameErrors['l1 error gradients'], c = '#4500c4', marker='o')
leg = "c1 = " + str( round(logPoly_l1Grad_h[0], 3) ) + "  c0 = " + str( round(logPoly_l1Grad_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_l1Grad_h(log(minHToPlot)) ), exp( p_l1Grad_h(log(maxHtoPlot)) )], c = "#3800a1", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors gradients and average edge length, circle two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGradl1_EdgeLengthCubic.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.plot(info_frameErrors['Edge Length'] , info_frameErrors['l2 error gradients'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Relative l2 errors gradients and average edge length, circle two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/NOLOG_ErrorsGrad_EdgeLengthCubic.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l2 error gradients'], c = '#543988', marker='o')
leg = "c1 = " + str( round(logPoly_l2Grad_n[0], 3) ) + "  c0 = " + str( round(logPoly_l2Grad_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l2Grad_n(log(minNToPlot)) ), exp( p_l2Grad_n(log(maxNtoPlot)) ) ] , c = "#452f6d", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l2 errors gradients and number of points in triangulation, circle two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGrad_nPointsCubic.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['l1 error gradients'], c = '#754ebd', marker='o')
leg = "c1 = " + str( round(logPoly_l1Grad_n[0], 3) ) + "  c0 = " + str( round(logPoly_l1Grad_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_l1Grad_n(log(minNToPlot)) ), exp( p_l1Grad_n(log(maxNtoPlot)) )], c = "#6442a2", linestyle='--', label = leg)
plt.legend()
plt.title("Relative l1 errors gradients and number of points in triangulation, circle two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ErrorsGradl1_nPointsCubic.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'], info_frameErrors['angle error gradients'], c = '#754ebd', marker='o')
leg = "c1 = " + str( round(logPoly_angleErr_h[0], 3) ) + "  c0 = " + str( round(logPoly_angleErr_h[1], 3) )
plt.loglog(  [minHToPlot, maxHtoPlot], [exp( p_angleErr_h(log(minHToPlot)) ), exp( p_angleErr_h(log(maxHtoPlot)) )], c = "#6442a2", linestyle='--', label = leg)
plt.legend()
plt.title("Relative angle errors gradients and average edge length, circle two indices of refraction")
plt.xlabel("Average edge length")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/AngleErrorsGrad_EdgeLengthCubic.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['angle error gradients'], c = '#754ebd', marker='o')
leg = "c1 = " + str( round(logPoly_angleErr_n[0], 3) ) + "  c0 = " + str( round(logPoly_angleErr_n[1], 3) )
plt.loglog(  [minNToPlot, maxNtoPlot], [exp( p_angleErr_n(log(minNToPlot)) ), exp( p_angleErr_n(log(maxNtoPlot)) )], c = "#6442a2", linestyle='--', label = leg)
plt.legend()
plt.title("Relative angle errors gradients and number of points in triangulation, circle two indices of refraction")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/AngleErrorsGrad_nPointsCubic.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length'], info_frameErrors['Time to solve (s)'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Average edge length and time taken to solve, circle two indices of refraction")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures//TestBaseSnow/EdgeLength_timesCubic.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Time to solve (s)'], info_frameErrors['l2 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l2 errors Eikonal, circle two indices of refraction")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_ErrorsCubic.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Time to solve (s)'], info_frameErrors['l1 error Eikonal'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l1 errors Eikonal, circle two indices of refraction")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
#plt.show(block = False)
if saveFigures:
    plt.savefig('/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errorsl1Cubic.png', dpi=my_dpi * 10)


table = {"Average h": info_frameErrors["Edge Length"], "Time taken (s)": info_frameErrors["Time to solve (s)"], 
              "l2 errors Eikonal": info_frameErrors["l2 error Eikonal"], "l1 errors Eikonal": info_frameErrors["l1 error Eikonal"],
              "l2 errors gradients": info_frameErrors["l2 error gradients"], "l1 errors gradients": info_frameErrors["l1 error gradients"],
              "angle error gradients": info_frameErrors['angle error gradients'],
              "Points in triangulation": info_frameErrors["nPoints"]}


print("\n\n\n\nNice table errors\n\n")
print(tabulate(table, headers="keys", tablefmt="latex"))

print("\n\n\n\nNice table coefficients from fitting\n\n")
print(tabulate(table_polynomialCoefs, headers = "keys", tablefmt = "latex"))



plt.show()






