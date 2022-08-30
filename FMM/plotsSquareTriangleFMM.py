# Script to generate plots from the square with inverted triangle and fmm

# SCRIPT TO VISUALIZE ERRORS (can I say this?)
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm 
from math import sqrt
import matplotlib.tri as tri
from scipy.optimize import NonlinearConstraint, minimize
import matplotlib.animation as animation
import tabulate
import pandas as pd
import colorcet as cc

colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('cet_diverging_cwm_80_100_c22')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

saveFigures = True
nx = 20*10
ny = 20*10
my_dpi=96

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
        s = 1.452
    else:
        s = 1
    return s
        
def exact_solution1(xi, yi):
    '''
    Computation of the exact solution for this exact domain STARTING AT xi = 0, yi = -2
    '''
    x = [xi, yi]
    x0 = [0,0]
    constrain_1 = NonlinearConstraint(constrainOnBoundary, 0, 0)
    if yi >= abs(xi):
        f1 = lambda xlam : IndexRefractionRegions(x0)*sqrt(  xlam[0]**2 + (xlam[1] + 2)**2  ) + IndexRefractionRegions(x)*sqrt(  (xlam[0]-xi)**2 + (xlam[1]-yi)**2 )
        opti_problem1 = minimize( f1, [-10, 10], constraints = constrain_1  )
        opti_problem2 = minimize( f1, [10, 10], constraints = constrain_1  )
        zi = min(f1(opti_problem1.x), f1(opti_problem2.x))
    else:
        zi = sqrt(  xi**2 + (yi + 2)**2  )
    return zi   
        
n = 0
averageH = []
errorNorm = []
nPointsH = []

# times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/Times.bin")

# Compute the analytic solution in a grid

xi, yi = np.meshgrid(np.linspace(-10, 10, nx), np.linspace(-10, 10, ny))
true_solGrid = np.zeros(xi.shape)



for i in range(ny):
    for j in range(nx):
        sol= exact_solution1(  xi[i, j], yi[i,j])
        true_solGrid[i, j] = sol
        
averageH_orig = []
errorNorm_orig = []
nPointsH_orig = []
averageH_artf = []
errorNorm_artf = []
nPointsH_artf = []
times_orig_vec = []
times_artf_vec = []

Hs = ["H0_1","H1_1", "H2_1","H3_1", "H4_1","H5_1", "H6_1","H7_1", "H8_1","H9_1", "H10_1", "H11_1", "H12_1","H13_1", "H14_1","H15_1", "H16_1","H17_1", "H18_1","H19_1", "H20_1"]
# Hs += ["H11_1", "H12_1", "H13_1","H14_1", "H15_1", "H16_1", "H17_1", "H18_1","H19_1", "H20_1"]
        
for stringPart in Hs:
    # We want to plot for each of the H's we're considering
    ######
    ######       FOR THE ORIGINAL TRIANGLES (IN MESH) UPDATES
    times_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_Times.bin")
    eik_vals_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_ComputedValues.bin")
    eik_coords_orig = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_MeshPoints.txt", delimiter=",")
    triangles_orig = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_Faces.txt", delimiter=",", dtype=np.int32)
    eik_grads_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_ComputedGradients.bin");
    eik_grads_orig = eik_grads_orig.reshape(len(eik_coords_orig), 2)
    eik_parents_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_Parents.bin", dtype=np.int32)
    eik_parents_orig = eik_parents_orig.reshape(len(eik_coords_orig), 2)
    eik_lambdas_orig = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_LambdasOpt.bin")

    exact_values_orig = []
    errorsAbs_orig = []
    errors_orig = []
    for i in range(len(eik_coords_orig)):
        xi_coords = eik_coords_orig[i, 0]
        yi_coords = eik_coords_orig[i, 1]
        sol = exact_solution1(xi_coords, yi_coords)
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
    # #Plot the absolute errorsAbs_orig in 2D
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    im2_2 = plt.imshow( errorsAbs_inter_orig, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    plt.title("Point wise absolute errors, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_2)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_PointErrors_AugustC.png', dpi=my_dpi * 10)

    # Signed point wise errors
    vmax = np.max( errorsAbs_inter_orig )
    vmin = -1*vmax
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-10,10) 
    ax.set_ylim(-10,10)
    im2_3 = plt.imshow( errors_inter_orig, cmap = colormap3, extent=[-10,10,-10,10], origin='lower', vmin = vmin, vmax = vmax  )
    plt.title("Signed point wise absolute errors, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_3)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_SignPointErrors_AugustC.png', dpi=my_dpi * 10)
    # The absolute errorsAbs_orig in 2D with the triangulation

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#ffffff')
    im2_4 = plt.imshow( errorsAbs_inter_orig, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    plt.title("Point wise absolute errors and triangulation, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_4)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_PointErrors_Mesh_AugustC.png', dpi=my_dpi * 10)

    #Now we can plot + plot the triangulation + dots on top
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    im2_5 = plt.contourf(xi, -yi, zi_lin, cmap = colormap2)
    plt.scatter(eik_coords_orig[:, 0], eik_coords_orig[:, 1], c = eik_vals_orig, cmap = colormap2)
    plt.triplot(eik_coords_orig[:, 0], eik_coords_orig[:, 1], triangles_orig, '-.', lw=0.2, c='#6800ff')
    plt.title("Linear interpolation, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_5)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_LinearInt_Mesh_AugustC.png', dpi=my_dpi * 10)
    
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    plt.title("Linear interpolation, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_6)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_LinearInt_AugustC.png', dpi=my_dpi * 10)

    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.axis('equal')
    ax = plt.gca()
    ax.set_xlim(-10,10)
    ax.set_ylim(-10,10)
    im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    plt.quiver(eik_coords_orig[:, 0], eik_coords_orig[:, 1], eik_grads_orig[:, 0], eik_grads_orig[:, 1])
    plt.title("Linear interpolation and computed eikonal gradient, test geometry just base " + stringPart + "original updates")
    plt.show(block = False)
    plt.colorbar(im2_13)
    if (saveFigures):
        plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_LinearInt_Grad_AugustC.png', dpi=my_dpi * 10)

    averageH_orig += [average_edge_length(eik_coords_orig, triangles_orig)]
    errorNorm_orig += [norm( errorsAbs_orig  )/norm( exact_values_orig )]
    nPointsH_orig += [len(eik_coords_orig)]
    times_orig_vec += [times_orig[0]]

    
    ######
    ######       FOR THE UPDATES WITH ARTIFICIAL TRIANGLES
    times_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_Times_ARTIFICIAL.bin")
    eik_vals_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_ComputedValues_ARTIFICIAL.bin")
    eik_coords_artf = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_MeshPoints.txt", delimiter=",")
    triangles_artf = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_Faces.txt", delimiter=",")
    eik_grads_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_ComputedGradients_ARTIFICIAL.bin");
    eik_grads_artf = eik_grads_artf.reshape(len(eik_coords_artf), 2)
    eik_parents_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_Parents_ARTIFICIAL.bin", dtype=np.int32)
    eik_parents_artf = eik_parents_artf.reshape(len(eik_coords_artf), 2)
    eik_lambdas_artf = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + stringPart + "/" + stringPart + "_LambdasOpt_ARTIFICIAL.bin")
    
    exact_values_artf = []
    errorsAbs_artf = []
    errors_artf = []
    for i in range(len(eik_coords_artf)):
        xi_coords = eik_coords_artf[i, 0]
        yi_coords = eik_coords_artf[i, 1]
        sol = exact_solution1(xi_coords, yi_coords)
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
    # #Plot the absolute errorsAbs_artf in 2D
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-10,10)
    # ax.set_ylim(-10,10)
    # im2_2 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    # plt.title("Point wise absolute errors, test geometry just base " + stringPart + "artificial triangles")
    # plt.show(block = False)
    # plt.colorbar(im2_2)
    # if (saveFigures):
    #     plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_PointErrors_ARTIFICIAL_AugustC.png', dpi=my_dpi * 10)

    # # Signed point wise errors
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-10,10) 
    # ax.set_ylim(-10,10)
    # im2_3 = plt.imshow( errors_inter_artf, cmap = colormap3, extent=[-10,10,-10,10], origin='lower', vmin = vmin, vmax = vmax  )
    # plt.title("Signed point wise absolute errors, test geometry just base " + stringPart + "artificial triangles")
    # plt.show(block = False)
    # plt.colorbar(im2_3)
    # if (saveFigures):
    #     plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_SignPointErrors_ARTIFICIAL_AugustC.png', dpi=my_dpi * 10)
    # # The absolute errorsAbs_artf in 2D with the triangulation

    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-10,10)
    # ax.set_ylim(-10,10)
    # plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#ffffff')
    # im2_4 = plt.imshow( errorsAbs_inter_artf, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    # plt.title("Point wise absolute errors and triangulation, test geometry just base " + stringPart + "artificial triangles")
    # plt.show(block = False)
    # plt.colorbar(im2_4)
    # if (saveFigures):
    #     plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_PointErrors_Mesh_ARTIFICIAL_AugustC.png', dpi=my_dpi * 10)

    # #Now we can plot + plot the triangulation + dots on top
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-10,10)
    # ax.set_ylim(-10,10)
    # im2_5 = plt.contourf(xi, -yi, zi_lin, cmap = colormap2)
    # plt.scatter(eik_coords_artf[:, 0], eik_coords_artf[:, 1], c = eik_vals_artf, cmap = colormap2)
    # plt.triplot(eik_coords_artf[:, 0], eik_coords_artf[:, 1], triangles_artf, '-.', lw=0.2, c='#6800ff')
    # plt.title("Linear interpolation, test geometry just base " + stringPart + "artificial triangles")
    # plt.show(block = False)
    # plt.colorbar(im2_5)
    # if (saveFigures):
    #     plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_LinearInt_Mesh_ARTIFICIAL_AugustC.png', dpi=my_dpi * 10)
    
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-10,10)
    # ax.set_ylim(-10,10)
    # im2_6 = plt.imshow( zi_linP, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    # plt.title("Linear interpolation, test geometry just base " + stringPart + "artificial triangles")
    # plt.show(block = False)
    # plt.colorbar(im2_6)
    # if (saveFigures):
    #     plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_LinearInt_ARTIFICIAL_AugustC.png', dpi=my_dpi * 10)

    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-10,10)
    # ax.set_ylim(-10,10)
    # im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-10,10,-10,10], origin='lower'  )
    # plt.quiver(eik_coords_artf[:, 0], eik_coords_artf[:, 1], eik_grads_artf[:, 0], eik_grads_artf[:, 1])
    # plt.title("Linear interpolation and computed eikonal gradient, test geometry just base " + stringPart + "artificial triangles")
    # plt.show(block = False)
    # plt.colorbar(im2_13)
    # if (saveFigures):
    #     plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/' + stringPart + "/" + stringPart + '_LinearInt_Grad_ARTIFICIAL_AugustC.png', dpi=my_dpi * 10)

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


# First we need to order these things so that the plots look nice

info_frameErrors = pd.DataFrame(  data = {'H': Hs, 'Edge Length original':  averageH_orig, 'Error original': errorNorm_orig,
                                          'nPoints': nPointsH_orig, 'Times original': times_orig_vec,
                                          'Edge Length artificial': averageH_artf, 'Error artificial': errorNorm_artf,
                                          'Times artificial': times_artf_vec}  )

# Sort them according to the average edge length original

info_frameErrors = info_frameErrors.sort_values( by = ['Edge Length original'], ignore_index = True )

print(info_frameErrors)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length original'] , info_frameErrors['Error original'], c = '#6800ff', linestyle='--', marker='o')
plt.title("l2 errors and average edge length, triangles in mesh")
plt.xlabel("Average edge length")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Errors_EdgeLength_orig_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['Error original'], c = '#6800ff', linestyle='--', marker='o')
plt.title("l2 errors and number of points in triangulation, triangles in mesh")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Errors_nPoints_orig_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length original'], info_frameErrors['Times original'], c = '#6800ff', linestyle='--', marker='o')
plt.title("Average edge length and time taken to solve, triangles in mesh")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/EdgeLength_Times_orig_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Times original'], info_frameErrors['Error original'], c = '#6800ff', linestyle='--', marker='o')
plt.title("Time taken to solve and l2 errors, triangles in mesh")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Times_Errors_orig_AugustC.png', dpi=my_dpi * 10)


table_orig = {"Average h": averageH_orig, "Time taken": times_orig_vec, "l2 errors": errorNorm_orig, "Points in triangulation": nPointsH_orig}

# with artificial triangles

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length artificial'] , info_frameErrors['Error artificial'], c = '#5993b3', linestyle='--', marker='o')
plt.title("l2 errors and average edge length, artificial triangles")
plt.xlabel("Average edge length")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Errors_EdgeLength_artf_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['Error artificial'], c = '#5993b3', linestyle='--', marker='o')
plt.title("l2 errors and number of points in triangulation, artificial triangles")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Errors_nPoints_artf_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length artificial'], info_frameErrors['Times artificial'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Average edge length and time taken to solve, artificial triangles")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/EdgeLength_Times_artf_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Times artificial'], info_frameErrors['Error artificial'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Time taken to solve and l2 errors, artificial triangles")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Times_Errors_artf_AugustC.png', dpi=my_dpi * 10)


table_artf = {"Average h": averageH_artf, "Time taken": times_artf_vec, "l2 errors": errorNorm_artf, "Points in triangulation": nPointsH_artf}



#### We plot the comparison

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length original'] , info_frameErrors['Error original'], c = '#5993b3', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.loglog(info_frameErrors['Edge Length artificial'] , info_frameErrors['Error artificial'], c = '#6800ff', linestyle='--', marker='o', label = 'Artificial triangles')
plt.title("l2 errors and average edge length, artificial triangles")
plt.xlabel("Average edge length")
plt.ylabel("Error")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Errors_EdgeLength_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['Error original'] , c = '#5993b3', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.loglog(info_frameErrors['nPoints'], info_frameErrors['Error artificial'], c = '#6800ff', linestyle='--', marker='o', label = 'Artificial triangles')
plt.title("l2 errors and number of points in triangulation, artificial triangles")
plt.xlabel("Number of points in triangulation")
plt.ylabel("Error")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Errors_nPoints_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length original'] , info_frameErrors['Times original'], c = '#5993b3', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.loglog(info_frameErrors['Edge Length artificial'] , info_frameErrors['Times artificial'], c = '#6800ff', linestyle='--', marker='o', label = 'Artificial triangles')
plt.title("Average edge length and time taken to solve, artificial triangles")
plt.ylabel("Time taken to solve (sec)")
plt.xlabel("Average edge length")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/EdgeLength_Times_AugustC.png', dpi=my_dpi * 10)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Times original'], info_frameErrors['Error original'], c = '#5993b3', linestyle='--', marker='o', label = 'Triangles in mesh')
plt.loglog(info_frameErrors['Times artificial'], info_frameErrors['Error artificial'], c = '#6800ff', linestyle='--', marker='o', label = 'Artificial triangles')
plt.title("Time taken to solve and l2 errors, artificial triangles")
plt.xlabel("Time taken to solve (sec)")
plt.ylabel("Error")
plt.legend()
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/Times_Errors_AugustC.png', dpi=my_dpi * 10)

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.loglog(info_frameErrors['Edge Length original'], info_frameErrors['nPoints'], c = '#5993b3', linestyle='--', marker='o')
plt.title("Points in mesh vs average edge length")
plt.xlabel("Points in mesh")
plt.ylabel("Average edge length")
plt.show(block = False)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/EdgeLength_nPoints_AugustC.png', dpi=my_dpi * 10)

print("\n\n\nTable with original method:\n")
print(tabulate(info_frameErrors, headers="keys", tablefmt="latex"))


