
# Similar to generate_BaseSnowTest but here we generate all the plots that go on a geometry, i.e.
# plots that are specific for each mesh. Although there is an option to plot for a list of meshes


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, log, exp, acos, atan
import matplotlib.animation as animation
from tabulate import tabulate
from matplotlib.patches import Arc
from analyticSol_circle import trueSolution
import matplotlib.tri as tri
import pandas as pd
import colorcet as cc
import matplotlib.colors as clr

# Previous useful definitions

colormap1 = plt.cm.get_cmap('BuPu')

colormap2  = clr.LinearSegmentedColormap.from_list('Retro',
                                                   [(0,    '#120c52'),
                                                    (0.25, '#0d0093'),
                                                    (0.60, '#7035c0'),
                                                    (1,    '#e800ff')], N=256)


colormap3  = clr.LinearSegmentedColormap.from_list('Retro_div',
                                                   [(0,    '#120c52'),
                                                    (0.5, '#ffffff'),
                                                    (1,    '#e800ff')], N=256)



spacingGrid = 10
nx_default = 36*spacingGrid
ny_default = 42*spacingGrid
my_dpi=96
eta1_default = 1.0
eta2_default = 1.452
x0_default = np.array([-15, -10])
center_default = np.array([0,0])
R_default = 10.0
eps = np.finfo(np.float64).resolution

# Now useful functions to plot

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

def generatePlotsOnGeometryCircle(H, xi, yi,
                                  eik_vals, eik_coords, eik_grads, triangles_points,
                                  x0 = None, center = None, R = None, eta1 = None, eta2 = None,
                                  true_solGrid = None, type_Sol = None, true_grads = None,
                                  errorsAbs_inter = None, errors_inter = None,
                                  zi_lin = None, point_errors_grads = None,
                                  saveFigures = False,
                                  show = False, 
                                  path_to_save = '/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/'):
    '''
    Function To Generate plots regarding the numerical solution to the Eikonal in the geometry
    with the circle.
    :param str H: name of the mesh to use, including the H
    :param ndarray xi: grid to plot
    :param ndarray yi: grid to plot
    :param ndarray eik_vals: computed Eikonal values on the mesh via the solver
    :param ndarray eik_coords: coordinates of the points in the mesh
    :param ndarray eik_grads: computed gradients on the mesh via the solver
    :param triangles_points: the faces on the mesh
    :param ndarray x0: source point
    :param ndarray center: center of the circle
    :param float R: radius of the circle
    :param float eta1: index of refraction outside the circle
    :param float eta2: index of refraction inside the circle
    :param ndarray true_solGrid: true solution on the xi,yi grid
    :param ndarray typeSol: true type of solution on the xi,yi grid
    :param ndarray true_grads: true gradients on the xi, yi grid
    :param ndarray errorsAbs_inter: absolute point wise errors of the eikonal
    :param ndarray errors_inter: signed point wise errors of the eikonal
    :param ndarray zi_lin: linear interpolation of the computed eikonal to the grid
    :param ndarray point_errors_grads: angle error of the gradients
    :param bool saveFigures: true if the figures should be saved locally
    :param bool show: if the figures should be shown after computed
    :param str path_to_save: path on which to save the figures generated
    '''
    # Depending on the information given. If no essential information is given we have to compute it
    if true_solGrid is None:
        if x0 is None:
            x0 = x0_default
            center = center_default
            R = R_default
            eta1 = eta1_default
            eta2 = eta2_default
        ny, nx = xi.shape
        true_solGrid = np.zeros(xi.shape)
        type_solution = np.zeros(xi.shape)
        for i in range(ny):
            for j in range(nx):
                sol, typeSol, trueGrad = trueSolution( xi[i,j], yi[i,j], x0, center, R, eta1, eta2  )
                true_solGrid[i,j] = sol
                type_solution[i,j] = typeSol
    if errorsAbs_inter is None:
        true_grads = np.zeros(eik_grads.shape)
        point_errors_grads = []
        # We compute the solution at the nodes of the triangulation
        for i in range(len(eik_coords)):
            xi_coords = eik_coords[i,0]
            yi_coords = eik_coords[i,1]
            solution, typeSol, gradient = trueSolution(xi_coords, yi_coords, x0, center, R, eta1, eta2)
            sol = solution
            true_grads[i,0], true_grads[i,1] = gradient[0], gradient[1]
            point_errors_grads += [ angle_error( true_grads[i, :], eik_grads[i, :]  ) ]
        point_errors_grads[0] = 0
        # Now we need to linearly interpolate the solver's output in order to compare it with the true solution on the grid
        triang = tri.Triangulation(eik_coords[:,0], eik_coords[:,1], triangles_points)
        interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
        zi_lin = interp_lin(xi, yi)
        # We can compute the errors and the absolute point wise errors
        errors_inter = true_solGrid - zi_lin
        errorsAbs_inter = abs( true_solGrid - zi_lin )
        
    # For the errors we need to get the range so that white is plotted exactly in the middle (error = 0)
    vMaxAbs = np.amax(errorsAbs_inter)
        
    # Now we have all the date (either it was given or we computed it here)
    # We can now plot

    # ABSOLUTE ERRORS (POINT-WISE)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_1 = plt.imshow( errorsAbs_inter, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors")
    plt.colorbar(im2_1)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig( path_to_save + stringPart + "/" + stringPart + '_PointErrors.png', dpi=my_dpi * 10)

    # SIGNED POINT WISE ERRORS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_2 = plt.imshow( errors_inter, cmap = colormap3, extent=[-18, 18, -18, 24], origin='lower', vmin=-vMaxAbs, vmax=vMaxAbs  )
    plt.title("Signed point wise errors")
    plt.colorbar(im2_2)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig( path_to_save + stringPart + "/" + stringPart + '_SignPointErrors.png', dpi=my_dpi * 10)

    # ABSOLUTE ERRORS WITH TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#ffffff')
    im2_3 = plt.imshow( errorsAbs_inter, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors and triangulation")
    plt.colorbar(im2_3)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_PointErrors_Mesh.png', dpi=my_dpi * 10)

    # SIGNED POINT WISE ERRORS WITH TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#ffffff')
    im2_4 = plt.imshow( errors_inter, cmap = colormap3,  extent=[-18, 18, -18, 24], origin='lower', vmin = -vMaxAbs, vmax = vMaxAbs)
    plt.title("Point wise absolute errors and triangulation")
    plt.colorbar(im2_4)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_SignedPointErrors_Mesh.png', dpi=my_dpi * 10)

    # LEVEL SETS OF COMPUTED SOLUTION (SOLVER)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_5 = plt.contourf(xi, yi, zi_lin, cmap = colormap2, levels = 20)
    plt.title("Level sets" )
    plt.colorbar(im2_5)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_LevelSets.png', dpi=my_dpi * 10)

    # LEVEL SETS OF COMPUTES SOLUTION (SOLVER) + TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_6 = plt.contourf(xi, yi, zi_lin, cmap = colormap2, levels = 30)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#6800ff')
    plt.title("Level sets and triangulation")
    plt.colorbar(im2_6)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_LevelSets_Mesh.png', dpi=my_dpi * 10)

    # LEVEL SETS OF COMPUTES SOLUTION (SOLVER) + COMPUTED GRADIENTS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_6 = plt.contourf(xi, yi, zi_lin, cmap = colormap2, levels = 30)
    plt.quiver(eik_coords[:, 0], eik_coords[:, 1], eik_grads[:, 0], eik_grads[:, 1])
    plt.title("Level sets and computed gradients")
    plt.colorbar(im2_6)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_LevelSets_Grad.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_8 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower'  )
    plt.title("Linear interpolation")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    plt.colorbar(im2_8)
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_LinearInt.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION + TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_9 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower'  )
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#6800ff')
    plt.title("Linear interpolation and triangulation")
    plt.colorbar(im2_9)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_LinearInt_Mesh.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION + COMPUTED GRADIETNS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_10 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower'  )
    plt.quiver(eik_coords[:, 0], eik_coords[:, 1], eik_grads[:, 0], eik_grads[:, 1])
    plt.title("Linear interpolation and computed gradients")
    plt.colorbar(im2_10)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_LinearInt_Grad.png', dpi=my_dpi * 10)

    # ERRORS IN GRADIENTS
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    im2_11 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = point_errors_grads, cmap = colormap2)
    plt.colorbar(im2_10)
    plt.title("Angle error in gradients")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_GradAngleErrors.png', dpi=my_dpi * 10)

    # ERRORS IN GRADIENTS + TRIANGULATION ON TOP
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    #plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
    im2_12 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = point_errors_grads, cmap = colormap2)
    plt.colorbar(im2_12)
    plt.title("Angle error in gradients and triangulation")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + stringPart + "/" + stringPart + '_GradAngleErrors.png', dpi=my_dpi * 10)

    if(show):
        plt.show()
        
