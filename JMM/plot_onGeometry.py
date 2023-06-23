
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

plt.ion()

# Previous useful definitions

colormap1 = plt.cm.get_cmap('BuPu')

colormap2  = "cet_linear_worb_100_25_c53_r"


colormap3  = "cet_diverging_cwm_80_100_c22"


colormap4  = "cet_linear_worb_100_25_c53"

colormap5  = "cet_linear_worb_100_25_c53"


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
                                  true_sol = None, type_sol = None, true_grads = None,
                                  true_solGrid = None, type_solGrid = None, true_gradsGrid = None,
                                  errorsAbs_inter = None, errors_inter = None,
                                  zi_lin = None, point_errors_eik = None, point_errors_grads = None,
                                  saveFigures = False,
                                  show = False, 
                                  path_to_save = '/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/', path_save_vals = "/Users/marianamartinez/Documents/Curvy-JMM/JMM/"):
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
    :param ndarray true_sol: true solution on the mesh
    :param ndarray type_sol: true type of solution on the mesh 
    :param ndarray true_grads: true gradients on the mesh
    :param ndarray true_solGrid: true solution on the xi,yi grid
    :param ndarray type_solGrid: true type of solution on the xi,yi grid
    :param ndarray true_gradsGrid: true gradients on the xi, yi grid
    :param ndarray errorsAbs_inter: absolute point wise errors of the eikonal
    :param ndarray errors_inter: signed point wise errors of the eikonal
    :param ndarray zi_lin: linear interpolation of the computed eikonal to the grid
    :param ndarray point_errors_eik: point wise errors for the eikonal (signed errors) on the mesh
    :param ndarray point_errors_grads: angle error of the gradients on the mesh
    :param bool saveFigures: true if the figures should be saved locally
    :param bool show: if the figures should be shown after computed
    :param str path_to_save: path on which to save the figures generated
    '''
    nPoints = len(eik_vals)
    # Depending on the information given. If no essential information is given we have to compute it
    if x0 is None:
        x0 = x0_default
    if center is None:
        center = center_default
    if R is None:
        R = R_default
    if eta1 is None:
        eta1 = eta1_default
    if eta2 is None:
        eta2 = eta2_default
    if true_solGrid is None and errors_inter is None:
        # IF WE DONT HAVE THE TRUE SOLUTION ON THE GRID AND WE DONT HAVE THE ERRORS
        ny, nx = xi.shape
        true_solGrid = np.zeros(xi.shape)
        true_gradsGrid = np.zeros( (xi.shape[0]*xi.shape[1], 2) )
        for i in range(ny):
            for j in range(nx):
                sol, typeSol, trueGrad = trueSolution( xi[i,j], yi[i,j], x0, center, R, eta1, eta2  )
                true_solGrid[i,j] = sol
    if true_sol is None and true_grads is None:
        # IF WE DONT HAVE THE TRUE SOLUTION ON THE MESH AND WE DONT HAVE THE POINT WISE ERRORS
        nPoints = len(eik_coords) # number of points on the mesh
        true_sol = np.zeros((nPoints))
        true_grads = np.zeros((nPoints, 2))
        point_errors_eik = np.zeros((nPoints))
        point_errors_grads = np.zeros((nPoints))
        for i in range(nPoints):
            sol, typeSol, trueGrad = trueSolution( eik_coords[i, 0], eik_coords[i,1], x0, center, R, eta1, eta2)
            true_sol[i] = sol
            true_grads[i, :] = trueGrad
            point_errors_eik[i] = sol - eik_vals[i]
            point_errors_grads[i] = angle_error( trueGrad, eik_grads[i, :])
        point_errors_grads[0] = 0 # source
        point_errors_eik[0] = 0
        # Save the true_sol and true_grads
        np.savetxt(path_save_vals + H + "/" + "H_trueEiks.txt", true_sol, delimiter = ", ", fmt = '%.0f' )
        np.savetxt(path_save_vals + H + "/" + "H_trueGrads.txt", true_grads, delimiter = ", ", fmt = '%.0f' )
    if point_errors_grads is None and true_grads is not None:
         point_errors_grads = np.zeros((nPoints))
         for i in range(nPoints):
             point_errors_grads[i] = angle_error(true_grads[i], eik_grads[i, :])
         point_errors_grads[0] = 0
    if point_errors_eik is None and true_sol is not None:
        point_errors_eik = np.zeros((nPoints))
        for i in range(nPoints):
            point_errors_eik = true_sol[i] - eik_vals[i]
    # We need the following in order to plot the triangulation
    triang = tri.Triangulation(eik_coords[:,0], eik_coords[:, 1], triangles_points)
    interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
    zi_lin = interp_lin(xi, yi)
    if errors_inter is None:
        # We can now compute the interpolation error since we have the true solution on the grid
        # and we can interpolate linearly the solution from the mesh to the grid
        # We can compute the errors and the absolute point wise errors
        errors_inter = true_solGrid - zi_lin
        errorsAbs_inter = abs( true_solGrid - zi_lin )
        
    # For the errors we need to get the range so that white is plotted exactly in the middle (error = 0)
    vMaxAbs = np.amax(errorsAbs_inter)
    print("Errors mean: ", norm(point_errors_eik)/norm(true_sol))
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
        plt.savefig( path_to_save + H + "/" + H + '_PointErrors.png', dpi=my_dpi * 10)

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
        plt.savefig( path_to_save + H + "/" + H + '_SignPointErrors.png', dpi=my_dpi * 10)

    # ABSOLUTE ERRORS WITH TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#454545')
    im2_3 = plt.imshow( errorsAbs_inter, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors and triangulation")
    plt.colorbar(im2_3)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_PointErrors_Mesh.png', dpi=my_dpi * 10)

    # SIGNED POINT WISE ERRORS WITH TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#454545')
    im2_4 = plt.imshow( errors_inter, cmap = colormap3,  extent=[-18, 18, -18, 24], origin='lower', vmin = -vMaxAbs, vmax = vMaxAbs)
    plt.title("Point wise absolute errors and triangulation")
    plt.colorbar(im2_4)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_SignedPointErrors_Mesh.png', dpi=my_dpi * 10)

    # LEVEL SETS OF COMPUTED SOLUTION (SOLVER)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_5 = plt.contourf(xi, yi, zi_lin, cmap = colormap2, levels = 30)
    plt.contour(xi, yi, zi_lin, colors = ["white"], levels = 30)
    plt.title("Level sets computed solution" )
    plt.colorbar(im2_5)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSets.png', dpi=my_dpi * 10)

    # LEVEL SETS OF EXACT SOLUTION
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_5 = plt.contourf(xi, yi, true_solGrid, cmap = colormap2, levels = 30 )
    plt.contour(xi, yi, true_solGrid, colors = ["white"], levels = 30)
    plt.title("Level sets exact solution" )
    plt.colorbar(im2_5)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSetsExact.png', dpi=my_dpi * 10)

    # LEVEL SETS OF COMPUTED SOLUTION (SOLVER) + TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_6 = plt.contourf(xi, yi, zi_lin, cmap = colormap2, levels = 30 )
    plt.contour(xi, yi, zi_lin, colors = ["white"], levels = 30)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#6800ff')
    plt.title("Level sets computed solution and triangulation")
    plt.colorbar(im2_6)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSets_Mesh.png', dpi=my_dpi * 10)

    # LEVEL SETS OF EXACT SOLUTION  + TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_6 = plt.contourf(xi, yi, true_solGrid, cmap = colormap2, levels = 30)
    plt.contour(xi, yi, true_solGrid, colors = ["white"], levels = 30)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#6800ff')
    plt.title("Level sets exact solution and triangulation")
    plt.colorbar(im2_6)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSets_MeshExact.png', dpi=my_dpi * 10)
    

    # LEVEL SETS OF COMPUTED SOLUTION (SOLVER) + COMPUTED GRADIENTS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_6 = plt.contourf(xi, yi, zi_lin, cmap = colormap2, levels = 30  )
    plt.contour(xi, yi, zi_lin, colors = ["white"], levels = 30)
    plt.quiver(eik_coords[:, 0], eik_coords[:, 1], eik_grads[:, 0], eik_grads[:, 1])
    plt.title("Level sets and computed gradients")
    plt.colorbar(im2_6)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSets_Grad.png', dpi=my_dpi * 10)

    # LEVEL SETS OF EXACT SOLUTION  + EXACT GRADIENTS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_6 = plt.contourf(xi, yi, true_solGrid, cmap = colormap2, levels = 30 )
    plt.contour(xi, yi, true_solGrid, colors = ["white"], levels = 30)
    plt.quiver(eik_coords[:, 0], eik_coords[:, 1], true_grads[:, 0], true_grads[:, 1])
    plt.title("Level sets and exact gradients")
    plt.colorbar(im2_6)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSets_GradExact.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_8 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower'  )
    plt.title("Linear interpolation computed solution")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    plt.colorbar(im2_8)
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LinearInt.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION + TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_9 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18, 18, -18, 24], origin='lower'  )
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#6800ff')
    plt.title("Linear interpolation and triangulation computed solution")
    plt.colorbar(im2_9)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LinearInt_Mesh.png', dpi=my_dpi * 10)

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
        plt.savefig(path_to_save + H + "/" + H + '_LinearInt_Grad.png', dpi=my_dpi * 10)

    # ERRORS IN GRADIENTS
    vmaxGrad = np.max(abs(point_errors_grads))
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    im2_11 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = point_errors_grads, cmap = colormap1, vmin = 0, vmax = vmaxGrad)
    plt.colorbar(im2_11)
    plt.title("Angle error in gradients")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_GradAngleErrors.png', dpi=my_dpi * 10)

    # ERRORS IN GRADIENTS + TRIANGULATION ON TOP
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
    im2_12 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = point_errors_grads, cmap = colormap1, vmin = 0, vmax = vmaxGrad)
    plt.colorbar(im2_12)
    plt.title("Angle error in gradients and triangulation")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_GradAngleErrorsTri.png', dpi=my_dpi * 10)

    # ERRORS IN EIKONAL
    vmaxGrad = np.max(abs(point_errors_grads))
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    im2_11 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = point_errors_eik, cmap = colormap1, vmin = 0, vmax = max(abs(point_errors_eik)))
    plt.colorbar(im2_11)
    plt.title("Point errors Eikonal")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_PointErrEik.png', dpi=my_dpi * 10)

    balEikGrads = abs(point_errors_eik)/np.max(abs(point_errors_eik)) + abs(point_errors_grads)/np.max(abs(point_errors_grads))
    balEikGrads = balEikGrads/2

    # ERRORS IN EIKONAL AND GRAD BALANCED
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#9affe8", lw = 0.3 )
    im2_11 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 5 + round(7500/len(eik_coords)), c = balEikGrads, cmap = colormap1, vmin = 0, vmax = max(abs(balEikGrads)))
    plt.colorbar(im2_11)
    plt.title("Average absolute point wise error Eikonal and gradient")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-18,18])
    ax.set_ylim([-18,24])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_PointErrAverage.png', dpi=my_dpi * 10)


    if(show):
        plt.show()
        
# With this function we can now just define another one that plots everything with just
# defining the H to plot

def plotEverthing_H(H, saveFigures = True):
    '''
    :param str H: name of the triangulation we want to plot everything for
    '''
    path_general = '/Users/marianamartinez/Documents/Curvy-JMM/JMM/'
    path_information = path_general + H + "/" + H 
    path_figures = "/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/"
    xi, yi = np.meshgrid(np.linspace(-18, 18, nx_default), np.linspace(-18, 24, ny_default))
    eik_vals = np.fromfile(path_information + "_ComputedValuesFast.bin")
    eik_coords = np.genfromtxt(path_information + "_MeshPoints.txt", delimiter=",")
    eik_grads = np.fromfile(path_information + "_ComputedGradientsFast.bin")
    eik_grads = eik_grads.reshape(len(eik_coords), 2)
    triangles_points = np.genfromtxt(path_information + "_Faces.txt", delimiter=",")
    true_sol = None
    true_solGrid = None
    true_grads = None
    type_Sol = None
    point_errors_grads = None
    point_errors_eik = None
    errors_inter = None
    errorsAbs_inter = None
    #true_sol = np.genfromtxt(path_information + "_true_values.txt", delimiter = ',')
    #true_solGrid = np.genfromtxt(path_general + "true_solGrid_" + str(nx_default) + "_" + str(ny_default) + ".txt", delimiter = ',')
    #true_grads = np.genfromtxt(path_information + "_true_grads.txt" , delimiter = ',')
    #type_Sol = np.genfromtxt(path_information + "_true_type.txt" , delimiter = ',', dtype = np.int32)
    # We can compute the point_errors_grads
    #point_errors_grads = np.zeros((len(eik_vals)))
    #for i in range(len(eik_coords)):
    #    point_errors_grads[i] =  angle_error( true_grads[i, :], eik_grads[i, :]  )
    #point_errors_grads[0] = 0
    #point_errors_eik = true_sol - eik_vals
    # Now we need to define the triangulation to be able to linearly interpolate
    triang = tri.Triangulation(eik_coords[:,0], eik_coords[:,1], triangles_points)
    interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
    zi_lin = interp_lin(xi, yi)
    # With the triangulation defined we can compute the signed and absolute errors
    #errors_inter = true_solGrid - zi_lin
    #errorsAbs_inter = abs( true_solGrid - zi_lin )

    generatePlotsOnGeometryCircle(H, xi, yi,
                                  eik_vals, eik_coords, eik_grads, triangles_points,
                                  true_sol = true_sol, type_sol = None, true_grads = true_grads,
                                  true_solGrid = true_solGrid, type_solGrid = None, true_gradsGrid = None,
                                  errorsAbs_inter = errorsAbs_inter, errors_inter = errors_inter,
                                  zi_lin = zi_lin, point_errors_eik = point_errors_eik,
                                  point_errors_grads = point_errors_grads,
                                  saveFigures = saveFigures,
                                  show = True, 
                                  path_to_save = "/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/")
