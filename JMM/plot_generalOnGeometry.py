
# Given a path of computed eikonals for a general domain, plot the results


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, log, exp, acos, atan
import matplotlib.tri as tri
import colorcet as cc
import matplotlib.colors as clr
from scipy.optimize import minimize_scalar
from pathlib import Path

plt.ion()


path_information = "/Users/marianamartinez/Documents/Curvy-JMM/JMM/EasyGeom/"
path_to_save =  "/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/EasyGeom/"

# Previous useful definitions

colormap1 = plt.cm.get_cmap('BuPu')

colormap2  = "cet_linear_worb_100_25_c53_r"


#colormap3 = "cet_diverging_linear_bjr_30_55_c53"
colormap3 = clr.LinearSegmentedColormap.from_list('Barbie_div',
                                                   [(0,    '#0008ff'),
                                                    (0.5, '#ffffff'),
                                                    (1,    '#ff00ec')], N=256)



colormap4  = "cet_linear_worb_100_25_c53"

colormap5  = "cet_linear_worb_100_25_c53"


spacingGrid = 10
nx_default = 36*spacingGrid
ny_default = 42*spacingGrid
my_dpi=96
eta1_default = 1.0
eta2_default = 1.0
x0_default = np.array([0, -1.5])
eps = np.finfo(np.float64).resolution




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


def optiLine(lam, objPoint):
    xOnLine = np.array( [lam, 0] )
    return eta1_default*norm( x0_default - xOnLine ) + eta2_default*norm( xOnLine - objPoint )


def solutionLineTwoIndices(points, x0 = x0_default):
    '''
    Solution to a square divided by a line - two segments
    '''
    direct_ind = np.where(points[:, 1] <= 0)
    direct_ind = direct_ind[0]
    trueEik = np.zeros( (len(points), 1))
    trueGrads = np.zeros( (len(points), 2))
    trueEik[direct_ind, 0] = eta1_default*norm( points[direct_ind, :] - x0, axis = 1)
    trueGrads[direct_ind, :] = eta1_default*((points[direct_ind, :] - x0)/trueEik[direct_ind])
    indirect_ind = np.where( points[:, 1] > 0)
    indirect_ind = indirect_ind[0]
    for i in indirect_ind:
        objPoint = points[i, :]
        tMin = lambda lam: optiLine(lam, objPoint)
        lamMin = minimize_scalar(tMin, None, [-2, 2], method = 'bounded').x
        trueEik[i, 0] = tMin(lamMin)
        infPoint = np.array( [lamMin, 0])
        trueGrads[i, :] = eta2_default*( (objPoint - infPoint)/norm(objPoint - infPoint)   )
    zInd = np.where( trueEik == 0 )
    zInd = zInd[0]
    trueGrads[zInd] = 0
    return trueEik, trueGrads
    

def solutionLineOneIndex(points, x0 = x0_default):
    '''
    Solution to a square divided by a line - two segments
    '''
    trueEik = np.zeros( (len(points), 1))
    trueEik[:, 0] = norm( points - x0, axis = 1)
    trueGrads = eta1_default*((points - x0)/trueEik)
    zInd = np.where( trueEik == 0 )
    zInd = zInd[0]
    trueGrads[zInd] = 0
    return trueEik, trueGrads


def generatePlotsOnGeometry(H, type = 0, x0 = None, show = True, saveFigures = True, path_information = path_information, eta1 = None, eta2 = None):
    '''
    Function To Generate plots regarding the numerical solution to the Eikonal in a given domain
    '''
    path_information = path_information + H + "/" + H
    xi, yi = np.meshgrid(np.linspace(-2, 2, nx_default), np.linspace(-2, 2, ny_default))
    eik_vals = np.fromfile(path_information + "_ComputedValuesFast.bin")
    eik_coords = np.genfromtxt(path_information + "_MeshPoints.txt", delimiter=",")
    eik_grads = np.fromfile(path_information + "_ComputedGradientsFast.bin")
    eik_grads = eik_grads.reshape(len(eik_coords), 2)
    triangles_points = np.genfromtxt(path_information + "_Faces.txt", delimiter=",")
    nPoints = len(eik_vals)
    # Depending on the information given. If no essential information is given we have to compute it
    if x0 is None:
        x0 = x0_default
    if eta1 is None:
        eta1 = eta1_default
    if eta2 is None:
        eta2 = eta2_default
    path_trueSolution = Path( path_information + "_true_values.txt" )
    path_trueGrads = Path( path_information + "_true_grads.txt" )
    path_trueSolGrid = Path( path_information + "_trueSolGrid.txt" )
    if path_trueSolution.is_file() and path_trueGrads.is_file():
        trueSol = np.genfromtxt(path_trueSolution, delimiter=",")
        trueGrads = np.genfromtxt(path_trueGrads, delimiter=",")
    else:
        if( type == 0):
            trueSol, trueGrads = solutionLineTwoIndices(eik_coords)
        elif( type == 1):
            trueSol, trueGrads = solutionLineOneIndex(eik_coords)
        np.savetxt(path_trueSolution, trueSol, delimiter =', ', fmt = '%.8f' )
        np.savetxt(path_trueGrads, trueGrads, delimiter =', ', fmt = '%.8f' )
    if path_trueSolGrid.is_file():
        true_solGrid = np.genfromtxt(path_trueSolGrid, delimiter=",")
    else:
        true_solGrid = np.zeros( xi.shape)
        if( type == 0):
            for i in range(ny_default):
                for j in range(nx_default):
                    p = np.array( [[xi[i,j], yi[i,j] ]] )
                    trueE, trueG = solutionLineTwoIndices(p)
                    true_solGrid[i,j] = trueE
        elif( type == 1):
            for i in range(ny_default):
                for j in range(nx_default):
                    p = np.array( [[xi[i,j], yi[i,j] ]] )
                    trueE, trueG = solutionLineOneIndex(p)
                    true_solGrid[i,j] = trueE
        np.savetxt(path_trueSolGrid, true_solGrid, delimiter =', ', fmt = '%.8f' )
    ### COMPUTE THE ERRORS
    # We need the following in order to plot the triangulation
    triang = tri.Triangulation(eik_coords[:,0], eik_coords[:, 1], triangles_points)
    interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
    zi_lin = interp_lin(xi, yi)
    errors_inter = true_solGrid - zi_lin
    errorsAbs_inter = abs( errors_inter )
    points_errors_grads = np.zeros((nPoints))
    points_errors_eik = np.zeros((nPoints))
    for i in range(nPoints):
        points_errors_grads[i] = angle_error( trueGrads[i, :], eik_grads[i, :])
        points_errors_eik[i] = trueSol[i] - eik_vals[i]

    indOutsideInitialBall = np.where( np.sqrt( np.power( eik_coords[:, 0] - x0[0], 2) + np.power( eik_coords[:, 1] - x0[1], 2)) > 0.5)
    vMaxAbs = np.amax(points_errors_eik[indOutsideInitialBall])
    vMaxAbsGrads = np.amax(points_errors_grads[indOutsideInitialBall])
        
    # ABSOLUTE ERRORS (POINT-WISE)
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_1 = plt.imshow( errorsAbs_inter, cmap = colormap2, extent = [-2, 2, -2, 2], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors")
    plt.colorbar(im2_1)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig( path_to_save + H + "/" + H + '_PointErrors.png', dpi=my_dpi * 10)

    # SIGNED POINT WISE ERRORS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_2 = plt.imshow( errors_inter, cmap = colormap3, extent = [-2, 2, -2, 2], origin='lower', vmin=-vMaxAbs, vmax=vMaxAbs  )
    plt.title("Signed point wise errors")
    plt.colorbar(im2_2)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig( path_to_save + H + "/" + H + '_SignPointErrors.png', dpi=my_dpi * 10)

    # ABSOLUTE ERRORS WITH TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#454545')
    im2_3 = plt.imshow( errorsAbs_inter, cmap = colormap2, extent = [-2, 2, -2, 2], origin='lower', vmin = 0, vmax = vMaxAbs  )
    plt.title("Point wise absolute errors and triangulation")
    plt.colorbar(im2_3)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_PointErrors_Mesh.png', dpi=my_dpi * 10)

    # SIGNED POINT WISE ERRORS WITH TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#000000')
    im2_4 = plt.imshow( errors_inter, cmap = colormap3,  extent = [-2, 2, -2, 2], origin='lower', vmin = -vMaxAbs, vmax = vMaxAbs)
    plt.title("Point wise absolute errors and triangulation")
    plt.colorbar(im2_4)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
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
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
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
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
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
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
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
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
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
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSets_Grad.png', dpi=my_dpi * 10)

    # LEVEL SETS OF EXACT SOLUTION  + EXACT GRADIENTS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_6 = plt.contourf(xi, yi, true_solGrid, cmap = colormap2, levels = 30 )
    plt.contour(xi, yi, true_solGrid, colors = ["white"], levels = 30)
    plt.quiver(eik_coords[:, 0], eik_coords[:, 1], trueGrads[:, 0], trueGrads[:, 1])
    plt.title("Level sets and exact gradients")
    plt.colorbar(im2_6)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LevelSets_GradExact.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_8 = plt.imshow( zi_lin, cmap = colormap2, extent = [-2, 2, -2, 2], origin='lower'  )
    plt.title("Linear interpolation computed solution")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    plt.colorbar(im2_8)
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LinearInt.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION + TRIANGULATION ON TOP
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_9 = plt.imshow( zi_lin, cmap = colormap2, extent = [-2, 2, -2, 2], origin='lower'  )
    plt.triplot(eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-.', lw=0.2, c='#6800ff')
    plt.title("Linear interpolation and triangulation computed solution")
    plt.colorbar(im2_9)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LinearInt_Mesh.png', dpi=my_dpi * 10)

    # SOLUTION (SOLVER) + LINEAR INTERPOLATION + COMPUTED GRADIETNS
    fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    im2_10 = plt.imshow( zi_lin, cmap = colormap2, extent = [-2, 2, -2, 2], origin='lower'  )
    plt.quiver(eik_coords[:, 0], eik_coords[:, 1], eik_grads[:, 0], eik_grads[:, 1])
    plt.title("Linear interpolation and computed gradients")
    plt.colorbar(im2_10)
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_LinearInt_Grad.png', dpi=my_dpi * 10)

    # ERRORS IN GRADIENTS
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    im2_11 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = points_errors_grads, cmap = colormap1, vmin = 0, vmax = vMaxAbsGrads)
    plt.colorbar(im2_11)
    plt.title("Angle error in gradients")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_GradAngleErrors.png', dpi=my_dpi * 10)

    # ERRORS IN GRADIENTS + TRIANGULATION ON TOP
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
    im2_12 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = points_errors_grads, cmap = colormap1, vmin = 0, vmax = vMaxAbsGrads)
    plt.colorbar(im2_12)
    plt.title("Angle error in gradients and triangulation")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_GradAngleErrorsTri.png', dpi=my_dpi * 10)

    # ERRORS IN EIKONAL
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    im2_11 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 2 + round(7500/len(eik_coords)), c = points_errors_eik, cmap = colormap1, vmin = 0, vmax = vMaxAbs)
    plt.colorbar(im2_11)
    plt.title("Point errors Eikonal")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_PointErrEik.png', dpi=my_dpi * 10)

    balEikGrads = abs(points_errors_eik)/np.max(abs(points_errors_eik)) + abs(points_errors_grads)/np.max(abs(points_errors_grads))
    balEikGrads = balEikGrads/2

    # ERRORS IN EIKONAL AND GRAD BALANCED
    fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi)
    plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#9affe8", lw = 0.3 )
    im2_11 = plt.scatter(eik_coords[:, 0], eik_coords[:, 1], s = 5 + round(7500/len(eik_coords)), c = balEikGrads, cmap = colormap1, vmin = 0, vmax = max(abs(balEikGrads)))
    plt.colorbar(im2_11)
    plt.title("Average absolute point wise error Eikonal and gradient")
    ax = plt.gca()
    ax.set_aspect('equal')
    ax.set_xlim([-2,2])
    ax.set_ylim([-2,2])
    if (saveFigures):
        plt.savefig(path_to_save + H + "/" + H + '_PointErrAverage.png', dpi=my_dpi * 10)


    if(show):
        plt.show()
        

