# Script to generate plots from the square with just a circle

# SCRIPT TO VISUALIZE ERRORS (can I say this?)
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm 
from math import sqrt
import matplotlib.tri as tri
from scipy.optimize import NonlinearConstraint, minimize
import matplotlib.animation as animation
import tabulate
from numpy import subtract as sb


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)

nx = 36*4
ny = 42*4
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
        s = 1.348
    else:
        s = 1
    return s
        

def whichRegion(xhat):
    '''
    This function determines in which of the 3 regions in the test geometry just base xhat is
    '''
    if (norm(xhat)< 10):
        reg = 3
    else:
        reg = 1
    return reg

def regA1Bool(xhat):
    '''
    This function determines if xhat can be accessed directly from x0
    '''
    ytan = (26*(7 + sqrt(14))*xhat[0])/(273 - 42*sqrt(13)) + (140*sqrt(13) + 130*sqrt(14))/(91 - 14*sqrt(13))
    if ( (xhat[1] >= ytan) | (xhat[1]<= -10) |  ((xhat[0]<=0) & (norm(xhat)>= 10) & ( xhat[1]<= 20/sqrt(13))  ) ):
        ans = True
    else:
        ans = False
    return ans

def consRayIntoCircle(z):
    '''
    Constrain such that the segment from x0 to p0 is outside the circle
    '''
    x0 = np.array([-15, -10])
    p0 = z[0:2]
    t0 = sb(x0, p0)
    return np.dot(t0, p0)

def consRayFromCircle(z, xhat):
    '''
    Constrain such that the segment from x0 to p1 is inside the circle
    '''
    p1 = z[2:4]
    t1 = sb(xhat, p1)
    return -1*np.dot(t1, p1)

def pointOnCircle(p):
    return 10 - norm(p)
        
def trueSolution(xi, yi):
    '''
    Analytic solution
    '''
    xhat = np.array([xi, yi])
    x0 = np.array([-15, -10])
    if(whichRegion(xhat)== 1): 
        if(regA1Bool(xhat)): # if xhat is directly accesible from x0 just via reg1
            tau = norm(sb(xhat, x0))
        else: # if xhat is in x0 but the ray has to go trough reg3 to get to xhat
            f_t1 = lambda z: norm(sb(x0, z[0:2])) + 1.348*norm(sb(z[0:2], z[2:4])) + norm(sb(z[2:4], xhat))
            cons1 = [{'type':'ineq', 'fun': lambda z: consRayIntoCircle(z) }, 
                     {'type':'ineq', 'fun': lambda z: consRayFromCircle(z, xhat) },
                     {'type':'eq', 'fun': lambda z: pointOnCircle(z[0:2]) },
                     {'type':'eq', 'fun': lambda z: pointOnCircle(z[2:4]) }]
            t1_opt = minimize( f_t1, [0, -10, 10, 0], constraints = cons1 )
            tau = f_t1(t1_opt.x)
    else: # x0 starts in reg1, goes into reg3 and that's it
        f_t2 = lambda z: norm(sb(x0, z[0:2])) + 1.348*norm(sb(z[0:2], xhat))
        cons2 = cons1 = [{'type':'ineq', 'fun': lambda z: consRayIntoCircle(z) }, # segment from x0 to pt0 goes into the circle
                         {'type':'eq', 'fun': lambda z: pointOnCircle(z) }] 
        t2_opt = minimize(f_t2, [0, -10], constraints = cons2)
        tau = f_t2(t2_opt.x)
    return tau


n = 0
averageH = []
errorNorm = []
nPointsH = []

# times = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/Times.bin")


# Compute the analytic solution in a grid

xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
true_solGrid = np.zeros(xi.shape)


for i in range(ny):
    for j in range(nx):
        true_solGrid[i, j] = trueSolution(  xi[i, j], yi[i,j]  )

# We plot the true solution

fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im1 = plt.imshow( true_solGrid, cmap = colormap2, extent=[-18,18,-18,24]  )
plt.title("Exact solution, test geometry just base")
plt.show(block = False)
plt.colorbar(im1)
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/ExactSolution.png', dpi=my_dpi * 10)


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_surface(xi, yi, true_solGrid, cmap=colormap2, linewidth=0, antialiased=False)
plt.show(block = False)





# ######################################################
# ######################################################
# ######################################################
# #### H1
# ## 1. Plot of the output

# eik_vals_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedValues.bin")
# eik_coords_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_MeshPoints.txt", delimiter=",")
# triangles_H1 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Faces.txt", delimiter=",")
# eik_grads_H1 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedGradients.bin");
# eik_grads_H1 = eik_grads_H1.reshape(len(eik_coords_H1), 2)

# # exact_values_H1 = []
# # errors_H1 = []
# # for i in range(len(eik_coords_H1)):
# #     xi = eik_coords_H1[i, 0]
# #     yi = eik_coords_H1[i, 1]
# #     sol = exact_solution1(xi, yi)
# #     exact_values_H1 += [sol]
# #     errors_H1 += [ abs( sol - eik_vals_H1[i] ) ]



# my_dpi=96
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_vals_H1, c= eik_vals_H1, cmap=colormap2)
# plt.title("Computed eikonal values, test geometry just base h1")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_ComputedValues.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], exact_values_H1, c= exact_values_H1, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base h1")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_ExactSolution.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], errors_H1, c = errors_H1, cmap=colormap2)
# # plt.title("Computed errors per point, test geometry just base h1")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointsPointErrors.gif', dpi=80, writer='Pillow')



# # We interpolate the solution on the triangles_H1 (so that we get a smooth plot + Sam´s idea)

# # Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]




# # We need a triangulation object thing
# triang = tri.Triangulation(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1)
# # To be able to use LinearTriInterpolator
# interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H1)
# zi_lin = interp_lin(xi, -yi+6)

# # # Contours of the errors_H1 in 3D and 2D
# # solution_interpolated = np.zeros(zi_lin.shape)
# # for i in range(len(xi)):
# #     for j in range(len(yi)):
# #         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# # errors_H1_abs = abs(zi_lin - solution_interpolated)

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im4 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Exact solution, test geometry just base h1")
# # plt.show(block = False)
# # plt.colorbar(im4)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_ExactSolution.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.contour3D(xi, yi, errors_H1_abs , 50, cmap=colormap2)
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# # ax.set_zlabel('errors_H1');
# # plt.title("3D point wise errors, test geometry just base h1")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors.gif', dpi=80, writer='Pillow')


# # Plot the absolute errors_H1 in 2D

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im6 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors, test geometry just base h1")
# # plt.show(block = False)
# # plt.colorbar(im6)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors.png', dpi=my_dpi * 10)



# # The absolute errors_H1 in 2D with the triangulation

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#04007e')
# # im7 = plt.imshow( errors_H1_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors and triangulation, test geometry just base h1")
# # plt.show(block = False)
# # plt.colorbar(im7)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_PointErrors_Mesh.png', dpi=my_dpi * 10)



# #Now we can plot + plot the triangulation + dots on top
# # This plots the contours (I think it looks horrible)
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im8 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
# plt.scatter(eik_coords_H1[:, 0], eik_coords_H1[:, 1], c = eik_vals_H1, cmap = colormap2)
# plt.triplot(eik_coords_H1[:, 0], eik_coords_H1[:, 1], triangles_H1, '-.', lw=0.5, c='#6800ff')
# plt.title("Linear interpolation, test geometry just base h1")
# plt.show(block = False)
# plt.colorbar(im8)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Mesh.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im9 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Linear interpolation, test geometry just base h1")
# plt.show(block = False)
# plt.colorbar(im9)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im10 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.quiver(eik_coords_H1[:, 0], eik_coords_H1[:, 1], eik_grads_H1[:, 0], eik_grads_H1[:, 1])
# plt.title("Linear interpolation and computed eikonal gradient, test geometry just base h1")
# plt.show(block = False)
# plt.colorbar(im10)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H1/H1_LinearInt_Grad.png', dpi=my_dpi * 10)



# # averageH += [average_edge_length(eik_coords_H1, triangles_H1)]
# # errorNorm += [norm( np.subtract(eik_vals_H1, exact_values_H1)  )/norm( exact_values_H1 )]
# # nPointsH += [len(eik_coords_H1)]


# ######################################################
# ######################################################
# ######################################################
# #### H2
# ## 1. Plot of the output

# eik_vals_H2 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedValues.bin")
# eik_coords_H2 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_MeshPoints.txt", delimiter=",")
# triangles_H2 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Faces.txt", delimiter=",")
# eik_grads_H2 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedGradients.bin");
# eik_grads_H2 = eik_grads_H2.reshape(len(eik_coords_H2), 2)

# # exact_values_H2 = []
# # errors_H2 = []
# # for i in range(len(eik_coords_H2)):
# #     xi = eik_coords_H2[i, 0]
# #     yi = eik_coords_H2[i, 1]
# #     sol = exact_solution1(xi, yi)
# #     exact_values_H2 += [sol]
# #     errors_H2 += [ abs( sol - eik_vals_H2[i] ) ]



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], eik_vals_H2, c= eik_vals_H2, cmap=colormap2)
# plt.title("Computed eikonal values, test geometry just base h2")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_ComputedValues.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], exact_values_H2, c= exact_values_H2, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base h2")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_ExactSolution.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], errors_H2, c = errors_H2, cmap=colormap2)
# # plt.title("Computed errors per point, test geometry just base h2")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_PointsPointErrors.gif', dpi=80, writer='Pillow')



# # We interpolate the solution on the triangles_H2 (so that we get a smooth plot + Sam´s idea)

# # Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]


# # We need a triangulation object thing
# triang = tri.Triangulation(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2)
# # To be able to use LinearTriInterpolator
# interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H2)
# zi_lin = interp_lin(xi, -yi+6)

# # # Contours of the errors_H2 in 3D and 2D
# # solution_interpolated = np.zeros(zi_lin.shape)
# # for i in range(len(xi)):
# #     for j in range(len(yi)):
# #         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# # errors_H2_abs = abs(zi_lin - solution_interpolated)

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar14 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Exact solution, test geometry just base h2")
# # plt.show(block = False)
# # plt.colorbar(im_bar14)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_ExactSolution.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.contour3D(xi, yi, errors_H2_abs , 50, cmap=colormap2)
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# # ax.set_zlabel('errors_H2');
# # plt.title("3D point wise errors, test geometry just base h2")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_PointErrors.gif', dpi=80, writer='Pillow')


# # # Plot the absolute errors_H2 in 2D

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar16 = plt.imshow( errors_H2_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors, test geometry just base h2")
# # plt.show(block = False)
# # plt.colorbar(im_bar16)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_PointErrors.png', dpi=my_dpi * 10)



# # # The absolute errors_H2 in 2D with the triangulation

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # plt.triplot(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2, '-.', lw=0.5, c='#04007e')
# # im_bar17 = plt.imshow( errors_H2_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors and triangulation, test geometry just base h2")
# # plt.show(block = False)
# # plt.colorbar(im_bar17)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_PointErrors_Mesh.png', dpi=my_dpi * 10)



# #Now we can plot + plot the triangulation + dots on top
# # This plots the contours (I think it looks horrible)
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar18 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
# plt.scatter(eik_coords_H2[:, 0], eik_coords_H2[:, 1], c = eik_vals_H2, cmap = colormap2)
# plt.triplot(eik_coords_H2[:, 0], eik_coords_H2[:, 1], triangles_H2, '-.', lw=0.5, c='#6800ff')
# plt.title("Linear interpolation, test geometry just base h2")
# plt.show(block = False)
# plt.colorbar(im_bar18)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_LinearInt_Mesh.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar19 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Linear interpolation, test geometry just base h2")
# plt.show(block = False)
# plt.colorbar(im_bar19)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_LinearInt.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar20 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.quiver(eik_coords_H2[:, 0], eik_coords_H2[:, 1], eik_grads_H2[:, 0], eik_grads_H2[:, 1])
# plt.title("Linear interpolation and computed eikonal gradient, test geometry just base h2")
# plt.show(block = False)
# plt.colorbar(im_bar20)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H2/H2_LinearInt_Grad.png', dpi=my_dpi * 10)



# # averageH += [average_edge_length(eik_coords_H2, triangles_H2)]
# # errorNorm += [norm( np.subtract(eik_vals_H2, exact_values_H2)  )/norm( exact_values_H2 )]
# # nPointsH += [len(eik_coords_H2)]


# ######################################################
# ######################################################
# ######################################################
# #### H3
# ## 1. Plot of the output

# eik_vals_H3 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_ComputedValues.bin")
# eik_coords_H3 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_MeshPoints.txt", delimiter=",")
# triangles_H3 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_Faces.txt", delimiter=",")
# eik_grads_H3 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_ComputedGradients.bin");
# eik_grads_H3 = eik_grads_H3.reshape(len(eik_coords_H3), 2)

# # exact_values_H3 = []
# # errors_H3 = []
# # for i in range(len(eik_coords_H3)):
# #     xi = eik_coords_H3[i, 0]
# #     yi = eik_coords_H3[i, 1]
# #     sol = exact_solution1(xi, yi)
# #     exact_values_H3 += [sol]
# #     errors_H3 += [ abs( sol - eik_vals_H3[i] ) ]



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], eik_vals_H3, c= eik_vals_H3, cmap=colormap2)
# plt.title("Computed eikonal values, test geometry just base h3")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_ComputedValues.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], exact_values_H3, c= exact_values_H3, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base h3")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_ExactSolution.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], errors_H3, c = errors_H3, cmap=colormap2)
# # plt.title("Computed errors per point, test geometry just base h3")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_PointsPointErrors.gif', dpi=80, writer='Pillow')



# # We interpolate the solution on the triangles_H3 (so that we get a smooth plot + Sam´s idea)

# # Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]


# # We need a triangulation object thing
# triang = tri.Triangulation(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3)
# # To be able to use LinearTriInterpolator
# interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H3)
# zi_lin = interp_lin(xi, -yi+6)

# # # Contours of the errors_H3 in 3D and 2D
# # solution_interpolated = np.zeros(zi_lin.shape)
# # for i in range(len(xi)):
# #     for j in range(len(yi)):
# #         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# # errors_H3_abs = abs(zi_lin - solution_interpolated)

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar24 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Exact solution, test geometry just base h3")
# # plt.show(block = False)
# # plt.colorbar(im_bar24)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_ExactSolution.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.contour3D(xi, yi, errors_H3_abs , 50, cmap=colormap2)
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# # ax.set_zlabel('errors_H3');
# # plt.title("3D point wise errors, test geometry just base h3")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_PointErrors.gif', dpi=80, writer='Pillow')


# # # Plot the absolute errors_H3 in 2D

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar26 = plt.imshow( errors_H3_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors, test geometry just base h3")
# # plt.show(block = False)
# # plt.colorbar(im_bar26)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_PointErrors.png', dpi=my_dpi * 10)



# # # The absolute errors_H3 in 2D with the triangulation

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # plt.triplot(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3, '-.', lw=0.5, c='#04007e')
# # im_bar27 = plt.imshow( errors_H3_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors and triangulation, test geometry just base h3")
# # plt.show(block = False)
# # plt.colorbar(im_bar27)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_PointErrors_Mesh.png', dpi=my_dpi * 10)



# #Now we can plot + plot the triangulation + dots on top
# # This plots the contours (I think it looks horrible)
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar28 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
# plt.scatter(eik_coords_H3[:, 0], eik_coords_H3[:, 1], c = eik_vals_H3, cmap = colormap2)
# plt.triplot(eik_coords_H3[:, 0], eik_coords_H3[:, 1], triangles_H3, '-.', lw=0.5, c='#6800ff')
# plt.title("Linear interpolation, test geometry just base h3")
# plt.show(block = False)
# plt.colorbar(im_bar28)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_LinearInt_Mesh.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar29 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Linear interpolation, test geometry just base h3")
# plt.show(block = False)
# plt.colorbar(im_bar29)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_LinearInt.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar30 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.quiver(eik_coords_H3[:, 0], eik_coords_H3[:, 1], eik_grads_H3[:, 0], eik_grads_H3[:, 1])
# plt.title("Linear interpolation and computed eikonal gradient, test geometry just base h3")
# plt.show(block = False)
# plt.colorbar(im_bar30)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H3/H3_LinearInt_Grad.png', dpi=my_dpi * 10)



# # averageH += [average_edge_length(eik_coords_H3, triangles_H3)]
# # errorNorm += [norm( np.subtract(eik_vals_H3, exact_values_H3)  )/norm( exact_values_H3 )]
# # nPointsH += [len(eik_coords_H3)]

# ######################################################
# ######################################################
# ######################################################
# #### H4
# ## 1. Plot of the output

# eik_vals_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_ComputedValues.bin")
# eik_coords_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_MeshPoints.txt", delimiter=",")
# triangles_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_Faces.txt", delimiter=",")
# eik_grads_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_ComputedGradients.bin");
# eik_grads_H4 = eik_grads_H4.reshape(len(eik_coords_H4), 2)

# # exact_values_H4 = []
# # errors_H4 = []
# # for i in range(len(eik_coords_H4)):
# #     xi = eik_coords_H4[i, 0]
# #     yi = eik_coords_H4[i, 1]
# #     sol = exact_solution1(xi, yi)
# #     exact_values_H4 += [sol]
# #     errors_H4 += [ abs( sol - eik_vals_H4[i] ) ]



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], eik_vals_H4, c= eik_vals_H4, cmap=colormap2)
# plt.title("Computed eikonal values, test geometry just base h4")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_ComputedValues.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], exact_values_H4, c= exact_values_H4, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base h4")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_ExactSolution.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], errors_H4, c = errors_H4, cmap=colormap2)
# # plt.title("Computed errors per point, test geometry just base h4")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_PointsPointErrors.gif', dpi=80, writer='Pillow')



# # We interpolate the solution on the triangles_H4 (so that we get a smooth plot + Sam´s idea)

# # Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]


# # We need a triangulation object thing
# triang = tri.Triangulation(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4)
# # To be able to use LinearTriInterpolator
# interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H4)
# zi_lin = interp_lin(xi, -yi+6)

# # # Contours of the errors_H4 in 3D and 2D
# # solution_interpolated = np.zeros(zi_lin.shape)
# # for i in range(len(xi)):
# #     for j in range(len(yi)):
# #         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# # errors_H4_abs = abs(zi_lin - solution_interpolated)

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar34 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Exact solution, test geometry just base h4")
# # plt.show(block = False)
# # plt.colorbar(im_bar34)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_ExactSolution.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.contour3D(xi, yi, errors_H4_abs , 50, cmap=colormap2)
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# # ax.set_zlabel('errors_H4');
# # plt.title("3D point wise errors, test geometry just base h4")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_PointErrors.gif', dpi=80, writer='Pillow')


# # # Plot the absolute errors_H4 in 2D

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar36 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors, test geometry just base h4")
# # plt.show(block = False)
# # plt.colorbar(im_bar36)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_PointErrors.png', dpi=my_dpi * 10)



# # # The absolute errors_H4 in 2D with the triangulation

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#04007e')
# # im_bar37 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors and triangulation, test geometry just base h4")
# # plt.show(block = False)
# # plt.colorbar(im_bar37)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_PointErrors_Mesh.png', dpi=my_dpi * 10)



# #Now we can plot + plot the triangulation + dots on top
# # This plots the contours (I think it looks horrible)
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar38 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
# plt.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], c = eik_vals_H4, cmap = colormap2)
# plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#6800ff')
# plt.title("Linear interpolation, test geometry just base h4")
# plt.show(block = False)
# plt.colorbar(im_bar38)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_LinearInt_Mesh.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar39 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Linear interpolation, test geometry just base h4")
# plt.show(block = False)
# plt.colorbar(im_bar39)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_LinearInt.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar40 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.quiver(eik_coords_H4[:, 0], eik_coords_H4[:, 1], eik_grads_H4[:, 0], eik_grads_H4[:, 1])
# plt.title("Linear interpolation and computed eikonal gradient, test geometry just base h4")
# plt.show(block = False)
# plt.colorbar(im_bar40)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H4/H4_LinearInt_Grad.png', dpi=my_dpi * 10)



# # averageH += [average_edge_length(eik_coords_H4, triangles_H4)]
# # errorNorm += [norm( np.subtract(eik_vals_H4, exact_values_H4)  )/norm( exact_values_H4 )]
# # nPointsH += [len(eik_coords_H4)]


# ######################################################
# ######################################################
# ######################################################
# #### H5
# ## 1. Plot of the output

# eik_vals_H5 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_ComputedValues.bin")
# eik_coords_H5 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_MeshPoints.txt", delimiter=",")
# triangles_H5 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_Faces.txt", delimiter=",")
# eik_grads_H5 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_ComputedGradients.bin");
# eik_grads_H5 = eik_grads_H5.reshape(len(eik_coords_H5), 2)

# # exact_values_H5 = []
# # errors_H5 = []
# # for i in range(len(eik_coords_H5)):
# #     xi = eik_coords_H5[i, 0]
# #     yi = eik_coords_H5[i, 1]
# #     sol = exact_solution1(xi, yi)
# #     exact_values_H5 += [sol]
# #     errors_H5 += [ abs( sol - eik_vals_H5[i] ) ]



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], eik_vals_H5, c= eik_vals_H5, cmap=colormap2)
# plt.title("Computed eikonal values, test geometry just base H5")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_ComputedValues.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], exact_values_H5, c= exact_values_H5, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base H5")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_ExactSolution.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], errors_H5, c = errors_H5, cmap=colormap2)
# # plt.title("Computed errors per point, test geometry just base H5")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_PointsPointErrors.gif', dpi=80, writer='Pillow')



# # We interpolate the solution on the triangles_H5 (so that we get a smooth plot + Sam´s idea)

# # Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]


# # We need a triangulation object thing
# triang = tri.Triangulation(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5)
# # To be able to use LinearTriInterpolator
# interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H5)
# zi_lin = interp_lin(xi, -yi+6)

# # # Contours of the errors_H5 in 3D and 2D
# # solution_interpolated = np.zeros(zi_lin.shape)
# # for i in range(len(xi)):
# #     for j in range(len(yi)):
# #         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# # errors_H5_abs = abs(zi_lin - solution_interpolated)

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar34 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Exact solution, test geometry just base H5")
# # plt.show(block = False)
# # plt.colorbar(im_bar34)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_ExactSolution.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.contour3D(xi, yi, errors_H5_abs , 50, cmap=colormap2)
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# # ax.set_zlabel('errors_H5');
# # plt.title("3D point wise errors, test geometry just base H5")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_PointErrors.gif', dpi=80, writer='Pillow')


# # # Plot the absolute errors_H5 in 2D

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar36 = plt.imshow( errors_H5_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors, test geometry just base H5")
# # plt.show(block = False)
# # plt.colorbar(im_bar36)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_PointErrors.png', dpi=my_dpi * 10)



# # # The absolute errors_H5 in 2D with the triangulation

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # plt.triplot(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5, '-.', lw=0.5, c='#04007e')
# # im_bar37 = plt.imshow( errors_H5_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors and triangulation, test geometry just base H5")
# # plt.show(block = False)
# # plt.colorbar(im_bar37)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_PointErrors_Mesh.png', dpi=my_dpi * 10)



# #Now we can plot + plot the triangulation + dots on top
# # This plots the contours (I think it looks horrible)
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar38 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
# plt.scatter(eik_coords_H5[:, 0], eik_coords_H5[:, 1], c = eik_vals_H5, cmap = colormap2)
# plt.triplot(eik_coords_H5[:, 0], eik_coords_H5[:, 1], triangles_H5, '-.', lw=0.5, c='#6800ff')
# plt.title("Linear interpolation, test geometry just base H5")
# plt.show(block = False)
# plt.colorbar(im_bar38)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_LinearInt_Mesh.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar39 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Linear interpolation, test geometry just base H5")
# plt.show(block = False)
# plt.colorbar(im_bar39)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_LinearInt.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar40 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.quiver(eik_coords_H5[:, 0], eik_coords_H5[:, 1], eik_grads_H5[:, 0], eik_grads_H5[:, 1])
# plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H5")
# plt.show(block = False)
# plt.colorbar(im_bar40)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H5/H5_LinearInt_Grad.png', dpi=my_dpi * 10)



# # averageH += [average_edge_length(eik_coords_H5, triangles_H5)]
# # errorNorm += [norm( np.subtract(eik_vals_H5, exact_values_H5)  )/norm( exact_values_H5 )]
# # nPointsH += [len(eik_coords_H5)]


# ######################################################
# ######################################################
# ######################################################
# #### H6
# ## 1. Plot of the output

# eik_vals_H6 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_ComputedValues.bin")
# eik_coords_H6 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_MeshPoints.txt", delimiter=",")
# triangles_H6 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_Faces.txt", delimiter=",")
# eik_grads_H6 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_ComputedGradients.bin");
# eik_grads_H6 = eik_grads_H6.reshape(len(eik_coords_H6), 2)

# # exact_values_H6 = []
# # errors_H6 = []
# # for i in range(len(eik_coords_H6)):
# #     xi = eik_coords_H6[i, 0]
# #     yi = eik_coords_H6[i, 1]
# #     sol = exact_solution1(xi, yi)
# #     exact_values_H6 += [sol]
# #     errors_H6 += [ abs( sol - eik_vals_H6[i] ) ]



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], eik_vals_H6, c= eik_vals_H6, cmap=colormap2)
# plt.title("Computed eikonal values, test geometry just base h6")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_ComputedValues.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], exact_values_H6, c= exact_values_H6, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base h6")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_ExactSolution.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], errors_H6, c = errors_H6, cmap=colormap2)
# # plt.title("Computed errors per point, test geometry just base h6")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_PointsPointErrors.gif', dpi=80, writer='Pillow')



# # We interpolate the solution on the triangles_H6 (so that we get a smooth plot + Sam´s idea)

# # Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]


# # We need a triangulation object thing
# triang = tri.Triangulation(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6)
# # To be able to use LinearTriInterpolator
# interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H6)
# zi_lin = interp_lin(xi, -yi+6)

# # # Contours of the errors_H6 in 3D and 2D
# # solution_interpolated = np.zeros(zi_lin.shape)
# # for i in range(len(xi)):
# #     for j in range(len(yi)):
# #         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# # errors_H6_abs = abs(zi_lin - solution_interpolated)

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar44 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Exact solution, test geometry just base h6")
# # plt.show(block = False)
# # plt.colorbar(im_bar44)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_ExactSolution.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.contour3D(xi, yi, errors_H6_abs , 50, cmap=colormap2)
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# # ax.set_zlabel('errors_H6');
# # plt.title("3D point wise errors, test geometry just base h6")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_PointErrors.gif', dpi=80, writer='Pillow')


# # Plot the absolute errors_H6 in 2D

# # fig = plt.figure(46)
# # plt.axis('equal')
# # im_bar46 = plt.imshow( errors_H6_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors, test geometry just base h6")
# # plt.show(block = False)
# # plt.colorbar(im_bar46)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_PointErrors.png', dpi=my_dpi * 10)



# # # The absolute errors_H6 in 2D with the triangulation

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # plt.triplot(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6, '-.', lw=0.5, c='#04007e')
# # im_bar47 = plt.imshow( errors_H6_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors and triangulation, test geometry just base h6")
# # plt.show(block = False)
# # plt.colorbar(im_bar47)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_PointErrors_Mesh.png', dpi=my_dpi * 10)



# #Now we can plot + plot the triangulation + dots on top
# # This plots the contours (I think it looks horrible)
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar48 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
# plt.scatter(eik_coords_H6[:, 0], eik_coords_H6[:, 1], c = eik_vals_H6, cmap = colormap2)
# plt.triplot(eik_coords_H6[:, 0], eik_coords_H6[:, 1], triangles_H6, '-.', lw=0.5, c='#6800ff')
# plt.title("Linear interpolation, test geometry just base h6")
# plt.show(block = False)
# plt.colorbar(im_bar48)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_LinearInt_Mesh.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar49 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Linear interpolation, test geometry just base h6")
# plt.show(block = False)
# plt.colorbar(im_bar49)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_LinearInt.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar50 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.quiver(eik_coords_H6[:, 0], eik_coords_H6[:, 1], eik_grads_H6[:, 0], eik_grads_H6[:, 1])
# plt.title("Linear interpolation and computed eikonal gradient, test geometry just base h6")
# plt.show(block = False)
# plt.colorbar(im_bar50)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H6/H6_LinearInt_Grad.png', dpi=my_dpi * 10)



# # averageH += [average_edge_length(eik_coords_H6, triangles_H6)]
# # errorNorm += [norm( np.subtract(eik_vals_H6, exact_values_H6)  )/norm( exact_values_H6 )]
# # nPointsH += [len(eik_coords_H6)]



# ######################################################
# ######################################################
# ######################################################
# #### H7
# ## 1. Plot of the output

# eik_vals_H7 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_ComputedValues.bin")
# eik_coords_H7 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_MeshPoints.txt", delimiter=",")
# triangles_H7 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_Faces.txt", delimiter=",")
# eik_grads_H7 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H7/H7_ComputedGradients.bin");
# eik_grads_H7 = eik_grads_H7.reshape(len(eik_coords_H7), 2)

# # exact_values_H7 = []
# # errors_H7 = []
# # for i in range(len(eik_coords_H7)):
# #     xi = eik_coords_H7[i, 0]
# #     yi = eik_coords_H7[i, 1]
# #     sol = exact_solution1(xi, yi)
# #     exact_values_H7 += [sol]
# #     errors_H7 += [ abs( sol - eik_vals_H7[i] ) ]



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H7[:, 0], eik_coords_H7[:, 1], eik_vals_H7, c= eik_vals_H7, cmap=colormap2)
# plt.title("Computed eikonal values, test geometry just base H7")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_ComputedValues.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H7[:, 0], eik_coords_H7[:, 1], exact_values_H7, c= exact_values_H7, cmap=colormap2)
# # plt.title("Exact solution, test geometry just base H7")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_ExactSolution.gif', dpi=80, writer='Pillow')


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.scatter(eik_coords_H7[:, 0], eik_coords_H7[:, 1], errors_H7, c = errors_H7, cmap=colormap2)
# # plt.title("Computed errors per point, test geometry just base H7")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_PointsPointErrors.gif', dpi=80, writer='Pillow')



# # We interpolate the solution on the triangles_H7 (so that we get a smooth plot + Sam´s idea)

# # Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]


# # We need a triangulation object thing
# triang = tri.Triangulation(eik_coords_H7[:, 0], eik_coords_H7[:, 1], triangles_H7)
# # To be able to use LinearTriInterpolator
# interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H7)
# zi_lin = interp_lin(xi, -yi+6)

# # # Contours of the errors_H7 in 3D and 2D
# # solution_interpolated = np.zeros(zi_lin.shape)
# # for i in range(len(xi)):
# #     for j in range(len(yi)):
# #         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# # errors_H7_abs = abs(zi_lin - solution_interpolated)

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar34 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Exact solution, test geometry just base H7")
# # plt.show(block = False)
# # plt.colorbar(im_bar34)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_ExactSolution.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # ax = plt.axes(projection='3d')
# # ax.contour3D(xi, yi, errors_H7_abs , 50, cmap=colormap2)
# # ax.set_xlabel('x')
# # ax.set_ylabel('y')
# # ax.set_zlabel('errors_H7');
# # plt.title("3D point wise errors, test geometry just base H7")
# # plt.show(block = False)
# # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_PointErrors.gif', dpi=80, writer='Pillow')


# # # Plot the absolute errors_H7 in 2D

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # im_bar36 = plt.imshow( errors_H7_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors, test geometry just base H7")
# # plt.show(block = False)
# # plt.colorbar(im_bar36)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_PointErrors.png', dpi=my_dpi * 10)



# # # The absolute errors_H7 in 2D with the triangulation

# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.axis('equal')
# # plt.triplot(eik_coords_H7[:, 0], eik_coords_H7[:, 1], triangles_H7, '-.', lw=0.5, c='#04007e')
# # im_bar37 = plt.imshow( errors_H7_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# # plt.title("Point wise absolute errors and triangulation, test geometry just base H7")
# # plt.show(block = False)
# # plt.colorbar(im_bar37)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_PointErrors_Mesh.png', dpi=my_dpi * 10)



# #Now we can plot + plot the triangulation + dots on top
# # This plots the contours (I think it looks horrible)
# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar38 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
# plt.scatter(eik_coords_H7[:, 0], eik_coords_H7[:, 1], c = eik_vals_H7, cmap = colormap2)
# plt.triplot(eik_coords_H7[:, 0], eik_coords_H7[:, 1], triangles_H7, '-.', lw=0.5, c='#6800ff')
# plt.title("Linear interpolation, test geometry just base H7")
# plt.show(block = False)
# plt.colorbar(im_bar38)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_LinearInt_Mesh.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar39 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Linear interpolation, test geometry just base H7")
# plt.show(block = False)
# plt.colorbar(im_bar39)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_LinearInt.png', dpi=my_dpi * 10)



# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar40 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.quiver(eik_coords_H7[:, 0], eik_coords_H7[:, 1], eik_grads_H7[:, 0], eik_grads_H7[:, 1])
# plt.title("Linear interpolation and computed eikonal gradient, test geometry just base H7")
# plt.show(block = False)
# plt.colorbar(im_bar40)
# plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/H7/H7_LinearInt_Grad.png', dpi=my_dpi * 10)



# # averageH += [average_edge_length(eik_coords_H7, triangles_H7)]
# # errorNorm += [norm( np.subtract(eik_vals_H7, exact_values_H7)  )/norm( exact_values_H7 )]
# # nPointsH += [len(eik_coords_H7)]


# # ######################################################
# # ######################################################
# # ######################################################
# # ######################################################
# # ################## ERRORS ############################
# # ################### EACH #############################
# # ####################  H  #############################
# # ######################################################
# # ######################################################


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.loglog(averageH, errorNorm, c = '#6800ff')
# # plt.title("l2 errors and average edge length")
# # plt.xlabel("Average edge length")
# # plt.ylabel("Error")
# # plt.show(block = False)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_EdgeLength.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.loglog(averageH, nPointsH, c = '#6800ff')
# # plt.title("l2 errors and number of points in triangulation")
# # plt.xlabel("Number of points in triangulation")
# # plt.ylabel("Error")
# # plt.show(block = False)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Errors_nPoints.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.loglog(averageH, times, c = '#6800ff')
# # plt.title("Average edge length and time taken to solve")
# # plt.ylabel("Time taken to solve (sec)")
# # plt.xlabel("Average edge length")
# # plt.show(block = False)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/EdgeLength_Times.png', dpi=my_dpi * 10)


# # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# # plt.loglog(times, errorNorm, c = '#6800ff')
# # plt.title("Time taken to solve and l2 errors")
# # plt.xlabel("Time taken to solve (sec)")
# # plt.ylabel("Error")
# # plt.show(block = False)
# # plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/Times_Errors.png', dpi=my_dpi * 10)


# # table_sqTr = {"Average h": averageH, "Time taken": times, "l2 errors": errorNorm, "Points in triangulation": nPointsH}

# # print(tabulate(table_sqTr, headers="keys", tablefmt="latex"))

plt.show()