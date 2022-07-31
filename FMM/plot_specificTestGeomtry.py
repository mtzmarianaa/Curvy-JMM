# SCRIPT TO VISUALIZE ERRORS (can I say this?)
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm 
from math import sqrt
import matplotlib.tri as tri
from scipy.optimize import NonlinearConstraint, minimize
import matplotlib.animation as animation

colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)

nx = 36*200
ny = 42*200
my_dpi=96


def rotate(angle):
    ax.view_init(azim=angle)

######################################################
######################################################
######################################################
#### H4
## 1. Plot of the output

eik_vals_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_ComputedValues.bin")
eik_coords_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_MeshPoints.txt", delimiter=",")
triangles_H4 = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_Faces.txt", delimiter=",")
eik_grads_H4 = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_ComputedGradients.bin");
eik_grads_H4 = eik_grads_H4.reshape(len(eik_coords_H4), 2)

# exact_values_H4 = []
# errors_H4 = []
# for i in range(len(eik_coords_H4)):
#     xi = eik_coords_H4[i, 0]
#     yi = eik_coords_H4[i, 1]
#     sol = exact_solution1(xi, yi)
#     exact_values_H4 += [sol]
#     errors_H4 += [ abs( sol - eik_vals_H4[i] ) ]



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
ax = plt.axes(projection='3d')
ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], eik_vals_H4, c= eik_vals_H4, cmap=colormap2)
plt.title("Computed eikonal values, test geometry H4")
plt.show(block = False)
rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
#rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_ComputedValues.gif', dpi=80, writer='imagemagick')





# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], exact_values_H4, c= exact_values_H4, cmap=colormap2)
# plt.title("Exact solution, test geometry H4")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# #rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_ExactSolution.gif', dpi=80, writer='imagemagick')


# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], errors_H4, c = errors_H4, cmap=colormap2)
# plt.title("Computed errors per point, test geometry H4")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# #rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_PointsPointErrors.gif', dpi=80, writer='imagemagick')



# We interpolate the solution on the triangles_H4 (so that we get a smooth plot + Sam´s idea)

# Since we are dealing with a square on [-5, 5], [-5, -5], [5, -5] , [5, 5]

xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4)
# To be able to use LinearTriInterpolator
interp_lin = tri.LinearTriInterpolator(triang, eik_vals_H4)
zi_lin = interp_lin(xi, -yi+6)

# # Contours of the errors_H4 in 3D and 2D
# solution_interpolated = np.zeros(zi_lin.shape)
# for i in range(len(xi)):
#     for j in range(len(yi)):
#         solution_interpolated[i, j] = exact_solution1(  xi[i, j], yi[i,j]  )
# errors_H4_abs = abs(zi_lin - solution_interpolated)

# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# im_bar44 = plt.imshow( solution_interpolated, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Exact solution, test geometry H4")
# plt.show(block = False)
# plt.colorbar(im_bar44)
# #plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_ExactSolution.png', dpi=my_dpi * 10)


# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# ax = plt.axes(projection='3d')
# ax.contour3D(xi, yi, errors_H4_abs , 50, cmap=colormap2)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('errors_H4');
# plt.title("3D point wise errors, test geometry H4")
# plt.show(block = False)
# rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
# #rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_PointErrors.gif', dpi=80, writer='imagemagick')


# Plot the absolute errors_H4 in 2D

# fig = plt.figure(46)
# plt.axis('equal')
# im_bar46 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Point wise absolute errors, test geometry H4")
# plt.show(block = False)
# plt.colorbar(im_bar46)
# #plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_PointErrors.png', dpi=my_dpi * 10)



# # The absolute errors_H4 in 2D with the triangulation

# fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
# plt.axis('equal')
# plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#04007e')
# im_bar47 = plt.imshow( errors_H4_abs, cmap = colormap2, extent=[-18,18,-18,24]  )
# plt.title("Point wise absolute errors and triangulation, test geometry H4")
# plt.show(block = False)
# plt.colorbar(im_bar47)
# #plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_PointErrors_Mesh.png', dpi=my_dpi * 10)



#Now we can plot + plot the triangulation + dots on top
# This plots the contours (I think it looks horrible)
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im_bar48 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2)
plt.scatter(eik_coords_H4[:, 0], eik_coords_H4[:, 1], c = eik_vals_H4, cmap = colormap2)
plt.triplot(eik_coords_H4[:, 0], eik_coords_H4[:, 1], triangles_H4, '-.', lw=0.5, c='#6800ff')
plt.title("Linear interpolation, test geometry H4")
plt.show(block = False)
plt.colorbar(im_bar48)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_LinearInt_Mesh.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im_bar49 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
plt.title("Linear interpolation, test geometry H4")
plt.show(block = False)
plt.colorbar(im_bar49)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_LinearInt.png', dpi=my_dpi * 10)



fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im_bar50 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
plt.quiver(eik_coords_H4[:, 0], eik_coords_H4[:, 1], eik_grads_H4[:, 0], eik_grads_H4[:, 1])
plt.title("Linear interpolation and computed eikonal gradient, test geometry H4")
plt.show(block = False)
plt.colorbar(im_bar50)
#plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/H4/H4_LinearInt_Grad.png', dpi=my_dpi * 10)



# averageH += [average_edge_length(eik_coords_H4, triangles_H4)]
# errorNorm += [norm( np.subtract(eik_vals_H4, exact_values_H4)  )/norm( exact_values_H4 )]
# nPointsH += [len(eik_coords_H4)]

plt.show()