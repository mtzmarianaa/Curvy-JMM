# SCRIPT TO VISUALIZE ERRORS (can I say this?)
from dataclasses import replace
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from numpy.linalg import norm 


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('magma')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

######################################################
## 1. Make the same kind of image plots as before, but compute the pointwise difference.

# Visualize my solver's output
x_FMM = np.fromfile('from_FMM_h001k8_bTree.bin')
x_FMM = x_FMM.reshape(513, 513)
colors = colormap1(x_FMM)
plt.figure(1)
plt.imshow(x_FMM, cmap=colormap2)
plt.colorbar(sm2)
plt.show(block=False)
plt.title("Result from FMM implementation")

# Visualize the point-wise errors in 2d
x_exact = np.fromfile('exact_h001k8.bin')
x_exact = x_exact.reshape(513, 513)
colors = colormap1(x_exact)
plt.figure(2)
plt.imshow(x_exact, cmap=colormap2)
plt.colorbar(sm2)
plt.show(block=False)
plt.title("Exact solution")


# Visuelize the point wise error between the two

error_points = abs(x_FMM - x_exact)
plt.figure(3)
plt.imshow(error_points, cmap=colormap2)
plt.colorbar(sm2)
plt.show(block=False)
plt.title("Absolute point wise errors")
# Seems suspicious

# Visualize in 3D?
plt.figure(4)
ax = plt.axes(projection='3d')
M, N = np.shape(error_points)
x_data = np.linspace(0, 0.2, 513)
y_data = np.linspace(0, 0.2, 513)
ax.contour3D(x_data, y_data, x_FMM - x_exact, 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Errors');
plt.title("3D point wise errors")

######################################################
## 2. Compute the relative error for varying problem size.
## there are 2*n + 1 nodes along each axis (so that h = 1/n), 
# then let n = 2^k for k = 1, 2, 3, 4, 5, 6, 7, 8.
## In this set up the starting point is always in the center of the grid

## Relative error with varying problem size (varying n), h = 0.01
# x_FMM_n = np.fromfile('from_FMM_h001k1.bin')
# x_exact_n = np.fromfile('exact_h001k1.bin')
# rel_err_n = [norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2)]
# x_FMM_n = np.fromfile('from_FMM_h001k2.bin')
# x_exact_n = np.fromfile('exact_h001k2.bin')
# rel_err_n.append(norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2))
# x_FMM_n = np.fromfile('from_FMM_h001k3.bin')
# x_exact_n = np.fromfile('exact_h001k3.bin')
# rel_err_n.append(norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2))
# x_FMM_n = np.fromfile('from_FMM_h001k4.bin')
# x_exact_n = np.fromfile('exact_h001k4.bin')
# rel_err_n.append(norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2))
# x_FMM_n = np.fromfile('from_FMM_h001k5.bin')
# x_exact_n = np.fromfile('exact_h001k5.bin')
# rel_err_n.append(norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2))
# x_FMM_n = np.fromfile('from_FMM_h001k6.bin')
# x_exact_n = np.fromfile('exact_h001k6.bin')
# rel_err_n.append(norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2))
# x_FMM_n = np.fromfile('from_FMM_h001k7.bin')
# x_exact_n = np.fromfile('exact_h001k7.bin')
# rel_err_n.append(norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2))
# x_FMM_n = np.fromfile('from_FMM_h001k8_bTree.bin')
# x_exact_n = np.fromfile('exact_h001k8_bTree.bin')
# rel_err_n.append(norm(x_FMM_n -x_exact_n, 2 )/norm(x_exact_n, 2))

# plt.figure(5)
# plt.plot([2*2**1+1 ,2*2**2 + 1,2*2**3 + 1,2*2**4 + 1,2*2**5+ 1,2*2**6 + 1,2*2**7 + 1,2*2**8 +1], rel_err_n, color = colormap2(130))
# plt.show(block=False)
# ax.set_xlabel('Relative error')
# ax.set_ylabel('k (from n=2^k, 2*n+1 number of points')
# plt.yscale('log')
# plt.xscale('log')
# plt.title("Relative error vs number of points")


# plt.figure(6)
# plt.plot([2**(-1), 2**(-2), 2**(-3), 2**(-4), 2**(-5), 2**(-6), 2**(-7), 2**(-8)], rel_err_n, color = colormap2(80))
# plt.show(block=False)
# ax.set_xlabel('Relative error')
# ax.set_ylabel('h')
# plt.yscale('log')
# plt.xscale('log')
# plt.title("Step size")


plt.show()

