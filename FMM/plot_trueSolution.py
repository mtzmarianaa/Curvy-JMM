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

plt.ion()


colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)



saveFigures = False
nx = 36*1
ny = 42*1
my_dpi=96
eta1 = 1.0
eta2 = 1.452
x0 = np.array([-15, -10])
center = np.array([0,0])
R = 10
eps = np.finfo(np.float64).resolution

# Meshing
xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
true_solGrid = np.zeros(xi.shape)
type_solution = np.zeros(xi.shape)


# Computing the true solution and gradients

coords = np.zeros((nx*ny, 2))
true_grads = np.zeros((nx*ny, 2))
k = 0

for i in range(ny):
    for j in range(nx):
        sol, typeSol, trueGrad  = trueSolution(  xi[i, j], yi[i,j], x0, center, R, eta1, eta2  )
        true_solGrid[i, j] = sol
        type_solution[i, j] = typeSol
        coords[k, 0] = xi[i,j]
        coords[k, 1] = yi[i,j]
        true_grads[k, 0] = trueGrad[0]
        true_grads[k, 1] = trueGrad[1]
        k += 1

# Plot the gradients with the true solution

fig = plt.figure(figsize = (800/my_dpi, 800/my_dpi), dpi = my_dpi )
plt.axis('equal')
im_grad = plt.imshow(true_solGrid, cmap = colormap2, extent = [-18, 18, -18, 24], origin = "lower")
plt.quiver( coords[:, 0], coords[:, 1], true_grads[:, 0], true_grads[:,1]  )
plt.title("True gradients and true solution")
plt.show(block = False)
plt.colorbar(im_grad)

plt.show()


