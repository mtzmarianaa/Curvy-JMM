# SCRIPT TO VISUALIZE ERRORS (can I say this?)
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)

######################################################
## 1. Make the same kind of image plots as before, but compute the pointwise difference.

# Visualize my solver's output
x_FMM = np.fromfile('from_FMM.bin')
x_FMM = x_FMM.reshape(200, 200)
colors = colormap1(x_FMM)
plt.figure(1)
plt.imshow(x_FMM, cmap=colormap2)
plt.colorbar(sm2)
plt.show(block=False)

# Visualize the point-wise errors in 2d
x_exact = np.fromfile('exact.bin')
x_exact = x_exact.reshape(200, 200)
colors = colormap1(x_exact)
plt.figure(2)
plt.imshow(x_exact, cmap=colormap2)
plt.colorbar(sm2)
plt.show(block=False)


# Visuelize the point wise error between the two

error_points = abs(x_FMM - x_exact)
plt.figure(3)
plt.imshow(error_points, cmap=colormap2)
plt.colorbar(sm2)
plt.show(block=False)
# Seems suspicious

# Visualize in 3D?
plt.figure(4)
ax = plt.axes(projection='3d')
M, N = np.shape(error_points)
x_data = np.linspace(0, 0.2, 200)
y_data = np.linspace(0, 0.2, 200)
ax.contour3D(x_data, y_data, x_FMM - x_exact, 50, cmap=colormap2)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Errors');


######################################################
## 2. Compute the relative error for varying problem size.





plt.show()

