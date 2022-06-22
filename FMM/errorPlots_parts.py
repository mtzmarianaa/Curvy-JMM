# SCRIPT TO VISUALIZE ERRORS FROM VARIOUS PARTS OF FMM
from cProfile import label
from dataclasses import replace
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from numpy.linalg import norm 
from math import sqrt


colormap1 = plt.cm.get_cmap('cubehelix')
sm1 = plt.cm.ScalarMappable(cmap=colormap1)
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
colormap3 = plt.cm.get_cmap('magma')
sm3 = plt.cm.ScalarMappable(cmap=colormap3)

# ERRORS FROM THE OPTIMIZATION PART TO FIND THE OPTIMUM LAMBDA 

# How the approximated lambda goes to the actual value
lambdas_secant = np.fromfile('Outputs/convSec.bin')
iteration_n = range(1, len(lambdas_secant)+1)
exact_lambda = [sqrt(2)/2]*len(lambdas_secant)
plt.figure(1)

plt.plot(iteration_n, lambdas_secant, color = colormap2(100), label = 'From method', linewidth = 0.8)
plt.plot(iteration_n, exact_lambda, color = colormap2(50), label = 'Exact solution', linewidth = 1)
plt.legend()
plt.xlabel('Iteration number')
plt.ylabel('Lambda')
plt.show(block=False)
plt.title("Convergence of secant method")

# Log plot of the absolute error
plt.figure(2)
plt.plot(iteration_n, abs(lambdas_secant - exact_lambda), color = colormap2(75), linewidth = 0.8 )
plt.xlabel('Iteration number')
plt.ylabel('Absolute error')
plt.yscale('log')
plt.xscale('log')
plt.show(block = False)
plt.title('Absolute error of secant method')

# plot of the absolute error
plt.figure(3)
plt.plot(iteration_n, abs(lambdas_secant - exact_lambda), color = colormap2(75), linewidth = 0.8 )
plt.xlabel('Iteration number')
plt.ylabel('Absolute error')
plt.show(block = False)
plt.title('Absolute error of secant method')

plt.show()
