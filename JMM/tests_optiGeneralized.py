

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from optiPython import blockCoordinateGradient, plotResults, fObj_noTops, gradient_TY, fObj_generalized, forwardPassUpdate
from analyticSol_circle import trueSolution # So that we can test all of this in an actual "true geometry
import colorcet as cc
import matplotlib.colors as clr
from mpl_toolkits.mplot3d import axes3d


# colormap2  = clr.LinearSegmentedColormap.from_list('Retro',
#                                                    [(0,    '#000000'),
#                                                     (0.1, '#2c3454'),
#                                                     (0.25, '#0033ff'),
#                                                     (0.60, '#00f3ff'),
#                                                     (1,    '#e800ff')], N=256)

colormap2 = "cet_linear_worb_100_25_c53_r"
colormap2_r = "cet_linear_worb_100_25_c53"

maxIter = 30
tol = 1e-10


###########################################
###########################################
###########################################
###########################################
## TEST OPTIPYTHON IN THE CIRCLE 
## Recall that we know the analytic solution to the circle problem

# Test 1: x0 and x1 are on the boundary, xHat=x2 is inside the circle



###########################################
print("\n\n\n x0 and x1 on the boundary close\n\n")

xSource = np.array([-15, -10])
center = np.array([0,0])
R = 10.0
eta1 = 1.0
eta2 = 1.452

t0 = 1.865*pi
t1 = 1.87*pi

# Define the points
x1 = np.array([10*cos(t0), 10*sin(t0)])
x2 = np.array([10*cos(t1), 10*sin(t1)])
x0 = np.array([x1[0] , x1[1] - 0.075])
xHat = x2

hc = norm(x1-x2)

# Their derivatives
B1B2_0 = np.array([-sin(t0), cos(t0)])
B1B2_0 = (B1B2_0/norm(B1B2_0))*sqrt(hc)
B1B2_1 = np.array([-sin(t1), cos(t1)])
B1B2_1 = (B1B2_1/norm(B1B2_1))*sqrt(hc)
B01 = x1-x0
B1 = np.copy(B01)
B02 = x2-x0
B2 = np.copy(B02)

listIndices = [eta1, eta1, eta2]
listxk = [x0, x1, x2]
listB0k = [B01, B02]
listBk = [B1, B2]
listBkBk1 = [B1B2_0, B1B2_1]

# Compute solutions

T0, type0, grad0 = trueSolution(x0[0], x0[1], xSource, center, R, eta1, eta2)
grad0 = np.array(grad0)
T1, type1, grad1 = trueSolution(x1[0], x1[1], xSource, center, R, eta1, eta2)
grad1 = np.array(grad1)

# Compute the true solution for x2

T2, type2, grad2 = trueSolution(x2[0], x2[1], xSource, center, R, eta1, eta2)

# Use blockCoordinateGradient

mu1 = 0.4
lam2 = 0.6
params0 = [mu1, lam2, 1.0]
r1 = 0.1
s1 = 0.4
indCrTop = [1]
paramsCrTop0 = [r1, s1]
indStTop = None
paramsStTop0 = None
f0 = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                      listxk, listB0k, listBk, listBkBk1,
                      indCrTop = indCrTop, paramsCrTop = paramsCrTop0,
                      indStTop = indStTop, paramsStTop = paramsStTop0)

itt.plotFann(x0, listB0k, listxk, listBk, params = params0, indCrTop = indCrTop, paramsCrTop = paramsCrTop0, listBkBk1 = listBkBk1)
plt.title("Initial parameters")

# Compute one step of the pass forward method

listCurvingInwards = [1]
gammas = [0.1]
theta_gamma = 1

params1, paramsCrTop1, paramsStTop1, gradParams1, gradCrTop1, gradStTop1 = forwardPassUpdate(params0, gammas, theta_gamma,
                                                                                             x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                             listIndices, listxk, listB0k,
                                                                                             listBk, listBkBk1, indCrTop,
                                                                                             paramsCrTop0, indStTop, paramsStTop0,
                                                                                             listCurvingInwards)
# Plot to see what happened

f1 = fObj_generalized(params1, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                      listxk, listB0k, listBk, listBkBk1,
                      indCrTop = indCrTop, paramsCrTop = paramsCrTop1,
                      indStTop = indStTop, paramsStTop = paramsStTop1)
itt.plotFann(x0, listB0k, listxk, listBk, params = params1, indCrTop = indCrTop, paramsCrTop = paramsCrTop1, listBkBk1 = listBkBk1)
plt.title("After one step of projected coordinate subgradient descent")


