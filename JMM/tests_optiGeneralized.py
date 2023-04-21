

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from optiPython import blockCoordinateGradient, plotResults, fObj_noTops, gradient_TY, fObj_generalized, forwardPassUpdate
import optiPython as oP
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
print("f0: ", f0)
itt.plotFann(x0, listB0k, listxk, listBk, params = params0, indCrTop = indCrTop, paramsCrTop = paramsCrTop0, listBkBk1 = listBkBk1)
ax = plt.gca()
ax.set_aspect("auto")
plt.title("Initial parameters")


# # TEST ALL THE PROJECTIONS THAT COULD BE DONE IN THIS SET UP


# # Test project r1 given mu1

# r1_proj = oP.project_rkGivenmuk(r1, mu1, x0, B01, x1, B1, x2, B2, B1B2_0, B1B2_1)
# paramsCrTop_proj = [r1_proj, s1]

# itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj,
#              indStTop = None, paramsStTop = None)
# ax = plt.gca()
# ax.set_aspect("auto")
# plt.title("Generalized triangle fan, points on side edge, project r1 given mu1, no projection necessary")


# # Test project lam2 given s1

# lam2_proj = oP.project_lamkGivenskM1(lam2, s1, x0, B02, x2, B2, B1B2_0, x1, B1B2_1)
# paramsProj = [mu1, lam2_proj, 1.0]

# itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = paramsProj,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop0,
#              indStTop = None, paramsStTop = None)
# ax = plt.gca()
# ax.set_aspect("auto")
# plt.title("Generalized triangle fan, points on side edge, project lam2 fiven s1, no projection needed")


# # Test project s1 given lam2

# s1_proj = oP.project_skGivenlamk1(s1, lam2, x0, B02, x2, B2, x1, B1B2_0, B1B2_1)
# paramsCrTop_proj = [r1, s1_proj]

# itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj,
#              indStTop = None, paramsStTop = None)
# ax = plt.gca()
# ax.set_aspect("auto")
# plt.title("Generalized triangle fan, points on side edge, project s1 given lam2, no projection needed")


# # Test project mu1 given r1

# mu1_proj = oP.project_mukGivenrk(mu1, r1, x0, B01, x1, B1, B1B2_0, x2, B1B2_1)
# params_proj = [mu1_proj, lam2, 1.0]

# itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params_proj,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop0,
#              indStTop = None, paramsStTop = None)
# ax = plt.gca()
# ax.set_aspect("auto")
# plt.title("Generalized triangle fan, points on side edge, project mu1 given r1, no projection needed")




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
print("f1: ", f1)
itt.plotFann(x0, listB0k, listxk, listBk, params = params1, indCrTop = indCrTop, paramsCrTop = paramsCrTop1, listBkBk1 = listBkBk1)
ax = plt.gca()
ax.set_aspect("auto")
plt.title("After one step of projected coordinate subgradient descent")



# Compute another step of the pass forward method

listCurvingInwards = [1]
gammas = [0.1]
theta_gamma = 1

params2, paramsCrTop2, paramsStTop2, gradParams2, gradCrTop2, gradStTop2 = forwardPassUpdate(params1, gammas, theta_gamma,
                                                                                             x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                             listIndices, listxk, listB0k,
                                                                                             listBk, listBkBk1, indCrTop,
                                                                                             paramsCrTop1, indStTop, paramsStTop1,
                                                                                             listCurvingInwards)
# Plot to see what happened

f2 = fObj_generalized(params2, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                      listxk, listB0k, listBk, listBkBk1,
                      indCrTop = indCrTop, paramsCrTop = paramsCrTop2,
                      indStTop = indStTop, paramsStTop = paramsStTop2)
print("f2: ", f2)
itt.plotFann(x0, listB0k, listxk, listBk, params = params2, indCrTop = indCrTop, paramsCrTop = paramsCrTop2, listBkBk1 = listBkBk1)
ax = plt.gca()
ax.set_aspect("auto")
plt.title("After two steps of projected coordinate subgradient descent")






# Compute another step of the pass forward method

listCurvingInwards = [1]
gammas = [0.1]
theta_gamma = 1

params3, paramsCrTop3, paramsStTop3, gradParams3, gradCrTop3, gradStTop3 = forwardPassUpdate(params2, gammas, theta_gamma,
                                                                                             x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                             listIndices, listxk, listB0k,
                                                                                             listBk, listBkBk1, indCrTop,
                                                                                             paramsCrTop2, indStTop, paramsStTop2,
                                                                                             listCurvingInwards)
# Plot to see what happened

f3 = fObj_generalized(params3, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                      listxk, listB0k, listBk, listBkBk1,
                      indCrTop = indCrTop, paramsCrTop = paramsCrTop3,
                      indStTop = indStTop, paramsStTop = paramsStTop3)
print("f3: ", f3)
itt.plotFann(x0, listB0k, listxk, listBk, params = params3, indCrTop = indCrTop, paramsCrTop = paramsCrTop3, listBkBk1 = listBkBk1)
ax = plt.gca()
ax.set_aspect("auto")
plt.title("After three steps of projected coordinate subgradient descent")




# Compute another step of the pass forward method

listCurvingInwards = [1]
gammas = [0.1]
theta_gamma = 1

params4, paramsCrTop4, paramsStTop4, gradParams4, gradCrTop4, gradStTop4 = forwardPassUpdate(params3, gammas, theta_gamma,
                                                                                             x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                             listIndices, listxk, listB0k,
                                                                                             listBk, listBkBk1, indCrTop,
                                                                                             paramsCrTop3, indStTop, paramsStTop3,
                                                                                             listCurvingInwards)
# Plot to see what happened

f4 = fObj_generalized(params4, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                      listxk, listB0k, listBk, listBkBk1,
                      indCrTop = indCrTop, paramsCrTop = paramsCrTop4,
                      indStTop = indStTop, paramsStTop = paramsStTop4)
print("f4: ", f4)
itt.plotFann(x0, listB0k, listxk, listBk, params = params4, indCrTop = indCrTop, paramsCrTop = paramsCrTop4, listBkBk1 = listBkBk1)
ax = plt.gca()
ax.set_aspect("auto")
plt.title("After 4 steps of projected coordinate subgradient descent")






# Compute another step of the pass forward method

listCurvingInwards = [1]
gammas = [0.1]
theta_gamma = 1

params5, paramsCrTop5, paramsStTop5, gradParams5, gradCrTop5, gradStTop5 = forwardPassUpdate(params4, gammas, theta_gamma,
                                                                                             x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                             listIndices, listxk, listB0k,
                                                                                             listBk, listBkBk1, indCrTop,
                                                                                             paramsCrTop4, indStTop, paramsStTop4,
                                                                                             listCurvingInwards)
# Plot to see what happened

f5 = fObj_generalized(params5, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                      listxk, listB0k, listBk, listBkBk1,
                      indCrTop = indCrTop, paramsCrTop = paramsCrTop5,
                      indStTop = indStTop, paramsStTop = paramsStTop5)
print("f5: ", f5)
itt.plotFann(x0, listB0k, listxk, listBk, params = params5, indCrTop = indCrTop, paramsCrTop = paramsCrTop5, listBkBk1 = listBkBk1)
ax = plt.gca()
ax.set_aspect("auto")
plt.title("After 5 steps of projected coordinate subgradient descent")




# Compute another step of the pass forward method

listCurvingInwards = [1]
gammas = [0.1]
theta_gamma = 1

params6, paramsCrTop6, paramsStTop6, gradParams6, gradCrTop6, gradStTop6 = forwardPassUpdate(params5, gammas, theta_gamma,
                                                                                             x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                             listIndices, listxk, listB0k,
                                                                                             listBk, listBkBk1, indCrTop,
                                                                                             paramsCrTop5, indStTop, paramsStTop5,
                                                                                             listCurvingInwards)
# Plot to see what happened

f6 = fObj_generalized(params6, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                      listxk, listB0k, listBk, listBkBk1,
                      indCrTop = indCrTop, paramsCrTop = paramsCrTop6,
                      indStTop = indStTop, paramsStTop = paramsStTop6)
print("f6: ", f6)
itt.plotFann(x0, listB0k, listxk, listBk, params = params6, indCrTop = indCrTop, paramsCrTop = paramsCrTop6, listBkBk1 = listBkBk1)
ax = plt.gca()
ax.set_aspect("auto")
plt.title("After 6 steps of projected coordinate subgradient descent")



