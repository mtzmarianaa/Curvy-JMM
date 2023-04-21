

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

mu1 = 0.01
lam2 = 0.95
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



# Compute nS steps

nS = 9

f_vals = np.empty((nS))
subGradNorms = np.empty((nS))
normChangeParams = np.empty((nS))

paramsk, paramsCrTopk, paramsStTopk, gradParamsk, gradCrTopk, gradStTopk = params1, paramsCrTop1, paramsStTop1, gradParams1, gradCrTop1, gradStTop1

for i in range(nS):
    paramskM1, paramsCrTopkM1, paramsStTopkM1, gradParamskM1, gradCrTopkM1, gradStTopkM1 = paramsk, paramsCrTopk, paramsStTopk, gradParamsk, gradCrTopk, gradStTopk
    paramsk, paramsCrTopk, paramsStTopk, gradParamsk, gradCrTopk, gradStTopk = forwardPassUpdate(paramskM1, gammas, theta_gamma,
                                                                                             x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                             listIndices, listxk, listB0k,
                                                                                             listBk, listBkBk1, indCrTop,
                                                                                             paramsCrTopkM1, indStTop, paramsStTopkM1,
                                                                                             listCurvingInwards)
    fk = fObj_generalized(paramsk, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                          listxk, listB0k, listBk, listBkBk1,
                          indCrTop = indCrTop, paramsCrTop = paramsCrTopk,
                          indStTop = indStTop, paramsStTop = paramsStTopk)
    itt.plotFann(x0, listB0k, listxk, listBk, params = paramsk, indCrTop = indCrTop, paramsCrTop = paramsCrTopk, listBkBk1 = listBkBk1)
    ax = plt.gca()
    ax.set_aspect("auto")
    plt.title("After " + str(i + 1) +" steps of projected coordinate subgradient descent")
    f_vals[i] = fk
    subGradNorms[i] = sqrt( norm(gradParamsk)**2 + norm(gradCrTopk)**2 + norm(gradStTopk)**2)
    normChangeParams[i] = sqrt( norm( paramskM1 - paramsk)**2 + norm( paramsCrTopkM1 - paramsCrTopk)**2 + norm(paramsStTopkM1 - paramsStTopk)**2 )

fig = plt.figure(figsize=(800/96, 800/96), dpi=96) 
plt.semilogy( range(0, nS), f_vals, c = "#394664", linewidth = 0.8)
plt.xlabel("Iteration")
plt.ylabel("Function value")
plt.title("Function value at each iteration")


fig = plt.figure(figsize=(800/96, 800/96), dpi=96) 
plt.semilogy( range(0, nS), subGradNorms, c = "#394664", linewidth = 0.8)
plt.xlabel("Iteration")
plt.ylabel("Norm of subgradients")
plt.title("Norm of subgradients at each iteration")


fig = plt.figure(figsize=(800/96, 800/96), dpi=96) 
plt.semilogy( range(0, nS), normChangeParams, c = "#394664", linewidth = 0.8)
plt.xlabel("Iteration")
plt.ylabel("Norm of change in parameters")
plt.title("Norm of change in parameters at each iteration")



rs, ss = np.meshgrid( np.linspace(0, 1, 200), np.linspace(0, 1, 200))
fk_grid = np.empty((200,200))

for i in range(200):
    for j in range(200):
        r1 = rs[i,j]
        s1 = ss[i,j]
        paramsCrTop = np.array([r1, s1])
        fk_grid[i,j] = fObj_generalized(paramsk, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                          listxk, listB0k, listBk, listBkBk1,
                          indCrTop = indCrTop, paramsCrTop = paramsCrTop,
                          indStTop = indStTop, paramsStTop = paramsStTopk)

fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
im = plt.imshow(fk_grid, cmap = colormap2, extent = [0,1,0,1], origin = "lower")
plt.scatter( paramsCrTop[0], paramsCrTop[1], c = "white", marker = "*", label = "optimum found")
#plt.contour(rs[0, :], ss[0, :], fk_grid, cmap = colormap2, extend = [0,1,0,1], origin = "lower", lower = 15)
plt.title("Level set of objective function")
plt.legend()




rs, ss = np.meshgrid( np.linspace(0, 1, 200), np.linspace(0, 1, 200))
fk_grid = np.empty((200,200))

for i in range(200):
    for j in range(200):
        r1 = rs[i,j]
        s1 = ss[i,j]
        paramsCrTop = np.array([r1, s1])
        fk_grid[i,j] = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                        listxk, listB0k, listBk, listBkBk1,
                                        indCrTop = indCrTop, paramsCrTop = paramsCrTop,
                                        indStTop = indStTop, paramsStTop = paramsStTop0)

fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
im = plt.imshow(fk_grid, cmap = colormap2, extent = [0,1,0,1], origin = "lower")
plt.scatter( paramsCrTop[0], paramsCrTop[1], c = "white", marker = "*", label = "optimum found")
plt.contour(rs[0, :], ss[:, 0], fk_grid, colors = ["white"], extend = [0,1,0,1], origin = "lower", levels = 15, linewidths = 0.7)
plt.title("Level set of objective function")
plt.colorbar(im)
plt.legend()


