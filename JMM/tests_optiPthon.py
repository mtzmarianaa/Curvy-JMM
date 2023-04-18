# Before going to C with the optimization ruting
# we test optiPython and intermediateTests here


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from optiPython import blockCoordinateGradient, plotResults, fObj_noTops, findRtan
from optiPython import gradient_TY, fObj_generalized, project_lamkGivenmuk1_noCr, project_mukGivenrk
from optiPython import project_mukGivenlamk1_noCr, project_rkGivenmuk, project_skGivenlamk1
from optiPython import project_lamkGivenskM1
from analyticSol_circle import trueSolution # So that we can test all of this in an actual "true geometry
import colorcet as cc
import matplotlib.colors as clr
from mpl_toolkits.mplot3d import axes3d
from scipy.optimize import root_scalar


colormap2  = clr.LinearSegmentedColormap.from_list('Retro',
                                                   [(0,    '#000000'),
                                                    (0.1, '#2c3454'),
                                                    (0.25, '#0033ff'),
                                                    (0.60, '#00f3ff'),
                                                    (1,    '#e800ff')], N=256)
maxIter = 30
tol = 1e-8
      
###########################################
## NAIVELY TEST OPTIPYTHON
## TESTS FOR THESE FUNCTIONS

# x0 = np.array([0.0,0.0])
# x1 = np.array([2, -0.2])
# x2 = np.array([1.5, 0.8])
# x3 = np.array([0.2, 1.2])
# xHat = np.array([-0.8, 0.7])
# B01 = np.array([2.2, 1])
# B01 = B01/norm(B01)
# B02 = np.array([1, 1.5])
# B02 = B02/norm(B02)
# B03 = np.array([0.2, 2])
# B03 = B03/norm(B03)
# B0Hat = np.array([-1, 0.4])
# B0Hat = B0Hat/norm(B0Hat)
# B1 = np.array([1, -0.6])
# B1 = B1/norm(B1)
# B2 = np.array([2, -0.2])
# B2 = B2/norm(B2)
# B3 = np.array([1, 2])
# B3 = B3/norm(B3)
# BHat = np.array([-1, 1])
# BHat = BHat/norm(BHat)


# # No tops = straight lines on the top of the triangles
# B1B2_0 = x2 - x1
# B1B2_1 = x2 - x1
# B2B3_0 = x3 - x2
# B2B3_1 = x3 - x2
# B3B4_0 = xHat - x3
# B3B4_1 = xHat - x3
# listBkBk1 = [B1B2_0, B1B2_1, B2B3_0, B2B3_1, B3B4_0, B3B4_1]



# mu1 = 0.5
# lam2 = 0.65
# mu2 = 0.5
# lam3 = 0.75
# mu3 = 0.5
# lam4 = 0.45
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4]

# xSource = np.array([ 1, -0.5 ])
# T0 = norm(xSource - x0)
# grad0 = (x0 - xSource)/T0
# T1 = norm(xSource - x1)
# grad1 = (x1 - xSource)/T1

# THat_true = norm(xHat - xSource)

# listIndices = [1.0, 1.0, 1.0, 1.0]
# listxk = [x0, x1, x2, x3, xHat]
# listB0k = [B01, B02, B03, B0Hat]
# listBk = [B1, B2, B3, BHat]


# # # # Another example

# print("Start test foward pass update \n\n")


# paramsOpt, listObjVals, listGrads, listChangefObj = blockCoordinateGradient(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol, plotSteps = False)

# plotResults(x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listB0k, listxk, listBk, params0, paramsOpt, listObjVals, listGrads, listChangefObj, trueSol = None)


# print("Value of objective function using old no tops function:\n")
# print(listObjVals[-1])

# fObj_gen = fObj_generalized(paramsOpt, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1 = listBkBk1,
#                             indCrTop = None, paramsCrTop = None, indStTop = None, paramsStTop = None)

# print("Value of objective function using new generalized function \n")
# print(fObj_gen)



# # mu1 = 0.5
# # lam2 = 0.65
# # mu2 = 0.5
# # lam3 = 0.75
# # mu3 = 0.5
# # lam4 = 0.45
# # params = [mu1, lam2, mu2, lam3, mu3, lam4]
# # listIndices = [3.5, 1, 0.5, 1.0]


# # paramsOpt, listObjVals, listGrads, listChangefObj = blockCoordinateGradient(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol, plotSteps = False)

# # plotResults(x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listB0k, listxk, listBk, params0, paramsOpt, listObjVals, listGrads, listChangefObj, trueSol = None)



###### Test the generalized triangle fan (with tops)
# Example for an unfeasible path in a generalized triangle fan

# x0 = np.array([0.0,0.0])
# x1 = np.array([2, -0.2])
# x2 = np.array([1.5, 0.8])
# x3 = np.array([0.2, 1.2])
# xHat = np.array([-0.8, 0.7])
# B01 = np.array([2.2, 1])
# B01 = B01/norm(B01)
# B02 = np.array([1, 1.5])
# B02 = B02/norm(B02)
# B03 = np.array([0.2, 2])
# B03 = B03/norm(B03)
# B0Hat = np.array([-1, 0.4])
# B0Hat = B0Hat/norm(B0Hat)
# B1 = np.array([1, -0.6])
# B1 = B1/norm(B1)
# B2 = np.array([2, -0.2])
# B2 = B2/norm(B2)
# B3 = np.array([1, 2])
# B3 = B3/norm(B3)
# BHat = np.array([-1, 1])
# BHat = BHat/norm(BHat)
# # For the triangle top
# B1B2_0 = np.array([-1,1])
# B1B2_0 = (B1B2_0/norm(B1B2_0))*0.3
# B1B2_1 = np.array([-1, 4])
# B1B2_1 = (B1B2_1/norm(B1B2_1))*0.7
# B2B3_0 = np.array([-3,1])
# B2B3_0 = (B2B3_0/norm(B2B3_0))*0.3
# B2B3_1 = np.array([-1, 4])
# B2B3_1 = (B2B3_1/norm(B2B3_1))*0.7
# B3B4_0 = np.array([-1,-1])
# B3B4_0 = (B3B4_0/norm(B3B4_0))*0.3
# B3B4_1 = np.array([-1.5, 1])
# B32B4_1 = (B3B4_1/norm(B3B4_1))*0.7


# listIndices = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# listxk = [x0, x1, x2, x3, xHat]
# listB0k = [B01, B02, B03, B0Hat]
# listBk = [B1, B2, B3, BHat]
# listBkBk1 = [B1B2_0, B1B2_1, B2B3_0, B2B3_1, B3B4_0, B3B4_1]


# xSource = np.array([ 1, -0.5 ])
# T0 = norm(xSource - x0)
# grad0 = (x0 - xSource)/T0
# T1 = norm(xSource - x1)
# grad1 = (x1 - xSource)/T1

# mu1 = 0.7
# lam2 = 0.8
# mu2 = 0.1
# lam3 = 0.1
# mu3 = 0.4
# lam4 = 0.5
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4, 1.0]

# indCrTop = [1]
# paramsCrTop = [0.3, 0.5]
# indStTop = [3]
# paramsStTop = [0.3, 0.7]

# itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop,
#          indStTop = indStTop, paramsStTop = paramsStTop)
# plt.title("Generalized triangle fan, creeping and going through tops")

# print("Value of objective function with these parameters: \n")
# fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1,
#                             indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = indStTop, paramsStTop = paramsStTop)

# print("    ", fObj_gen)


######## Another example with the generalized triangle fan

# Testing if the objective function which considers side edges can compute the perimeter of the
# triangle fan accordingly

# x0 = np.array([0.0,0.0])
# x1 = np.array([2, -0.2])
# x2 = np.array([1.5, 0.8])
# x3 = np.array([0.2, 1.2])
# xHat = np.array([-0.8, 0.7])
# B01 = np.array([2.2, 1])
# B01 = B01/norm(B01)
# B02 = np.array([1, 1.5])
# B02 = B02/norm(B02)
# B03 = np.array([0.2, 2])
# B03 = B03/norm(B03)
# B0Hat = np.array([-1, 0.4])
# B0Hat = B0Hat/norm(B0Hat)
# B1 = np.array([1, -0.6])
# B1 = B1/norm(B1)
# B2 = np.array([2, -0.2])
# B2 = B2/norm(B2)
# B3 = np.array([1, 2])
# B3 = B3/norm(B3)
# BHat = np.array([-1, 1])
# BHat = BHat/norm(BHat)
# # For the triangle top
# B1B2_0 = x2 - x1
# B1B2_1 = x2 - x1
# B2B3_0 = x3 - x2
# B2B3_1 = x3 - x2
# B3B4_0 = xHat - x3
# B3B4_1 = xHat - x3
# listBkBk1 = [B1B2_0, B1B2_1, B2B3_0, B2B3_1, B3B4_0, B3B4_1]




# listIndices = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
# listxk = [x0, x1, x2, x3, xHat]
# listB0k = [B01, B02, B03, B0Hat]
# listBk = [B1, B2, B3, BHat]

# mu1 = 0
# lam2 = 1
# mu2 = 1
# lam3 = 1
# mu3 = 1
# lam4 = 1
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4, 1.0]

# indCrTop = [1]
# paramsCrTop = [0, 1]
# indStTop = [3]
# paramsStTop = [0, 1]

# itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop,
#              indStTop = None, paramsStTop = None)
# plt.title("Generalized triangle fan, creeping and going through tops")

# print("Value of objective function with these parameters: \n")
# fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1,
#                             indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = indStTop, paramsStTop = paramsStTop)

# print("    ", fObj_gen)

# print("The value of this function should be just the permiter:")

# print(T0 + norm(x1-x0)+norm(x2-x1)+norm(x3-x2)+norm(xHat-x3))



##################################
#### Test the project_lamkGivenmuk1_noCr function on a single curvy triangle


x0 = np.array([0.0,0.0])
x1 = np.array([2, -0.2])
x2 = np.array([1.5, 0.8])
B01 = np.array([2.2, 1])
B01 = B01/norm(B01)
B02 = np.array([1, 1.5])
B02 = B02/norm(B02)
B1 = np.array([1, -0.6])
B1 = B1/norm(B1)
B2 = np.array([2, -0.2])
B2 = B2/norm(B2)
# For the triangle top
B1B2_0 = np.array([-1, 0.6])
B1B2_0 = (B1B2_0/norm(B1B2_0))*4
B1B2_1 = np.array([2, -0.2])
B1B2_1 = (B1B2_1/norm(B1B2_1))*0.5


listIndices = [1.0, 1.0, 1.0]
listxk = [x0, x1, x2]
listB0k = [B01, B02]
listBk = [B1, B2]
listBkBk1 = [B1B2_0, B1B2_1]


xSource = np.array([ 1, -0.5 ])
T0 = norm(xSource - x0)
grad0 = (x0 - xSource)/T0
T1 = norm(xSource - x1)
grad1 = (x1 - xSource)/T1


#############################################################
# Project lam2 using mu

mu1 = 0.7
lam2 = 0.8
params0 = [mu1, lam2, 1.0]


itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = None, paramsCrTop = None,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, no points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = None, paramsCrTop = None, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


lam2_proj = project_lamkGivenmuk1_noCr(mu1, lam2, x0, B01, x1, B1, B02, x2, B2, B1B2_0, B1B2_1)

paramsProj = [mu1, lam2_proj, 1.0]

# We also want to plot the a_tan used for this
z1 = itt.hermite_boundary(mu1, x0, B01, x1, B1)
rPass = lambda r: findRtan(r, x1, x2, B1B2_0, B1B2_1, z1)
rootTan = root_scalar(rPass, bracket = [0,1])
r_tan = rootTan.root
a_tan = itt.hermite_boundary(r_tan, x1, B1B2_0, x2, B1B2_1)

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = paramsProj,
         indCrTop = None, paramsCrTop = None,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, no points on side edge, lam2 projected back")
plt.scatter(a_tan[0], a_tan[1], marker = "+", c = "black", label = "tangent point")
plt.legend()

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(paramsProj, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = None, paramsCrTop = None, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)

#############################################################
# Project mu1 using lam2
mu1_proj = project_mukGivenlamk1_noCr(mu1, lam2, x0, B01, x1, B1, B02, x2, B2, B1B2_0, B1B2_1)
paramsProj2 = [mu1_proj, lam2, 1.0]

# We also want to plot a_tan used for this
yk1 = itt.hermite_boundary(lam2, x0, B02, x2, B2)
rPass = lambda r: findRtan(r, x1, x2, B1B2_0, B1B2_1, yk1)
rootTan = root_scalar(rPass, method = "secant", x0 = 0.4, x1 = 0.5)
r_tan = rootTan.root
a_tan = itt.hermite_boundary(r_tan, x1, B1B2_0, x2, B1B2_1)


itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = paramsProj2,
         indCrTop = None, paramsCrTop = None,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, no points on side edge, mu1 projected back")
plt.scatter(a_tan[0], a_tan[1], marker = "+", c = "black", label = "tangent point")
plt.legend()

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(paramsProj2, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = None, paramsCrTop = None, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)

#############################################################
# Project back rk given muk

# Test rMax

indCrTop = [1]
r1 = 0.94
s1 = 0.8
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


r1_proj = project_rkGivenmuk(r1, mu1, x0, B01, x1, B1, x2, B2, B1B2_0, B1B2_1)
paramsCrTop_proj = [r1_proj, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project r1 given mu1, rMax")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


# Test rMin

mu1 = 0
lam2 = 0.8
params0 = [mu1, lam2, 1.0]

indCrTop = [1]
r1 = 0
s1 = 0.8
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


r1_proj = project_rkGivenmuk(r1, mu1, x0, B01, x1, B1, x2, B2, B1B2_0, B1B2_1)
paramsCrTop_proj = [r1_proj, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project r1 given mu1, rMin")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)




#############################################################
# Project back sk given lamk1

# Test sMin

indCrTop = [1]
r1 = 0.1
s1 = 0.15
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


s1_proj = project_skGivenlamk1(s1, lam2, x0, B02, x2, B2, x1, B1B2_0, B1B2_1)
paramsCrTop_proj = [r1, s1_proj]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project s1 given lam2, sMin")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


# Test sMax

mu1 = 0.9
lam2 = 0
params0 = [mu1, lam2, 1.0]

indCrTop = [1]
r1 = 0.1
s1 = 0.95
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


s1_proj = project_skGivenlamk1(s1, lam2, x0, B02, x2, B2, x1, B1B2_0, B1B2_1)
paramsCrTop_proj = [r1, s1_proj]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project s1 given lam2, sMax")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop_proj, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


#############################################################
# Project back muk given rk


# Test muMin

mu1 = 0.1
lam2 = 0.8
params0 = [mu1, lam2, 1.0]

indCrTop = [1]
r1 = 0.1
s1 = 0.5
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)

mu1_proj = project_mukGivenrk(mu1, r1, x0, B01, x1, B1, B1B2_0, x2, B1B2_1)

params_proj = [mu1_proj, lam2, 1.0]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params_proj,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project mu1 given r1, muMin")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params_proj, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)



# Test muMax

mu1 = 0.9
lam2 = 0.95
params0 = [mu1, lam2, 1.0]

indCrTop = [1]
r1 = 0.7
s1 = 0.8 
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)

mu1_proj = project_mukGivenrk(mu1, r1, x0, B01, x1, B1, B1B2_0, x2, B1B2_1)

params_proj = [mu1_proj, lam2, 1.0]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params_proj,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project mu1 given r1, muMax")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params_proj, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


#############################################################
# Project back lamk given skM1

# Test lamMin

mu1 = 0.9
lam2 = 0
params0 = [mu1, lam2, 1.0]

indCrTop = [1]
r1 = 0.1
s1 = 0.9
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)

lam2_proj = project_lamkGivenskM1(lam2, s1, x0, B02, x2, B2, B1B2_0, x1, B1B2_1)

paramsProj = [mu1, lam2_proj, 1.0]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = paramsProj,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project lam2 fiven s1, lamMin")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(paramsProj, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)


# Test lamMax

mu1 = 0.9
lam2 = 0.8
params0 = [mu1, lam2, 1.0]

indCrTop = [1]
r1 = 0.1
s1 = 0.4
paramsCrTop = [r1, s1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, unfeasible")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)

lam2_proj = project_lamkGivenskM1(lam2, s1, x0, B02, x2, B2, B1B2_0, x1, B1B2_1)

paramsProj = [mu1, lam2_proj, 1.0]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = paramsProj,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, points on side edge, project lam2 fiven s1, lamMax")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(paramsProj, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = None, paramsStTop = None)

print("    ", fObj_gen)












