################ TESTS FOR THE GENERALIZED OPTIMIZATION METHOD BUT JUST USING TWO CURVY TRIANGLES
################ IN THE TRIANGLE FAN

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from optiPython import blockCoordinateGradient, plotResults, fObj_noTops, gradient_TY, fObj_generalized, forwardPassUpdate, blockCoordinateGradient_generalized
import optiPython as oP
from analyticSol_circle import trueSolution # So that we can test all of this in an actual "true geometry
import colorcet as cc
import matplotlib.colors as clr
from mpl_toolkits.mplot3d import axes3d


colormap2 = "cet_linear_worb_100_25_c53_r"
colormap2_r = "cet_linear_worb_100_25_c53"

maxIter = 30
tol = 1e-10



# #################################################################################################################################
# #################################################################################################################################
# #################################################################################################################################
# #################################################################################################################################
# ## JUST CR ON X1X2


# ###########################################

######
### x1x2 are on the boundary of the circle, they are the only points on the boundary
### the index of refraction outside the circle is 1 (i.e. for the internal triangles)
## the index of refraction inside the circle is 1.452 (i.e. for the triangle on top of x1x2)


print("\n\n\n x0 and x1 on the boundary close, x1 to the left of x0\n\n")

xSource = np.array([-15, -10])
center = np.array([0,0])
R = 10.0
eta1 = 1.0
eta2 = 1.452

t0 = 1.85*pi
t1 = 1.875*pi

# Define the points
x1 = np.array([10*cos(t0), 10*sin(t0)])
x2 = np.array([10*cos(t1), 10*sin(t1)])
x3 = np.array([x2[0] + 0.3, x2[1] - 0.15])
x0 = np.array([x1[0] + 0.1 , x1[1] - 0.2])
xHat = x3

hc = norm(x1-x2)

# Their derivatives
B1B2_0 = np.array([-sin(t0), cos(t0)])
B1B2_0 = (B1B2_0/norm(B1B2_0))*sqrt(hc)
B1B2_1 = np.array([-sin(t1), cos(t1)])
B1B2_1 = (B1B2_1/norm(B1B2_1))*sqrt(hc)
B2B3_0 = x3 - x2
B2B3_1 = x3 - x2
B01 = x1-x0
B1 = np.copy(B01)
B02 = x2-x0
B2 = np.copy(B02)
B03 = x3 - x0
B3 = np.copy(B03)

listIndices = [eta1, eta1, eta1, eta2, eta1]
listxk = [x0, x1, x2, x3]
listB0k = [B01, B02, B03]
listBk = [B1, B2, B3]
listBkBk1 = [B1B2_0, B1B2_1, B2B3_0, B2B3_1]
listCurvingInwards = [1,0]

# Compute solutions

T0, type0, grad0 = trueSolution(x0[0], x0[1], xSource, center, R, eta1, eta2)
grad0 = np.array(grad0)
T1, type1, grad1 = trueSolution(x1[0], x1[1], xSource, center, R, eta1, eta2)
grad1 = np.array(grad1)

# Compute the true solution for x3

T3, type3, grad3 = trueSolution(x3[0], x3[1], xSource, center, R, eta1, eta2)

# Use blockCoordinateGradient

mu1 = 0.4
lam2 = 0.5
mu2 = 0.1
lam3 = 0.7
params0 = [mu1, lam2, mu2, lam3, 1.0]
r1 = 0.54432577
s1 = 0.63893313
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
ax.set_aspect("equal")
plt.title("Initial parameters")


# Run opti

paramsk, paramsCrTopk, paramsStTopk, gradParamsk, gradCrTopk, gradStTopk, listObjVals,listGradNorms, listChangefObj, listChangeParams = blockCoordinateGradient_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1, indCrTop, paramsCrTop0, indStTop, paramsStTop0, listCurvingInwards, plotSteps = False, maxIter = 10)



# Plot

oP.plotResults(x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listB0k, listxk,
                listBk, params0, paramsk, listObjVals, listGradNorms, listChangefObj,
                listChangeParams = listChangeParams, trueSol = None, contours = True,
                listBkBk1 = listBkBk1, indCrTop = indCrTop, paramsCrTop0 = paramsCrTop0, 
                indStTop = indStTop, paramsStTop0 = paramsStTop0, 
                paramsCrTop = paramsCrTopk, paramsStTop = paramsStTopk)







######
### x1x2 are on the boundary of the circle, they are the only points on the boundary
### the index of refraction outside the circle is 1.452 (i.e. for the internal triangles)
## the index of refraction inside the circle is 1 (i.e. for the triangle on top of x1x2)


print("\n\n\n x0 and x1 on the boundary close, x1 to the left of x0\n\n")

xSource = np.array([-15, -10])
center = np.array([0,0])
R = 10.0
eta1 = 1.452
eta2 = 1.0

t0 = 1.85*pi
t1 = 1.875*pi

# Define the points
x1 = np.array([10*cos(t0), 10*sin(t0)])
x2 = np.array([10*cos(t1), 10*sin(t1)])
x3 = np.array([x2[0] + 0.3, x2[1] - 0.15])
x0 = np.array([x1[0] + 0.1 , x1[1] - 0.2])
xHat = x3

hc = norm(x1-x2)

# Their derivatives
B1B2_0 = np.array([-sin(t0), cos(t0)])
B1B2_0 = (B1B2_0/norm(B1B2_0))*sqrt(hc)
B1B2_1 = np.array([-sin(t1), cos(t1)])
B1B2_1 = (B1B2_1/norm(B1B2_1))*sqrt(hc)
B2B3_0 = x3 - x2
B2B3_1 = x3 - x2
B01 = x1-x0
B1 = np.copy(B01)
B02 = x2-x0
B2 = np.copy(B02)
B03 = x3 - x0
B3 = np.copy(B03)

listIndices = [eta1, eta1, eta1, eta2, eta1]
listxk = [x0, x1, x2, x3]
listB0k = [B01, B02, B03]
listBk = [B1, B2, B3]
listBkBk1 = [B1B2_0, B1B2_1, B2B3_0, B2B3_1]
listCurvingInwards = [1,0]

# Compute solutions

T0, type0, grad0 = trueSolution(x0[0], x0[1], xSource, center, R, eta1, eta2)
grad0 = np.array(grad0)
T1, type1, grad1 = trueSolution(x1[0], x1[1], xSource, center, R, eta1, eta2)
grad1 = np.array(grad1)

# Use blockCoordinateGradient

mu1 = 0.4
lam2 = 0.5
mu2 = 0.1
lam3 = 0.7
params0 = [mu1, lam2, mu2, lam3, 1.0]
r1 = 0.54432577
s1 = 0.63893313
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
ax.set_aspect("equal")
plt.title("Initial parameters")


# Run opti

paramsk, paramsCrTopk, paramsStTopk, gradParamsk, gradCrTopk, gradStTopk, listObjVals,listGradNorms, listChangefObj, listChangeParams = blockCoordinateGradient_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1, indCrTop, paramsCrTop0, indStTop, paramsStTop0, listCurvingInwards, plotSteps = False, maxIter = 25)



# Plot

oP.plotResults(x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listB0k, listxk,
                listBk, params0, paramsk, listObjVals, listGradNorms, listChangefObj,
                listChangeParams = listChangeParams, trueSol = None, contours = True,
                listBkBk1 = listBkBk1, indCrTop = indCrTop, paramsCrTop0 = paramsCrTop0, 
                indStTop = indStTop, paramsStTop0 = paramsStTop0, 
                paramsCrTop = paramsCrTopk, paramsStTop = paramsStTopk)









