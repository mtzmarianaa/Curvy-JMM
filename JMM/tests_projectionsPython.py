## Testing all the different projections on the triangle fan
## from optiPython


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
import intermediateTests as itt
from optiPython import blockCoordinateGradient, plotResults, fObj_noTops, findRtan
from optiPython import gradient_TY, fObj_generalized, project_lamkGivenmuk1_noCr, project_mukGivenrk
from optiPython import project_mukGivenlamk1_noCr, project_rkGivenmuk, project_skGivenlamk1
from optiPython import project_lamkGivenskM1
import colorcet as cc
import matplotlib.colors as clr
from scipy.optimize import root_scalar


maxIter = 30
tol = 1e-8


######################################################################################################
######################################################################################################
################################## TYPE 4 TRIANGLE FAN ##################################
######################################################################################################
######################################################################################################


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

mu1 = 0.9
lam2 = 0
params0 = [mu1, lam2, 1.0]

indCrTop = [1]
r1 = 0.05
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

mu1 = 0.85
lam2 = 0.87
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



