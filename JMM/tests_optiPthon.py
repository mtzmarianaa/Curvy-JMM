# Before going to C with the optimization ruting
# we test optiPython and intermediateTests here


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from optiPython import blockCoordinateGradient, plotResults, fObj_noTops, findRtan
from optiPython import gradient_TY, fObj_genera
         indStTop = indStTop, paramsStTop = paramsStTop)
plt.title("Generalized triangle fan, creeping and going through tops")

print("Value of objective function with these parameters: \tesn")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = indStTop, paramsStTop = paramsStTop)

print("    ", fObj_gen)


##Another example with the generalized triangle fan

##Testing if the objective function which considers side edges can compute the perimeter of the
##triangle fan accordingly

x0 = np.array([0.0,0.0])
x1 = np.array([2, -0.2])
x2 = np.array([1.5, 0.8])
x3 = np.array([0.2, 1.2])
xHat = np.array([-0.8, 0.7])
B01 = np.array([2.2, 1])
B01 = B01/norm(B01)
B02 = np.array([1, 1.5])
B02 = B02/norm(B02)
B03 = np.array([0.2, 2])
B03 = B03/norm(B03)
B0Hat = np.array([-1, 0.4])
B0Hat = B0Hat/norm(B0Hat)
B1 = np.array([1, -0.6])
B1 = B1/norm(B1)
B2 = np.array([2, -0.2])
B2 = B2/norm(B2)
B3 = np.array([1, 2])
B3 = B3/norm(B3)
BHat = np.array([-1, 1])
BHat = BHat/norm(BHat)
# For the triangle top
B1B2_0 = x2 - x1
B1B2_1 = x2 - x1
B2B3_0 = x3 - x2
B2B3_1 = x3 - x2
B3B4_0 = xHat - x3
B3B4_1 = xHat - x3
listBkBk1 = [B1B2_0, B1B2_1, B2B3_0, B2B3_1, B3B4_0, B3B4_1]




listIndices = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
listxk = [x0, x1, x2, x3, xHat]
listB0k = [B01, B02, B03, B0Hat]
listBk = [B1, B2, B3, BHat]

mu1 = 0
lam2 = 1
mu2 = 1
lam3 = 1
mu3 = 1
lam4 = 1
params0 = [mu1, lam2, mu2, lam3, mu3, lam4, 1.0]

indCrTop = [1]
paramsCrTop = [0, 1]
indStTop = [3]
paramsStTop = [0, 1]

itt.plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
         indCrTop = indCrTop, paramsCrTop = paramsCrTop,
             indStTop = None, paramsStTop = None)
plt.title("Generalized triangle fan, creeping and going through tops")

print("Value of objective function with these parameters: \n")
fObj_gen = fObj_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1,
                            indCrTop = indCrTop, paramsCrTop = paramsCrTop, indStTop = indStTop, paramsStTop = paramsStTop)

print("    ", fObj_gen)

print("The value of this function should be just the permiter:")

print(T0 + norm(x1-x0)+norm(x2-x1)+norm(x3-x2)+norm(xHat-x3))





