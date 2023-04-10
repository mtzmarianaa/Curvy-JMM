# Before going to C with the optimization ruting
# we test optiPython and intermediateTests here


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from optiPython import blockCoordinateGradient, plotResults, fObj_noTops, gradient_TY
from analyticSol_circle import trueSolution # So that we can test all of this in an actual "true geometry
import colorcet as cc
import matplotlib.colors as clr
from mpl_toolkits.mplot3d import axes3d


colormap2  = clr.LinearSegmentedColormap.from_list('Retro',
                                                   [(0,    '#000000'),
                                                    (0.1, '#2c3454'),
                                                    (0.25, '#0033ff'),
                                                    (0.60, '#00f3ff'),
                                                    (1,    '#e800ff')], N=256)

      
###########################################
## NAIVELY TEST OPTIPYTHON
## TESTS FOR THESE FUNCTIONS

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


mu1 = 0.15
lam2 = 0.15
mu2 = 0.13
lam3 = 0.17
mu3 = 0.15
lam4 = 1
params = [mu1, lam2, mu2, lam3, mu3, lam4]

xSource = np.array([ 1, -0.5 ])
T0 = norm(xSource - x0)
grad0 = (x0 - xSource)/T0
T1 = norm(xSource - x1)
grad1 = (x1 - xSource)/T1

THat_true = norm(xHat - xSource)

listIndices = [1.0, 1.0, 1.0, 1.0]
listxk = [x0, x1, x2, x3, xHat]
listB0k = [B01, B02, B03, B0Hat]
listBk = [B1, B2, B3, BHat]


# # Compute the projected gradient descent

maxIter = 20
tol = 1e-8

# Another example
print("Start test foward pass update \n\n")

mu1 = 0.1
lam2 = 0.5
mu2 = 0.5
lam3 = 0.9
mu3 = 0.9
lam4 = 0
params = [mu1, lam2, mu2, lam3, mu3, lam4]



paramsOpt, listObjVals, listGrads, listChangefObj = blockCoordinateGradient(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol, plotSteps = False)

plotResults(x0, listB0k, listxk, listBk, params, paramsOpt, listObjVals, listGrads, listChangefObj)

# # Another example

print("Start test foward pass update \n\n")


mu1 = 0.5
lam2 = 0.65
mu2 = 0.5
lam3 = 0.75
mu3 = 0.5
lam4 = 0.45
params = [mu1, lam2, mu2, lam3, mu3, lam4]
listIndices = [1.0, 1.5, 2, 1.0]


paramsOpt, listObjVals, listGrads, listChangefObj = blockCoordinateGradient(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol, plotSteps = False)

plotResults(x0, listB0k, listxk, listBk, params, paramsOpt, listObjVals, listGrads, listChangefObj)


print("Start test foward pass update \n\n")


mu1 = 0.5
lam2 = 0.65
mu2 = 0.5
lam3 = 0.75
mu3 = 0.5
lam4 = 0.45
params = [mu1, lam2, mu2, lam3, mu3, lam4]
listIndices = [3.5, 1, 0.5, 1.0]


paramsOpt, listObjVals, listGrads, listChangefObj = blockCoordinateGradient(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol, plotSteps = False)

plotResults(x0, listB0k, listxk, listBk, params, paramsOpt, listObjVals, listGrads, listChangefObj)

plt.close("all")

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

t0 = 14*pi/8
t1 = 15*pi/8

# Define the points
x0 = np.array([10*cos(t0), 10*sin(t0)])
x1 = np.array([10*cos(t1), 10*sin(t1)])
x2 = np.array([x1[0] - 1.5, x1[1]])

# Their derivatives
B01 = np.array([-sin(t0), cos(t0)])
B1 = np.array([-sin(t1), cos(t1)])
B02 = (x2 - x0)/norm(x2-x0) # Because x2 is NOT on the boundary
B2 = np.copy(B02) # Because x2 is not on the boundary

listIndices = [eta2, eta2]
listxk = [x0, x1, x2]
listB0k = [B01, B02]
listBk = [B1, B2]

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
params0 = [mu1, lam2]
f0 = fObj_noTops(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

itt.plotFann(x0, listB0k, listxk, listBk, params = params0)
plt.title("Initial parameters")

paramsOpt, listObjVals, listGrads, listChangefObj, listIterates, listGradIterations = blockCoordinateGradient(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, maxIter, tol, plotSteps = False, saveIterates = True)

plotResults(x0, listB0k, listxk, listBk, params0, paramsOpt, listObjVals, listGrads, listChangefObj, trueSol = T2)




# Since this is a just a 2d problem we can plot the objective function and see how well
# this models the situation and if the optimization method converges to the minimum of
# our model

mus1 , lams2 = np.meshgrid( np.linspace(-0.1, 1.1, 200), np.linspace(-0.1, 1.1, 200) )
fObjMesh = np.empty(mus1.shape)

for i in range(200):
    for j in range(200):
        params = [mus1[i,j], lams2[i,j]]
        fObjMesh[i, j] = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

# plot it
# 2d
fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
plt.imshow( fObjMesh, cmap = colormap2, extent = [-0.1, 1.1, -0.1, 1.1], origin = "lower")
for i in range(len(listIterates)):
    pointk = listIterates[i]
    plt.scatter(pointk[0], pointk[1], marker = "o", s = 0.75, color = "#999999")
plt.scatter(params0[0], params0[1], marker = "o", color = "white", label = "starting point")
plt.scatter(paramsOpt[0], paramsOpt[1], marker = "*", color = "white", label = "optimum found")
plt.title("Objective function")
plt.xlabel("mu1")
plt.ylabel("lam2")
plt.legend()

# 3d
ax = plt.figure(figsize=(800/96, 800/96), dpi=96).add_subplot(projection='3d')
ax.plot_surface(mus1, lams2, fObjMesh, cmap = 'cet_linear_bmy_10_95_c78', edgecolor='royalblue', lw=0.5, rstride=8, cstride=8, alpha=0.3)
for i in range(len(listIterates)):
    pointk = listIterates[i]
    fk = fObj_noTops(pointk, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    ax.scatter(pointk[0], pointk[1], fk, c = "#565656", marker = "o", s = 2)
ax.scatter(paramsOpt[0], paramsOpt[1], listObjVals[-1], label = "optimum found", c = "black", marker = "*")
ax.scatter(params0[0], params0[1], f0, label = "starting point", c = "black")
ax.legend()




###########################################
print("\n\n\n x0 and x1 on the boundary far\n\n")

xSource = np.array([-15, -10])
center = np.array([0,0])
R = 10.0
eta1 = 1.0
eta2 = 1.452

t0 = 12*pi/8
t1 = 15*pi/8

# Define the points
x0 = np.array([10*cos(t0), 10*sin(t0)])
x1 = np.array([10*cos(t1), 10*sin(t1)])
x2 = np.array([x1[0] - 3.5, x1[1]])

# Their derivatives
B01 = np.array([-sin(t0), cos(t0)])
B1 = np.array([-sin(t1), cos(t1)])
B02 = (x2 - x0)/norm(x2-x0) # Because x2 is NOT on the boundary
B2 = np.copy(B02) # Because x2 is not on the boundary

listIndices = [eta2, eta2]
listxk = [x0, x1, x2]
listB0k = [B01, B02]
listBk = [B1, B2]

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
params0 = [mu1, lam2]
f0 = fObj_noTops(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

itt.plotFann(x0, listB0k, listxk, listBk, params = params0)
plt.title("Initial parameters")

paramsOpt, listObjVals, listGrads, listChangefObj, listIterates, listGradIterations = blockCoordinateGradient(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, maxIter, tol, saveIterates = True)

plotResults(x0, listB0k, listxk, listBk, params0, paramsOpt, listObjVals, listGrads, listChangefObj, trueSol = T2)


# Since this is a just a 2d problem we can plot the objective function and see how well
# this models the situation and if the optimization method converges to the minimum of
# our model

mus1 , lams2 = np.meshgrid( np.linspace(-0.1, 1.1, 200), np.linspace(-0.1, 1.1, 200) )
fObjMesh = np.empty(mus1.shape)

for i in range(200):
    for j in range(200):
        params = [mus1[i,j], lams2[i,j]]
        fObjMesh[i, j] = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

# plot it
# 2d
fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
plt.imshow( fObjMesh, cmap = colormap2, extent = [-0.1, 1.1, -0.1, 1.1], origin = "lower")
for i in range(len(listIterates)):
    pointk = listIterates[i]
    plt.scatter(pointk[0], pointk[1], marker = "o", s = 0.75, color = "#999999")
plt.scatter(params0[0], params0[1], marker = "o", color = "white", label = "starting point")
plt.scatter(paramsOpt[0], paramsOpt[1], marker = "*", color = "white", label = "optimum found")
plt.title("Objective function")
plt.xlabel("mu1")
plt.ylabel("lam2")
plt.legend()


# 3d
ax = plt.figure(figsize=(800/96, 800/96), dpi=96).add_subplot(projection='3d')
ax.plot_surface(mus1, lams2, fObjMesh, cmap = 'cet_linear_bmy_10_95_c78', edgecolor='royalblue', lw=0.5, rstride=8, cstride=8, alpha=0.3)
for i in range(len(listIterates)):
    pointk = listIterates[i]
    fk = fObj_noTops(pointk, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    ax.scatter(pointk[0], pointk[1], fk, c = "#565656", marker = "o", s = 2)
ax.scatter(paramsOpt[0], paramsOpt[1], listObjVals[-1], label = "optimum found", c = "black", marker = "*")
ax.scatter(params0[0], params0[1], f0, label = "starting point", c = "black")
ax.legend()




###########################################
print("\n\n\n x0 and x1 on the boundary far far\n\n")

xSource = np.array([-15, -10])
center = np.array([0,0])
R = 10.0
eta1 = 1.0
eta2 = 1.452

t0 = 10*pi/8
t1 = 15*pi/8

# Define the points
x0 = np.array([10*cos(t0), 10*sin(t0)])
x1 = np.array([10*cos(t1), 10*sin(t1)])
x2 = np.array([x1[0] - 3.5, x1[1]])

# Their derivatives
B01 = np.array([-sin(t0), cos(t0)])
B1 = np.array([-sin(t1), cos(t1)])
B02 = (x2 - x0)/norm(x2-x0) # Because x2 is NOT on the boundary
B2 = np.copy(B02) # Because x2 is not on the boundary

listIndices = [eta2, eta2]
listxk = [x0, x1, x2]
listB0k = [B01, B02]
listBk = [B1, B2]

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
params0 = [mu1, lam2]
f0 = fObj_noTops(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

itt.plotFann(x0, listB0k, listxk, listBk, params = params0)
plt.title("Initial parameters")

paramsOpt, listObjVals, listGrads, listChangefObj, listIterates, listGradIterations = blockCoordinateGradient(params0, x0, T0, grad0, x1, T1, grad1, x2, listIndices, listxk, listB0k, listBk, maxIter, tol, saveIterates = True)

plotResults(x0, listB0k, listxk, listBk, params0, paramsOpt, listObjVals, listGrads, listChangefObj, trueSol = T2)



# Since this is a just a 2d problem we can plot the objective function and see how well
# this models the situation and if the optimization method converges to the minimum of
# our model

mus1 , lams2 = np.meshgrid( np.linspace(-0.1, 1.1, 200), np.linspace(-0.1, 1.1, 200) )
fObjMesh = np.empty(mus1.shape)

for i in range(200):
    for j in range(200):
        params = [mus1[i,j], lams2[i,j]]
        fObjMesh[i, j] = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

# plot it
# 2d
fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
plt.imshow( fObjMesh, cmap = colormap2, extent = [-0.1, 1.1, -0.1, 1.1], origin = "lower")
for i in range(len(listIterates)):
    pointk = listIterates[i]
    plt.scatter(pointk[0], pointk[1], marker = "o", s = 0.75, color = "#999999")
plt.scatter(params0[0], params0[1], marker = "o", color = "white", label = "starting point")
plt.scatter(paramsOpt[0], paramsOpt[1], marker = "*", color = "white", label = "optimum found")
plt.title("Objective function")
plt.xlabel("mu1")
plt.ylabel("lam2")
plt.legend()

# 3d

ax = plt.figure(figsize=(800/96, 800/96), dpi=96).add_subplot(projection='3d')
ax.plot_surface(mus1, lams2, fObjMesh, cmap = 'cet_linear_bmy_10_95_c78', edgecolor='royalblue', lw=0.5, rstride=8, cstride=8, alpha=0.3)
for i in range(len(listIterates)):
    pointk = listIterates[i]
    fk = fObj_noTops(pointk, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    ax.scatter(pointk[0], pointk[1], fk, c = "#565656", marker = "o", s = 2)
ax.scatter(paramsOpt[0], paramsOpt[1], listObjVals[-1], label = "optimum found", c = "black", marker = "*")
ax.scatter(params0[0], params0[1], f0, label = "starting point", c = "black")
ax.legend()


