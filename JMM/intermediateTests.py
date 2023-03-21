import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi

plt.ion()
my_dpi = 96
a = 800
b = 800

def archlength(a, b):
     chordLength = norm(a-b)
     return 20*np.arcsin(chordLength/20)

def archlength_small(a, b):
      chordLength = norm(a-b)
      return 2*np.arcsin(chordLength/2)

def gradientBoundary(param, xFrom, Bfrom, xTo, Bto):
     '''
     Tangent to the boundary (interpolated using Hermite)
     '''
     return 3*(2*xFrom + Bfrom - 2*xTo + Bto)*param**2 + 2*(-3*xFrom - 2*Bfrom + 3*xTo - Bto)*param + Bfrom

def hermite_boundary(param, xFrom, Bfrom, xTo, Bto):
     '''
     Hermite interpolation of the boundary
     '''
     return (2*xFrom + Bfrom - 2*xTo + Bto)*param**3 + (-3*xFrom - 2*Bfrom + 3*xTo - Bto)*param**2 + Bfrom*param + xFrom

def arclengthSimpson(mu, lam, xFrom, Bfrom, xTo, Bto):
     '''
     arclength along a boundary from xLam to xMu
     '''
     Bmu = gradientBoundary(mu, xFrom, Bfrom, xTo, Bto)
     Blam = gradientBoundary(lam, xFrom, Bfrom, xTo, Bto)
     B_mid = gradientBoundary((mu + lam)/2, xFrom, Bfrom, xTo, Bto)
     return (norm(Bmu) + 4*norm(B_mid) + norm(Blam))/6


def tMuMin(mu, xA, xR, BR, xHat, BHat):
     xMu = hermite_boundary(mu, xR, BR, xHat, BHat)
     Bmu = gradientBoundary(mu, xR, BR, xHat, BHat)
     Bmu_perp = np.array([Bmu[1], -Bmu[0]])
     return np.dot(xB - xMu, Bmu_perp)

def tMu(mu, xA, xB, xR, BR, xHat, BHat):
     xMu = hermite_boundary(mu, xR, BR, xHat, BHat)
     Bmu = gradientBoundary(mu, xR, BR, xHat, BHat)
     Bmu_perp = np.array([Bmu[1], -Bmu[0]])
     return np.dot(xMu - xA, Bmu_perp)/np.dot(xB - xA, Bmu_perp)
 

def plotFan3(x0, B01, B02, B03, B0Hat, x1, B1, x2, B2, x3, B3, xHat, BHat,
             mu1 = None, lam2 = None, mu2 = None, lam3 = None,
             mu3 = None, lam4 = None):
     '''
     Plots a triangle fan with 3 intermediate points to xHat
     '''
     axMin = min(x0[0], x1[0], x2[0], x3[0], xHat[0],
                 x0[1], x1[1], x2[1], x3[1], xHat[1])
     axMax = max(x0[0], x1[0], x2[0], x3[0], xHat[0],
                 x0[1], x1[1], x2[1], x3[1], xHat[1])
     fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
     # Plot the triangles in the triangulation (straight lines)
     plt.plot([x0[0], x1[0], x2[0], x0[0]], [x0[1], x1[1], x2[1], x0[1]], linewidth = 0.5, alpha = 0.8, c = "#aaaaaa")
     plt.plot([x0[0], x3[0], x2[0]], [x0[1], x3[1], x2[1]], linewidth = 0.5, alpha = 0.8, c = "#aaaaaa")
     plt.plot([x0[0], xHat[0], x3[0]], [x0[1], xHat[1], x3[1]], linewidth = 0.5, alpha = 0.8, c = "#aaaaaa")
     # Compute the curvy boundaries
     params = np.linspace(0, 1, 100)
     x0x1 = np.empty((100, 2))
     x0x2 = np.empty((100, 2))
     x0x3 = np.empty((100, 2))
     x0xHat =  np.empty((100,2))
     for i in range(100):
          x0x1[i, :] = hermite_boundary(params[i], x0, B01, x1, B1)
          x0x2[i, :] = hermite_boundary(params[i], x0, B02, x2, B2)
          x0x3[i, :] = hermite_boundary(params[i], x0, B03, x3, B3)
          x0xHat[i, :] = hermite_boundary(params[i], x0, B0Hat, xHat, BHat)
     # Plot the curvy boundaries
     plt.plot(x0x1[:, 0], x0x1[:,1], linewidth = 1, c = "#483e64")
     plt.plot(x0x2[:, 0], x0x2[:,1], linewidth = 1, c = "#483e64")
     plt.plot(x0x3[:, 0], x0x3[:,1], linewidth = 1, c = "#483e64")
     plt.plot(x0xHat[:, 0], x0xHat[:,1], linewidth = 1, c = "#483e64")
     # The points we are interested in
     plt.scatter(x0[0], x0[1], s = 15, c = "#0027ff", label = "x0")
     plt.scatter(x1[0], x1[1], s = 7, c = "#001871", label = "x1")
     plt.scatter(x2[0], x2[1], s = 7, c = "#001871", label = "x2")
     plt.scatter(x3[0], x3[1], s = 7, c = "#001871", label = "x3")
     plt.scatter(xHat[0], xHat[1], s = 10, c = "#8000ff", label = "xHat")
     # If given the information plot the points on the boundaries
     if(mu1 is not None):
          xmu1 = hermite_boundary(mu1, x0, B01, x1, B1)
          #print("xmu1:", xmu1)
          plt.scatter(xmu1[0], xmu1[1], s = 5, c = "#117783")
     if(mu2 is not None):
          xmu2 = hermite_boundary(mu2, x0, B02, x2, B2)
          #print("xmu2:", xmu2)
          plt.scatter(xmu2[0], xmu2[1], s = 5, c = "#117783")
     if(mu3 is not None):
          xmu3 = hermite_boundary(mu3, x0, B03, x3, B3)
          #print("xmu3: ", xmu3)
          plt.scatter(xmu3[0], xmu3[1], s = 5, c = "#117783")
     if(lam2 is not None):
          xlam2 = hermite_boundary(lam2, x0, B02, x2, B2)
          #print("xlam2: ", xlam2)
          plt.scatter(xlam2[0], xlam2[1], s = 5, c = "#561183")
     if(lam3 is not None):
          xlam3 = hermite_boundary(lam3, x0, B03, x3, B3)
          #print("xlam3: ", xlam3)
          plt.scatter(xlam3[0], xlam3[1], s = 5, c = "#561183")
     if(lam4 is not None):
          xlam4 = hermite_boundary(lam4, x0, B0Hat, xHat, BHat)
          #print("xlam4: ", xlam4)
          plt.scatter(xlam4[0], xlam4[1], s = 5, c = "#561183")
     if(mu1 is not None and lam2 is not None):
          plt.plot([xmu1[0], xlam2[0]], [xmu1[1], xlam2[1]], linewidth = 1, c = "#0024ff")
     if(mu2 is not None and lam3 is not None):
          plt.plot([xmu2[0], xlam3[0]], [xmu2[1], xlam3[1]], linewidth = 1, c = "#0024ff")
     if(mu3 is not None and lam4 is not None):
          plt.plot([xmu3[0], xlam4[0]], [xmu3[1], xlam4[1]], linewidth = 1, c = "#0024ff")
          
     plt.legend()
     plt.xlim(axMin - abs(0.2*axMax), axMax + abs(0.2*axMax))
     plt.ylim(axMin - abs(0.2*axMax), axMax + abs(0.2*axMax))
     ax = plt.gca()
     ax.set_aspect("equal")



x0 = np.array([0.0,0.0])
x1 = np.array([2, -0.2])
x2 = np.array([1.5, 0.8])
x3 = np.array([0.2, 1.2])
xHat = np.array([-0.8, 0.2])
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
BHat = np.array([-1, -0.4])
BHat = BHat/norm(BHat)


def drawOneCurvyTriangle(x0, B0k, B0k1, xk, xk1, Bk, Bk1, BkTop = None, Bk1Top = None, title = None, tangents = False):
     '''
     Draw one curvy triangle x0 xk x(k+1)
     '''
     fig = plt.figure(figsize=(a/my_dpi, b/my_dpi), dpi=my_dpi)
     # Plot the triangles in the triangulation (straight lines)
     plt.plot([x0[0], xk[0], xk1[0], x0[0]], [x0[1], xk[1], xk1[1], x0[1]], linewidth = 0.5, alpha = 0.8, c = "#aaaaaa")
     # Compute the curvy boundaries
     params = np.linspace(0,1,100)
     x0xk = np.empty((100,2))
     x0xk1 = np.empty((100,2))
     if( BkTop is not None and Bk1Top is not None ):
          # if the curvy triangle has also a curvy top
          xkxk1 = np.empty((100,2))
          for i in range(100):
               xkxk1[i, :] = hermite_boundary(params[i], xk, BkTop, xk1, Bk1Top)
     for i in range(100):
          x0xk[i, :] = hermite_boundary(params[i], x0, B0k, xk, Bk)
          x0xk1[i, :] = hermite_boundary(params[i], x0, B0k1, xk1, Bk1)
     # Plot the curvy boundaries
     plt.plot(x0xk[:, 0], x0xk[:, 1], linewidth = 2, c = "#483e64")
     plt.plot(x0xk1[:, 0], x0xk1[:,1], linewidth = 2, c = "#483e64")
     if(BkTop is not None and Bk1Top is not None):
          plt.plot(xkxk1[:,0], xkxk1[:,1], linewdith = 2, c = "#3e6451")
     # Plot the points we are interested in
     plt.scatter(x0[0], x0[1], s = 20, c = "#0027ff", label = "x0")
     plt.scatter([xk[0], xk1[0]], [xk[1], xk1[1]],  s = 12, c = "#001871", label = "xk, x(k+1)")
     # If specified, plot the tangents
     if(tangents):
          plt.arrow(x0[0], x0[1], B0k[0], B0k[1], linewidth = 0.2, width = 0.005, color = "#858cb7", label = "B0k")
          plt.arrow(x0[0], x0[1], B0k1[0], B0k1[1], linewidth = 0.2, width = 0.005, color = "#858cb7", label = "B0k1")
          plt.arrow(xk[0], xk[1], Bk[0], Bk[1], linewidth = 0.2, width = 0.005, color = "#858cb7", label = "Bk")
          plt.arrow(xk1[0], xk1[1], Bk1[0], Bk1[1], linewidth = 0.2, width = 0.005, color = "#858cb7", label = "Bk1")
     if(title is not None):
          plt.title(title)
     plt.legend()
     ax = plt.gca()
     ax.set_aspect("equal")
     
## For type 1
x0 = np.array([0.0,0.0])
xk = np.array([2, -0.2])
xk1 = np.array([1.5, 0.8])
B0k = np.array([2.2, 1])
B0k = B0k/norm(B0k)
B0k1 = np.array([1, 1.5])
B0k1 = B0k1/norm(B0k1)
Bk = np.array([1, -0.6])
Bk = Bk/norm(Bk)
Bk1 = np.array([2, -0.2])
Bk1 = Bk1/norm(Bk1)

##Plot feasible region
drawOneCurvyTriangle(x0, B0k, B0k1, x1, x2, Bk, Bk1, title = "Type 1", tangents = False)
lambdasFeas = np.linspace(0.252635, 1, 100)
regFes = np.empty((100,2))
yk1 = hermite_boundary(0.252635, x0, B0k1, xk1, Bk1)
for i in range(100):
     regFes[i, :] = hermite_boundary(lambdasFeas[i], x0, B0k1, xk1, Bk1)
plt.plot(regFes[:, 0], regFes[:, 1], alpha = 0.5, linewidth = 5, c = "#8bffff", label = "lambdas feasible")
plt.plot([1.71137263, 0.26964767], [-0.10010285, 0.2495371], linewidth = 0.8, c = "#656565", linestyle = "--")
plt.scatter(yk1[0], yk1[1],  c = "#8bffff", label = "lambdaMin")
plt.scatter(1.71137263, -0.10010285, s = 15, c = "#09c700", label = "ykPrime")
plt.legend()
plt.title("Type 1, feasible region")

## For type 2
x0 = np.array([0.0,0.0])
xk = np.array([2, -0.2])
xk1 = np.array([1.5, 0.8])
B0k = np.array([0, -1])
B0k = B0k/norm(B0k)
B0k1 = np.array([0.4, 0.1])
B0k1 = B0k1/norm(B0k1)
Bk = np.array([1,1])
Bk = Bk/norm(Bk)
Bk1 = np.array([1, 2])
Bk1 = Bk1/norm(Bk1)
drawOneCurvyTriangle(x0, B0k, B0k1, x1, x2, Bk, Bk1, title = "Type 2", tangents = False)
ykPrime = hermite_boundary(0.1, x0, B0k, xk, Bk)
yk1 = hermite_boundary(0.589507, x0, B0k1, xk1, Bk1)
lambdasFeas = np.linspace(0, 0.589507, 100)
regFes = np.empty((100,2))
for i in range(100):
     regFes[i, :] = hermite_boundary(lambdasFeas[i], x0, B0k1, xk1, Bk1)
plt.plot(regFes[:, 0], regFes[:, 1], alpha = 0.5, linewidth = 5, c = "#8bffff", label = "lambdas feasible")
plt.plot([ykPrime[0], yk1[0]], [ykPrime[1], yk1[1]], linewidth = 0.8, c = "#656565", linestyle = "--")
plt.scatter(yk1[0], yk1[1],  c = "#8bffff", label = "lambdaMax")
plt.scatter(ykPrime[0], ykPrime[1], s = 15, c = "#09c700", label = "ykPrime")
plt.legend()
plt.title("Type 2, feasible region")


## For type 3
# x0 = np.array([0.0,0.0])
# xk = np.array([2, -0.2])
# xk1 = np.array([1.5, 0.8])
# B0k = np.array([0, -1])
# B0k = B0k/norm(B0k)
# B0k1 = np.array([1, 1.5])
# B0k1 = B0k1/norm(B0k1)
# Bk = np.array([1,1])
# Bk = Bk/norm(Bk)
# Bk1 = np.array([2, -0.2])
# Bk1 = Bk1/norm(Bk1)

## For type 4
x0 = np.array([0.0,0.0])
xk = np.array([2, -0.2])
xk1 = np.array([1.5, 0.8])
B0k = np.array([0.4, 0.1])
B0k = B0k/norm(B0k)
B0k1 = np.array([0.4, 0.1])
B0k1 = B0k1/norm(B0k1)
Bk = np.array([1, -1])
Bk = Bk/norm(Bk)
Bk1 = np.array([1, 2])
Bk1 = Bk1/norm(Bk1)

drawOneCurvyTriangle(x0, B0k, B0k1, x1, x2, Bk, Bk1, title = "Type 4", tangents = False)
lamMin = 0.141814
lamMax = 0.900304
ykPrime = hermite_boundary(0.4, x0, B0k, xk, Bk)
yk1Min = hermite_boundary(lamMin, x0, B0k1, xk1, Bk1)
yk1Max = hermite_boundary(lamMax, x0, B0k1, xk1, Bk1)
lambdasFeas = np.linspace(lamMin, lamMax, 100)
regFes = np.empty((100,2))
for i in range(100):
     regFes[i, :] = hermite_boundary(lambdasFeas[i], x0, B0k1, xk1, Bk1)
plt.plot(regFes[:, 0], regFes[:, 1], alpha = 0.5, linewidth = 5, c = "#8bffff", label = "lambdas feasible")
plt.plot([ykPrime[0], yk1Max[0]], [ykPrime[1], yk1Max[1]], linewidth = 0.8, c = "#656565", linestyle = "--")
plt.plot([ykPrime[0], yk1Min[0]], [ykPrime[1], yk1Min[1]], linewidth = 0.8, c = "#656565", linestyle = "--")
plt.scatter(yk1Max[0], yk1Max[1],  c = "#8bffff", label = "lambdaMax")
plt.scatter(yk1Min[0], yk1Min[1],  c = "#8bffff", label = "lambdaMin")
plt.scatter(ykPrime[0], ykPrime[1], s = 15, c = "#09c700", label = "ykPrime")
plt.legend()
plt.title("Type 4, feasible region")

