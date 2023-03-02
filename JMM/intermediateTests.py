import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi

plt.ion()
my_dpi = 96

def archlength(a, b):
     chordLength = norm(a-b)
     return 20*np.arcsin(chordLength/20)

def archlength_small(a, b):
      chordLength = norm(a-b)
      return 2*np.arcsin(chordLength/2)

def gradientBoundary(param, xFrom, Bfrom, xTo, Bto):
     '''
     Hermite interpolation for the boundary
     '''
     return 3*param**2*(2*xFrom + Bfrom - 2*xTo + Bto) + 2*param*(-3*xFrom - 2*Bfrom + 3*xTo - Bto) + Bfrom

def hermite_boundary(param, xFrom, Bfrom, xTo, Bto):
     return (2*xFrom + Bfrom - 2*xTo + Bto)*param**3 + (-3*xFrom - 2*Bfrom + 3*xTo - Bto)*param**2 + Bfrom*param + xFrom

 

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
     plt.scatter(x0[0], x0[1], s = 10, c = "#0036ff", label = "x0")
     plt.scatter(x1[0], x1[1], s = 7, c = "#001871", label = "x1")
     plt.scatter(x2[0], x2[1], s = 7, c = "#001871", label = "x2")
     plt.scatter(x3[0], x3[1], s = 7, c = "#001871", label = "x3")
     plt.scatter(xHat[0], xHat[1], s = 10, c = "#8000ff", label = "xHat")
     # If given the information plot the points on the boundaries
     if(mu1 is not None):
          xmu1 = hermite_boundary(mu1, x0, B01, x1, B1)
          plt.scatter(xmu1[0], xmu1[1], s = 5, c = "#117783")
     if(mu2 is not None):
          xmu2 = hermite_boundary(mu2, x0, B02, x2, B2)
          plt.scatter(xmu2[0], xmu2[1], s = 5, c = "#117783")
     if(mu3 is not None):
          xmu3 = hermite_boundary(mu3, x0, B03, x3, B3)
          plt.scatter(xmu3[0], xmu3[1], s = 5, c = "#117783")
     if(lam2 is not None):
          xlam2 = hermite_boundary(lam2, x0, B02, x2, B2)
          plt.scatter(xlam2[0], xlam2[1], s = 5, c = "#561183")
     if(lam3 is not None):
          xlam3 = hermite_boundary(lam3, x0, B03, x3, B3)
          plt.scatter(xlam3[0], xlam3[1], s = 5, c = "#561183")
     if(lam4 is not None):
          xlam4 = hermite_boundary(lam4, x0, B0Hat, xHat, BHat)
          plt.scatter(xlam4[0], xlam4[1], s = 5, c = "#561183")
     if(mu1 is not None and lam2 is not None):
          plt.plot([xmu1[0], xlam2[0]], [xmu1[1], xlam2[1]], linewidth = 1, c = "#525460")
     if(mu2 is not None and lam3 is not None):
          plt.plot([xmu2[0], xlam3[0]], [xmu2[1], xlam3[1]], linewidth = 1, c = "#525460")
     if(mu3 is not None and lam4 is not None):
          plt.plot([xmu3[0], xlam4[0]], [xmu3[1], xlam4[1]], linewidth = 1, c = "#525460")
          
     plt.legend()
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
