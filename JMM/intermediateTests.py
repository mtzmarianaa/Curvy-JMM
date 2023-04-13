import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin

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

def secondDer_Boundary(param, xFrom, Bfrom, xTo, Bto):
     '''
     d2B/dparam2
     '''
     return 6*(2*xFrom + Bfrom - 2*xTo + Bto)*param + 2*(-3*xFrom - 2*Bfrom + 3*xTo - Bto)

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
     plt.plot(x0x1[:, 0], x0x1[:,1], linewidth = 1, c = "#aaaaaa")
     plt.plot(x0x2[:, 0], x0x2[:,1], linewidth = 1, c = "#aaaaaa")
     plt.plot(x0x3[:, 0], x0x3[:,1], linewidth = 1, c = "#aaaaaa")
     plt.plot(x0xHat[:, 0], x0xHat[:,1], linewidth = 1, c = "#aaaaaa")
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

def evaluateCreepingRay(paramMin, paramMax, xFrom, BFrom, xTo, BTo, nEvals = 100):
     # Returns a (nEvals, 2) array of points on the boundary
     paramsToEval = np.linspace(paramMin, paramMax, nEvals)
     xFromxTo = np.empty((nEvals,2))
     for i in range(nEvals):
          xFromxTo[i, :] = hermite_boundary(paramsToEval[i], xFrom, BFrom, xTo, BTo)
     return xFromxTo

def plotFann(x0, listB0k, listxk, listBk, listBkBk1 = None, params = None, indCrTop = None, paramsCrTop = None, indStTop = None, paramsStTop = None):
     '''
     Plots a triangle fan with n regions (i.e. n triangles)
     '''
     # Set indStTop if indCrTop is given
     if(paramsCrTop is not None and paramsStTop is None):
          indStTop = [0]
          paramsStTop = [0,0]
     if(paramsStTop is not None and paramsCrTop is None):
          indCrTop = [0]
          paramsCrTop = [0,0]
     axMin = np.min(listxk) # To set the ax limits
     axMax = np.max(listxk) # To set the ax limits
     nRegions = len(listxk) - 2 # Number of triangles
     fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
     # Plot the triangles in the triangulation (straight lines)
     x0 = listxk[0]
     x1 = listxk[1]
     # Plot "base" straight line
     plt.plot([x0[0], x1[0]], [x0[1], x1[1]], linewidth = 0.4, alpha = 0.6, c = "#aaaaaa")
     # Plot the "base" curvy line
     params_inter = np.linspace(0,1, 100)
     x0x1 = np.empty((100, 2))
     for i in range(100):
          x0x1[i, :] = hermite_boundary(params_inter[i], x0, listB0k[0], x1, listBk[0])
     plt.plot(x0x1[:, 0], x0x1[:, 1], linewidth = 0.8, alpha = 0.8, c = "#616276")
     # Plot x0
     plt.scatter(x0[0], x0[1], s = 15, c = "#0027ff", label = "x0")
     plt.scatter(x1[0], x1[1], s = 7, c = "#001871", label = "x1")
     # Plot the rest of the edges of the triangles
     for j in range (nRegions):
          # j from 0 to nRegions - 1
          # PLOT THE STRAIGHT LINE TRIANGLES
          xkM1 = listxk[j+1]
          xk = listxk[j+2]
          # Plot the top of the triangle
          plt.plot([x0[0], xk[0]], [x0[1], xk[1]], linewidth = 0.4, alpha = 0.6, c = "#aaaaaa")
          # Plot the side of the triangle
          plt.plot([xkM1[0], xk[0]], [xkM1[1], xk[1]], linewidth = 0.4, alpha = 0.6, c = "#aaaaaa")
          # PLOT THE CURVY TRIANGLES
          x0xk = np.empty((100, 2))
          for i in range(100):
               x0xk[i, :] = hermite_boundary(params_inter[i], x0, listB0k[j+1], xk, listBk[j+1])
          plt.plot(x0xk[:, 0], x0xk[:, 1], linewidth = 0.8, alpha = 0.8, c = "#616276")
          plt.scatter(xk[0], xk[1], s = 7, c = "#001871", label = "x" + str(j+2))
     # If given tangents for the triangle tops then plot those curve edges as well
     if( listBkBk1 is not None ):
          # Meaning that we have a general triangle fan
          for j in range(nRegions):
               k = 2*j
               xkM1 = listxk[j+1]
               xk = listxk[j+2]
               #Plot the curvy top
               xkxk1 = np.empty((100,2))
               for i in range(100):
                    xkxk1[i, :] = hermite_boundary(params_inter[i], xkM1, listBkBk1[k], xk, listBkBk1[k+1])
               plt.plot(xkxk1[:,0], xkxk1[:,1], linewidth = 0.8, alpha = 0.8, c = "#616276")
     # Now, if given information about points on the curvy triangles plot straight lines
     if(params is not None and paramsCrTop is None and paramsStTop is None):
          # This means that we only have a standard update, no creeping or going trought boundaries on the top
          for j in range(nRegions):
               k = 2*j
               muk = params[k]
               lamk1 = params[k+1]
               B0k = listB0k[j]
               xk = listxk[j+1]
               Bk = listBk[j]
               B0k1 = listB0k[j+1]
               xk1 = listxk[j+2]
               Bk1 = listBk[j+1]
               xmuk = hermite_boundary(muk, x0, B0k, xk, Bk)
               xlamk1 = hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
               # Plot a straight line between them
               plt.plot([xmuk[0], xlamk1[0]], [xmuk[1], xlamk1[1]], linewidth = 1.2, c = "#0024ff")
               # Scatter them
               plt.scatter([xmuk[0], xlamk1[0]], [xmuk[1], xlamk1[1]], s = 5, c = "#117783")
               if (j>0):
                    crRay = evaluateCreepingRay(min(params[k-1], params[k]), max(params[k-1], params[k]), x0, B0k, xk, Bk)
                    plt.plot(crRay[:,0], crRay[:,1], linewidth = 1.2, c = "#0024ff")
     elif(params is not None and (paramsCrTop is not None or paramsStTop is not None)):
          # This means that we have standard updates and creeping or going through boundaries on the top
          currentCrTop = 0
          currentStTop = 0
          for j in range(nRegions):
               k = 2*j
               muk = params[k]
               lamk1 = params[k+1]
               B0k = listB0k[j]
               xk = listxk[j+1]
               Bk = listBk[j]
               BkBk1_0 = listBkBk1[k] # grad of hkk1 at xk
               B0k1 = listB0k[j+1]
               xk1 = listxk[j+2]
               Bk1 = listBk[j+1]
               BkBk1_1 = listBkBk1[k+1] # grad of hkk1 at xk1
               # Compute the points
               xmuk = hermite_boundary(muk, x0, B0k, xk, Bk)
               xlamk1 = hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
               nTop = j + 1 # Number of boundary on the top that we are considering
               if (nTop == indCrTop[currentCrTop]):
                    # This means that there is creeping along this triangle top
                    rk = paramsCrTop[2*currentCrTop]
                    sk = paramsCrTop[2*currentCrTop + 1]
                    ak = hermite_boundary(rk, xk, BkBk1_0, xk1, BkBk1_1)
                    bk = hermite_boundary(sk, xk, BkBk1_0, xk1, BkBk1_1)
                    # Since there is creeping then we just need to plot straight lines
                    # between xmuk - ak and bk - xlamk1
                    plt.plot([xmuk[0], ak[0]], [xmuk[1], ak[1]], linewidth = 1.2, c = "#0024ff")
                    plt.scatter([xmuk[0], ak[0]], [xmuk[1], ak[1]], s = 5, c = "#117783")
                    plt.plot([bk[0], xlamk1[0]], [bk[1], xlamk1[1]], linewidth = 1.2, c = "#0024ff")
                    plt.scatter([bk[0], xlamk1[0]], [bk[1], xlamk1[1]], s = 5, c = "#117783")
                    # Plot the creeping update
                    crRay = evaluateCreepingRay(min(rk, sk), max(rk, sk), xk, BkBk1_0, xk1, BkBk1_1)
                    plt.plot(crRay[:,0], crRay[:,1], linewidth = 1.2, c = "#0024ff")
                    if( currentCrTop < len(indCrTop)-1 ):
                         currentCrTop += 1
               elif(nTop == indStTop[currentStTop]):
                    # This means that there is no creeping along this triangle top, it goes into the region
                    rk = paramsStTop[2*currentStTop]
                    sk = paramsStTop[2*currentStTop + 1]
                    ak = hermite_boundary(rk, xk, BkBk1_0, xk1, BkBk1_1)
                    bk = hermite_boundary(sk, xk, BkBk1_0, xk1, BkBk1_1)
                    # Now we need to plot straight lines between
                    # xmuk - ak, ak - bk, bk - xlamk1
                    plt.plot([xmuk[0], ak[0]], [xmuk[1], ak[1]], linewidth = 1.2, c = "#0024ff")
                    plt.plot([ak[0], bk[0]], [ak[1], bk[1]], linewidth = 1.2, c = "#0024ff")
                    plt.scatter([xmuk[0], ak[0]], [xmuk[1], ak[1]], s = 5, c = "#117783")
                    plt.plot([bk[0], xlamk1[0]], [bk[1], xlamk1[1]], linewidth = 1.2, c = "#0024ff")
                    plt.scatter([bk[0], xlamk1[0]], [bk[1], xlamk1[1]], s = 5, c = "#117783")
                    if (j>0):
                         crRay = evaluateCreepingRay(min(params[k-1], params[k]), max(params[k-1], params[k]), x0, B0k, xk, Bk)
                         plt.plot(crRay[:,0], crRay[:,1], linewidth = 1.2, c = "#0024ff")
                    if( currentStTop < len(indStTop)-1 ):
                         currentStTop += 1
               else:
                    # Meaning that there is no creeping or straight line on the triangle top
                    # We just need to plot straight lines between xmuk - xlamk1
                    plt.plot([xmuk[0], xlamk1[0]], [xmuk[1], xlamk1[1]], linewidth = 1.2, c = "#0024ff")
                    plt.scatter([xmuk[0], xlamk1[0]], [xmuk[1], xlamk1[1]], s = 5, c = "#117783")
                    if (j>0):
                         crRay = evaluateCreepingRay(min(params[k-1], params[k]), max(params[k-1], params[k]), x0, B0k, xk, Bk)
                         plt.plot(crRay[:,0], crRay[:,1], linewidth = 1.2, c = "#0024ff")
          # Plot the creeping that happens between xnRegions and xHat
          crRay = evaluateCreepingRay(params[2*nRegions - 1], 1, x0, B0k1, xk1, Bk1)
          plt.plot(crRay[:,0], crRay[:,1], linewidth = 1.2, c = "#0024ff")
     plt.legend()
     #plt.xlim(axMin - abs(0.2*axMax), axMax + abs(0.2*axMax))
     #plt.ylim(axMin - abs(0.2*axMax), axMax + abs(0.2*axMax))
     ax = plt.gca()
     ax.set_aspect("equal")
               
          


# x0 = np.array([0.0,0.0])
# x1 = np.array([2, -0.2])
# x2 = np.array([1.5, 0.8])
# x3 = np.array([0.2, 1.2])
# xHat = np.array([-0.8, 0.2])
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
# BHat = np.array([-1, -0.4])
# BHat = BHat/norm(BHat)


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
          plt.plot(xkxk1[:,0], xkxk1[:,1], linewidth = 2, c = "#3e6451")
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
     
# ## For type 1
# x0 = np.array([0.0,0.0])
# xk = np.array([2, -0.2])
# xk1 = np.array([1.5, 0.8])
# B0k = np.array([2.2, 1])
# B0k = B0k/norm(B0k)
# B0k1 = np.array([1, 1.5])
# B0k1 = B0k1/norm(B0k1)
# Bk = np.array([1, -0.6])
# Bk = Bk/norm(Bk)
# Bk1 = np.array([2, -0.2])
# Bk1 = Bk1/norm(Bk1)

# ##Plot feasible region
# drawOneCurvyTriangle(x0, B0k, B0k1, x1, x2, Bk, Bk1, title = "Type 1", tangents = False)
# lambdasFeas = np.linspace(0.252635, 1, 100)
# regFes = np.empty((100,2))
# yk1 = hermite_boundary(0.252635, x0, B0k1, xk1, Bk1)
# for i in range(100):
#      regFes[i, :] = hermite_boundary(lambdasFeas[i], x0, B0k1, xk1, Bk1)
# plt.plot(regFes[:, 0], regFes[:, 1], alpha = 0.5, linewidth = 5, c = "#8bffff", label = "lambdas feasible")
# plt.plot([1.71137263, 0.26964767], [-0.10010285, 0.2495371], linewidth = 0.8, c = "#656565", linestyle = "--")
# plt.scatter(yk1[0], yk1[1],  c = "#8bffff", label = "lambdaMin")
# plt.scatter(1.71137263, -0.10010285, s = 15, c = "#09c700", label = "ykPrime")
# plt.legend()
# plt.title("Type 1, feasible region")

# ## For type 2
# x0 = np.array([0.0,0.0])
# xk = np.array([2, -0.2])
# xk1 = np.array([1.5, 0.8])
# B0k = np.array([0, -1])
# B0k = B0k/norm(B0k)
# B0k1 = np.array([0.4, 0.1])
# B0k1 = B0k1/norm(B0k1)
# Bk = np.array([1,1])
# Bk = Bk/norm(Bk)
# Bk1 = np.array([1, 2])
# Bk1 = Bk1/norm(Bk1)
# drawOneCurvyTriangle(x0, B0k, B0k1, x1, x2, Bk, Bk1, title = "Type 2", tangents = False)
# ykPrime = hermite_boundary(0.1, x0, B0k, xk, Bk)
# yk1 = hermite_boundary(0.589507, x0, B0k1, xk1, Bk1)
# lambdasFeas = np.linspace(0, 0.589507, 100)
# regFes = np.empty((100,2))
# for i in range(100):
#      regFes[i, :] = hermite_boundary(lambdasFeas[i], x0, B0k1, xk1, Bk1)
# plt.plot(regFes[:, 0], regFes[:, 1], alpha = 0.5, linewidth = 5, c = "#8bffff", label = "lambdas feasible")
# plt.plot([ykPrime[0], yk1[0]], [ykPrime[1], yk1[1]], linewidth = 0.8, c = "#656565", linestyle = "--")
# plt.scatter(yk1[0], yk1[1],  c = "#8bffff", label = "lambdaMax")
# plt.scatter(ykPrime[0], ykPrime[1], s = 15, c = "#09c700", label = "ykPrime")
# plt.legend()
# plt.title("Type 2, feasible region")


# ## For type 3
# # x0 = np.array([0.0,0.0])
# # xk = np.array([2, -0.2])
# # xk1 = np.array([1.5, 0.8])
# # B0k = np.array([0, -1])
# # B0k = B0k/norm(B0k)
# # B0k1 = np.array([1, 1.5])
# # B0k1 = B0k1/norm(B0k1)
# # Bk = np.array([1,1])
# # Bk = Bk/norm(Bk)
# # Bk1 = np.array([2, -0.2])
# # Bk1 = Bk1/norm(Bk1)

# ## For type 4
# x0 = np.array([0.0,0.0])
# xk = np.array([2, -0.2])
# xk1 = np.array([1.5, 0.8])
# B0k = np.array([0.4, 0.1])
# B0k = B0k/norm(B0k)
# B0k1 = np.array([0.4, 0.1])
# B0k1 = B0k1/norm(B0k1)
# Bk = np.array([1, -1])
# Bk = Bk/norm(Bk)
# Bk1 = np.array([1, 2])
# Bk1 = Bk1/norm(Bk1)

# drawOneCurvyTriangle(x0, B0k, B0k1, x1, x2, Bk, Bk1, title = "Type 4", tangents = False)
# lamMin = 0.141814
# lamMax = 0.900304
# ykPrime = hermite_boundary(0.4, x0, B0k, xk, Bk)
# yk1Min = hermite_boundary(lamMin, x0, B0k1, xk1, Bk1)
# yk1Max = hermite_boundary(lamMax, x0, B0k1, xk1, Bk1)
# lambdasFeas = np.linspace(lamMin, lamMax, 100)
# regFes = np.empty((100,2))
# for i in range(100):
#      regFes[i, :] = hermite_boundary(lambdasFeas[i], x0, B0k1, xk1, Bk1)
# plt.plot(regFes[:, 0], regFes[:, 1], alpha = 0.5, linewidth = 5, c = "#8bffff", label = "lambdas feasible")
# plt.plot([ykPrime[0], yk1Max[0]], [ykPrime[1], yk1Max[1]], linewidth = 0.8, c = "#656565", linestyle = "--")
# plt.plot([ykPrime[0], yk1Min[0]], [ykPrime[1], yk1Min[1]], linewidth = 0.8, c = "#656565", linestyle = "--")
# plt.scatter(yk1Max[0], yk1Max[1],  c = "#8bffff", label = "lambdaMax")
# plt.scatter(yk1Min[0], yk1Min[1],  c = "#8bffff", label = "lambdaMin")
# plt.scatter(ykPrime[0], ykPrime[1], s = 15, c = "#09c700", label = "ykPrime")
# plt.legend()
# plt.title("Type 4, feasible region")


# # Plot generalized triangle fan

# print("Generalized triangle fan")

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


# mu1 = 0.15
# lam2 = 0.15
# mu2 = 0.13
# lam3 = 0.17
# mu3 = 0.15
# lam4 = 1
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
# listBkBk1 = [B1B2_0, B1B2_1, B2B3_0, B2B3_1, B3B4_0, B3B4_1]



# plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = None)
# plt.title("Generalized triangle fan")

# # Plot triangle fan without points on the tops

# mu1 = 0.7
# lam2 = 0.4
# mu2 = 0.6
# lam3 = 0.4
# mu3 = 0.1
# lam4 = 1
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4]


# plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0)
# plt.title("Generalized triangle fan, no points on triangle tops")

# # Plot triangle fan with points on the tops

# mu1 = 0.7
# lam2 = 0.4
# mu2 = 0.6
# lam3 = 0.4
# mu3 = 0.1
# lam4 = 1
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4]

# indCrTop = [1]
# paramsCrTop = [0.3, 0.5]

# plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop)
# plt.title("Generalized triangle fan, just creeping on tops")



# # Plot triangle fan with points on the tops

# mu1 = 0.7
# lam2 = 0.8
# mu2 = 0.1
# lam3 = 0.1
# mu3 = 0.4
# lam4 = 0.5
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4]

# indCrTop = [1]
# paramsCrTop = [0.3, 0.5]
# indStTop = [3]
# paramsStTop = [0.3, 0.7]

# plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop,
#          indStTop = indStTop, paramsStTop = paramsStTop)
# plt.title("Generalized triangle fan, creeping and going through tops")

# ############
# # Plot unfeasible paths


# mu1 = 0.7
# lam2 = 0.4
# mu2 = 0.6
# lam3 = 1
# mu3 = 0.8
# lam4 = 1
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4]

# indCrTop = [1]
# paramsCrTop = [0.3, 0.5]

# plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop)
# plt.title("Generalized triangle fan, just creeping on tops, unfeasible path")



# # Plot triangle fan with points on the tops

# mu1 = 0.7
# lam2 = 0.8
# mu2 = 0.1
# lam3 = 0.1
# mu3 = 0.4
# lam4 = 0.9
# params0 = [mu1, lam2, mu2, lam3, mu3, lam4]

# indCrTop = [1]
# paramsCrTop = [0.3, 0.5]
# indStTop = [3]
# paramsStTop = [0.3, 0.7]

# plotFann(x0, listB0k, listxk, listBk, listBkBk1 = listBkBk1, params = params0,
#          indCrTop = indCrTop, paramsCrTop = paramsCrTop,
#          indStTop = indStTop, paramsStTop = paramsStTop)
# plt.title("Generalized triangle fan, creeping and going through tops, unfeasible path")


