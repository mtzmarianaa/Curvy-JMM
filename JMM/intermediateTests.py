import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi

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
 


 
