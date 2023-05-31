#minimization

import numpy as np
from numpy.linalg import norm
from scipy.optimize import minimize
from math import sqrt

def fToMin(param):
    lam = param[0]
    mu = param[1]
    x0 = np.array([-2, 1])
    x1 = np.array( [0, 1] )
    x2 = np.array( [0, 2] )
    xHat = np.array( [-1, 3] )
    ind1 = 1.0
    ind2 = 1.453
    T0 = sqrt(2)
    T1 = sqrt(2)
    xmu = np.add( (1-mu)*x0, mu*x2 )
    xlam = np.add( (1-lam)*x0, lam*x1  )
    return lam*(T1-T0) + T0 + ind1*norm( np.subtract(xmu, xlam)  ) + ind2*norm( np.subtract(xHat, xmu)  )

x0 = [1,1]
res = minimize(fToMin, x0, bounds=( (0,1), (0,1) )  )
print(res.x)
print(res)
