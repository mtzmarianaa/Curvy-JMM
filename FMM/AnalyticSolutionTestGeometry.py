from scipy.optimize import NonlinearConstraint, minimize
import numpy as np
from numpy.linalg import norm 
from math import sqrt

def whichRegion(xhat):
    '''
    This function determines in which of the 3 regions in the test geometry xhat is
    '''
    if (norm(xhat)< 10):
        reg = 3
    elif(sqrt( xhat[0]**2 + (xhat[1]-5*sqrt(2))**2  )< 5*sqrt(2)):
        reg = 2
    else:
        reg = 1
    return reg

def regA1(xhat):
    '''
    This function determines if xhat can be accessed directly from x0
    '''
    ytan = (26*(7 + sqrt(14))*xhat[0])/(273 - 42*sqrt(13)) + (140*sqrt(13) + 130*sqrt(14))/(91 - 14*sqrt(13))
    if ( xhat[1] >= ytan | xhat[1]<= -10 |  (xhat[0]<=0 & norm(xhat)>= 10 & xhat[1]<= 20/sqrt(13)) ):
        ans = True
    else:
        ans = False
    return ans

def segmentBetweenPoints(lam, a, b):
    '''
    Parametrization of the segment from a to b (using lam)
    '''
    return (1-lam)*a + lam*b

def inCt(x):
    '''
    Function that determines if x is on Ct (top of the snowman)
    '''
    if ( sqrt( x[0]**2 + (x[1]-5*sqrt(2))**2  )== 5*sqrt(2) & norm(x)>=10 ):
        ans = True
    else:
        ans = False
    return ans

def inCb(x):
    '''
    Function that determines if x is on Cb (bottom of the snowman AND NOT THE ARCH BETWEEN reg2 and reg3)
    '''
    if (norm(x)== 10 & x[1]< 5*sqrt(2)):
        ans = True
    else:
        ans = False
    return ans

def inArch(x):
    '''
    Function that determines if x is on the arch that separates reg2 from reg3
    '''
    if(norm(x)==10 & x[1]> 5*sqrt(2)):
        ans = True
    else:
        ans = False
    return ans

def SnowSolution(xhat):
    '''
    This is the analytic solution to the Eikonal equation with the snowman domain starting at [-15, -10]
    '''
    start = np.array(  [-15, -10]  )
    regHat = whichRegion(xhat) # we need to know in which region xhat is
    if(regHat == 1):
        if( regA1(xhat) ):
            eik = norm( np.subtract(start, xhat) )
        else:
            t1 =