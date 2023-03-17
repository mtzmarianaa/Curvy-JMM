# Before going into C with the optimization rutine
# we test the projected coordinate gradient descent here
# to make sure it works or it makes sense to try this approach

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi
import intermediateTests as itt
from scipy.optimize import root_scalar # to project back blocks [mu_k, lambda_k+1]

def arclengthSimpson(mu, lam, xFrom, Bfrom, xTo, Bto):
     '''
     arclength along a boundary from xLam to xMu
     '''
     Bmu = itt.gradientBoundary(mu, xFrom, Bfrom, xTo, Bto)
     Blam = itt.gradientBoundary(lam, xFrom, Bfrom, xTo, Bto)
     B_mid = itt.gradientBoundary((mu + lam)/2, xFrom, Bfrom, xTo, Bto)
     return (norm(Bmu) + 4*norm(B_mid) + norm(Blam))*(abs(mu - lam)/6)

def hermite_interpolationT(param, x0, T0, grad0, x1, T1, grad1):
    '''
    Hermite interpolation of the eikonal
    '''
    sumGrads = (param**3 - 2*param**2 + param)*grad0 + (param**3 - param**2)*grad1
    return (2*param**3 - 3*param**2 + 1)*T0 + (-2*param**3 + 3*param**2)*T1 + np.dot(x1 - x0, sumGrads)

def der_hermite_interpolationT(param, x0, T0, grad0, x1, T1, grad1):
    '''
    derivative with respecto to param of the Hermite interpolation of the eikonal
    '''
    sumGrads = (3*param**2 - 4*param + 1)*grad0 + (3*param**2 - 2*param)*grad1
    return (6*param**2 - 6*param)*T0 + (-6*param**2 + 6*param)*T1 + np.dot(x1 - x0, sumGrads)

def t1(lam, x0, xk1, B0k1, Bk1, zk, Bk_mu):
    '''
    This function is useful to solve for lamMin
    '''
    yk1 = itt.hermite_boundary(lam, x0, B0k1, xk1, Bk1)
    return Bk_mu[0]*(yk1[1] - zk[1]) - Bk_mu[1]*(yk1[0] - zk[0])

def t2(lam, x0, xk1, B0k1, Bk1, zk):
    '''
    This function is useful to solve for lamMax
    '''
    yk1 = itt.hermite_boundary(lam, x0, B0k1, xk1, Bk1)
    Bk_lam = itt.gradientBoundary(lam , x0, B0k1, xk1, Bk1)
    return Bk_lam[0]*(yk1[1] - zk[1]) - Bk_lam[1]*(yk1[0] - zk[0])

# Find a root finding method

def fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
    '''
    Objective function of an update without the tops on a triangle fan
    params = [mu1, lam2, mu2, lam3, mu3, ..., lambda_n, mu_n, lambda_n1] # length 2n
    listIndices = [eta1, eta2, ..., eta_n, eta_n1]       # length n+1
    listxk = [x0, x1, x2, ..., xn, xn1]                  # length n+2
    listB0k = [B01, B02, B03, ..., B0n, B0n1]            # length n+1
    listBk = [B1, B2, B3, ..., Bn, Bn1]                  # length n+1
    '''
    n = len(listxk) - 2
    muk = params[0]
    etak = listIndices[0]
    Bk = listBk[0]
    B0k = listB0k[0]
    zk = itt.hermite_boundary(muk, x0, B0k, x1, Bk)
    sum = hermite_interpolationT(muk, x0, T0, grad0, x1, T1, grad1)
    for i in range(1, n):
        k = 2*i -1 # starts in k = 1, all the way to k = 2n-3
        mukPrev = muk
        muk = params[k+1] # starts in mu2 ends in mu_n
        lamk = params[k] # starts in lam2 ends in lamn
        etaPrev = etak
        etak = listIndices[i]
        etaMin = min(etaPrev, etak)
        Bk = listBk[i]
        B0k = listB0k[i]
        xk = listxk[i+1]
        zkPrev = zk
        zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
        yk = itt.hermite_boundary(lamk, x0, B0k, xk, Bk)
        sum += etaPrev*norm(yk - zkPrev) + etaMin*arclengthSimpson(muk, lamk, x0, B0k, xk, Bk)
    # now we need to add the last segment
    mukPrev = muk
    lamk = params[2*n-1] # last one
    etaPrev = etak
    etak = listIndices[n]
    etaMin = min(etak, etaPrev)
    Bk = listBk[n]
    B0k = listB0k[n]
    xk = listxk[n+1]
    zkPrev = zk
    yk = itt.hermite_boundary(lamk, x0, B0k, xk, Bk)
    sum += etaPrev*norm(yk - zkPrev) + etaMin*arclengthSimpson(lamk, 1, x0, B0k, xk, Bk)
    return sum

##########
## THESE ARE THE AUXILIARY FUNCTIONS FOR THE BLOCK COORDINATE PROJECTED GRADIENT DESCENT

def partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B01_mu, y2, z1):
    der_hermite_inter = der_hermite_interpolationT(mu1, x0, T0, grad0, x1, T1, grad1)
    return der_hermite_inter - np.dot(B01_mu, y2 - z1)/norm(y2 - z1)

def partial_fObj_mu(muk, etakM1, B0k_muk, yk, zkM1, etaMin):
    return etakM1*np.dot(-B0k_muk, yk - zkM1)/norm(yk - zkM1) + etaMin*norm(B0k_muk)

def partial_fObj_lambda(lambdak, etakM1, B0k_lamk, yk, zkM1, etaMin):
    return etakM1*np.dot(B0k_lamk, yk - zkM1)/norm(yk - zkM1) - etaMin*norm(B0k_lamk)


def project_block(muk, lamk1, Bk_muk, Bk1_lamk1, yk1, zk, x0, xk1, B0k1, Bk1):
    '''
    Project a block [mu_k, lambda_k+1] such that it is feasible
    '''
    Nk_muk = np.array([-Bk_muk[1], Bk_muk[0]])
    Nk1_lamk1 = np.array([-Bk1_lamk1[1], Bk1_lamk1[0]])
    dotTestMin = np.dot( yk1 - zk, Nk_muk )
    dotTestMax = np.dot( yk1 - zk, Nk1_lamk1 )
    print("       dotTestMin: ", dotTestMin, "  dotTestMax: ", dotTestMax)
    # Test if lamk < lamMin
    if( dotTestMin < 0):
        # Means that we need to find lambdaMin
        tMin = lambda lam: t1(lam, x0, xk1, B0k1, Bk1, zk, Bk_muk)
        # Find the root of tMin
        rootMin = root_scalar(tMin, bracket=[0, 1])
        lamk1 = rootMin.root
        print("       lambda < lambdaMin")
    if( dotTestMax < 0):
        # Means that we need to find lambdaMax
        tMax = lambda lam: t2(lam, x0, xk1, B0k1, Bk1, zk)
        rootMax = root_scalar(tMax, bracket=[0, 1])
        lamk1 = rootMax.root
        print("       lambda > lambdaMax")
    if( muk < 0):
        muk = 0
    elif( muk > 1):
        muk = 1
    if( lamk1 < 0):
        lamk1 = 0
    elif( lamk1 > 1):
        lamk1 = 1
    return muk, lamk1

def backTr_block(alpha0, k, dmuk, dlamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
    '''
    Backtracking to find the step size for block coordinate gradint descent
    '''
    i = 0
    alpha = alpha0
    f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    params_test = params
    params_test[k] = params[k] - alpha*dmuk
    params_test[k+1] = params[k+1] - alpha*dlamk1
    f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    while( f_before <= f_after and i < 25):
        alpha = alpha*0.5
        params_test[k] = params[k] - alpha*dmuk
        params_test[k+1] = params[k+1] - alpha*dlamk1
        f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        i += 1
    return alpha

def forwardPassUpdate(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
    '''
    Forward pass for the block coordinate gradient descent goes from k = 0 to k = 2n-1
    '''
    print("  At block iteration ", -1, "  muk = ", params[0], "   lamk1 = ", params[1])
    print("  At block iteration ", -1, " params:", params)
    itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params)
    muk = params[0]
    B0k_muk = itt.gradientBoundary(muk, x0, listB0k[0], listxk[1], listBk[0])
    zk = itt.hermite_boundary(muk, x0, listB0k[0], x1, listBk[0])
    lamk1 = params[1]
    yk1 = itt.hermite_boundary(lamk1, x0, listB0k[1], listxk[2], listBk[1])
    B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, listB0k[1], listxk[2], listBk[1])
    etaMin = min(listIndices[0], listIndices[1])
    # Compute direction
    partial_muk = partial_fObj_mu1(muk, x0, T0, grad0, x1, T1, grad1, B0k_muk, yk1, zk)
    partial_lamk1 = partial_fObj_lambda(lamk1, listIndices[0], B0k1_lamk1, yk1, zk, etaMin)
    # Compute step size
    alpha = backTr_block(1, 0, partial_muk, partial_lamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    # Update this block
    muk = muk - alpha*partial_muk
    lamk1 = lamk1 - alpha*partial_lamk1
    # Project back so that it is feasible
    muk, lamk1 = project_block(muk, lamk1, B0k_muk, B0k1_lamk1, yk1, zk, x0, listxk[2], listB0k[1], listBk[1])
    params[0] = muk
    params[1] = lamk1
    print("  At block iteration ", 0, "  muk = ", params[0], "   lamk1 = ", params[1])
    print("  At block iteration ", 0, " params:", params)
    itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params)
    # Start the forward loop
    n = len(listxk) - 2
    for j in range(1, n):
        k = 2*j
        muk = params[k]
        B0k_muk = itt.gradientBoundary(muk, x0, listB0k[j], listxk[j+1], listBk[j])
        zk = itt.hermite_boundary(muk, x0, listB0k[j], x1, listBk[j])
        lamk1 = params[k+1]
        yk1 = itt.hermite_boundary(lamk1, x0, listB0k[j+1], listxk[j+2], listBk[j+1])
        B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, listB0k[j+1], listxk[j+2], listBk[j+1])
        etaMin = min(listIndices[j], listIndices[j+1])
        # Compute direction
        partial_muk = partial_fObj_mu(muk, listIndices[j], B0k_muk, yk1, zk, etaMin)
        partial_lamk1 = partial_fObj_lambda(lamk1, listIndices[j], B0k1_lamk1, yk1, zk, etaMin)
        # Compute step size
        alpha = backTr_block(1, k, partial_muk, partial_lamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        # Update this block
        muk = muk - alpha*partial_muk
        lamk1 = lamk1 - alpha*partial_lamk1
        # Project back so that it is feasible
        muk, lamk1 = project_block(muk, lamk1, B0k_muk, B0k1_lamk1, yk1, zk, x0, listxk[j+2], listB0k[j+1], listBk[j+1])
        params[k] = muk
        params[k+1] = lamk1
        print("  At block iteration ", j, "  muk = ", params[k], "   lamk1 = ", params[k+1])
        print("  At block iteration ", j, " params:", params)
        itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params)
    return params
        
    

    
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

# params = [mu1, lam2, mu2, lam3, mu3, lam4, mu4, ..., lambda_n, mu_n, lambda_n1]
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


# Starting to plot everything

itt.plotFan3(x0, B01, B02, B03, B0Hat, x1, B1, x2, B2, x3, B3, xHat, BHat, mu1 = mu1, lam2 = lam2, mu2 = mu2, lam3 = lam3, mu3 = mu3, lam4 = lam4)

plt.scatter(xSource[0], xSource[1], s = 20, c = "#ff00f2")
plt.plot([xSource[0], xHat[0]], [xSource[1], xHat[1]], linestyle = ':', c = "#ab00a3")


f_init = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)


##### Test the projection

lam2 = 0
mu1 = 1
Bk_muk = itt.gradientBoundary(mu1, x0, B01, x1, B1)
Bk1_lamk1 = itt.gradientBoundary(lam2, x0, B02, x2, B2)
yk1 = itt.hermite_boundary(lam2, x0, B02, x2, B2)
zk = itt.hermite_boundary(mu1, x0, B01, x1, B1)

itt.plotFan3(x0, B01, B02, B03, B0Hat, x1, B1, x2, B2, x3, B3, xHat, BHat, mu1 = mu1, lam2 = lam2, mu2 = mu2, lam3 = lam3, mu3 = mu3, lam4 = lam4)

mu1, lam2 = project_block(mu1, lam2, Bk_muk, Bk1_lamk1, yk1, zk, x0, x2, B02, B2)

itt.plotFan3(x0, B01, B02, B03, B0Hat, x1, B1, x2, B2, x3, B3, xHat, BHat, mu1 = mu1, lam2 = lam2, mu2 = mu2, lam3 = lam3, mu3 = mu3, lam4 = lam4)


##### Test the forward and backward updates

print("Start test foward pass update \n\n")

mu1 = 1
lam2 = 0
mu2 = 1
lam3 = 0
mu3 = 1
lam4 = 1
params = [mu1, lam2, mu2, lam3, mu3, lam4]

paramsUp = forwardPassUpdate(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

print(paramsUp)


