
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


####### Projections

## Finding lamk1Min and lamk1Max given muk

def t1(lam, x0, xk1, B0k1, Bk1, zk, Bk_mu):
    '''
    This function is useful to solve for lamk1Min
    given muk
    '''
    yk1 = itt.hermite_boundary(lam, x0, B0k1, xk1, Bk1)
    return Bk_mu[0]*(yk1[1] - zk[1]) - Bk_mu[1]*(yk1[0] - zk[0])

def t2(lam, x0, xk1, B0k1, Bk1, zk):
    '''
    This function is useful to solve for lamk1Max
    given muk
    '''
    yk1 = itt.hermite_boundary(lam, x0, B0k1, xk1, Bk1)
    Bk_lam = itt.gradientBoundary(lam , x0, B0k1, xk1, Bk1)
    return Bk_lam[0]*(yk1[1] - zk[1]) - Bk_lam[1]*(yk1[0] - zk[0])

## Find mukMin and mukMAx given lamk1

def t3(mu, x0, xk, B0k, Bk, yk1):
     '''
     This function is useful to solve for mukMax
     given lamk1
     '''
     zk = itt.hermite_boundary(mu, x0, B0k, xk, Bk)
     Bk_mu = itt.gradientBoundary(mu, x0, B0k, xk, Bk)
     return Bk_mu[0]*(yk1[1] - zk[1]) - Bk_mu[1]*(yk1[0] - zk[0])
     

def t4(mu, x0, xk, B0k, Bk, yk1, Bk_lam):
     '''
     This function is useful to solve for mukMin
     given lamk1
     '''
     zk = itt.hermite_boundary(mu, x0, B0k, xk, Bk)
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

def partial_L_muk(muk, lamk, B0k_muk, secondDer_B0k_muk, B0k_halves, secondDer_B0khalves_muk, B0k_lamk):
     '''
     partial of the approximation of the arc length with respect to muk
     '''
     if( abs(muk - lamk) <= 1e-14):
          return 0
     else:
          normB0k_muk = norm(B0k_muk)
          normB0k_halves = norm(B0k_halves)
          sk = get_sk(muk, lamk)
          firstPart = sk/6*(normB0k_muk + 4*normB0k_halves + norm(B0k_lamk) )
          secondPart = (abs(muk - lamk)/6)*(np.dot(secondDer_B0k_muk, B0k_muk)/normB0k_muk + 2*np.dot(secondDer_B0khalves_muk, B0k_halves)/normB0k_halves)
          return firstPart + secondPart

def partial_L_lamk(muk, lamk, B0k_muk, B0k_halves, secondDer_B0khalves_lamk, B0k_lamk, secondDer_B0k_lamk):
     '''
     partial of the approximation of the arc length with respect to lamk
     '''
     if( abs(muk - lamk) <= 1e-14):
          return 0
     else:
          normB0k_halves = norm(B0k_halves)
          normB0k_lamk = norm(B0k_lamk)
          sk = get_sk(muk, lamk)
          firstPart = -sk/6*(norm(B0k_muk) + 4*normB0k_halves + normB0k_lamk )
          secondPart = (abs(muk - lamk)/6)*(2*np.dot(secondDer_B0khalves_lamk, B0k_halves)/normB0k_halves + np.dot(secondDer_B0k_lamk, B0k_lamk)/normB0k_lamk )
     return firstPart + secondPart
     

def partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B01_mu, y2, z1):
    der_hermite_inter = der_hermite_interpolationT(mu1, x0, T0, grad0, x1, T1, grad1)
    return der_hermite_inter - np.dot(B01_mu, y2 - z1)/norm(y2 - z1)

def partial_fObj_muk(muk, lamk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etakM1):
     '''
     Partial of the objective function (triangle fan no tops) with respect to muk
     '''
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
     secondDer_B0k_muk = itt.secondDer_Boundary(muk, x0, B0k, xk, Bk)
     B0k_halves = itt.gradientBoundary((muk + lamk)/2, x0, B0k, xk, Bk)
     secondDer_B0khalves_muk = itt.secondDer_Boundary((muk + lamk)/2, x0, B0k, xk, Bk)
     B0k_lamk = itt.gradientBoundary(lamk, x0, B0k, xk, Bk)
     parL_muk = partial_L_muk(muk, lamk, B0k_muk, secondDer_B0k_muk, B0k_halves, secondDer_B0khalves_muk, B0k_lamk)
     etaMin = min(etak, etakM1)
     return etakM1*np.dot(-B0k_muk, yk1 - zk)/norm(yk1 - zk) + etaMin*parL_muk


def partial_fObj_lambdak1(muk, muk1, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etak1):
     '''
     Partial of the objective function (triangle fan no tops) with respect to lamk1
     '''
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_muk1 = itt.gradientBoundary(muk1, x0, B0k1, xk1, Bk1)
     secondDer_B0k1_lamk1 = itt.secondDer_Boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_halves = itt.gradientBoundary( (muk1 + lamk1)/2, x0, B0k1, xk1, Bk1)
     secondDer_B0k1halves_lamk1 = itt.secondDer_Boundary((muk1 + lamk1)/2, x0, B0k, xk, Bk)
     B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
     perL_lamk1 = partial_L_lamk(muk1, lamk1, B0k1_muk1, B0k1_halves, secondDer_B0k1halves_lamk1, B0k1_lamk1, secondDer_B0k1_lamk1)
     etaMin = min(etak, etak1)
     return etak*np.dot(B0k1_lamk1, yk1 - zk)/norm(yk1 - zk) + etaMin*perL_lamk1

def partial_fObj_mu(muk, etakM1, B0k_muk, yk1, zk, etaMin, sk):
    return etakM1*np.dot(-B0k_muk, yk1 - zk)/norm(yk1 - zk) + sk*etaMin*norm(B0k_muk)

def partial_fObj_lambda(lambdak, etakM1, B0k_lamk, yk, zkM1, etaMin, sk):
    return etakM1*np.dot(B0k_lamk, yk - zkM1)/norm(yk - zkM1) - sk*etaMin*norm(B0k_lamk)


def project_lamk1Givenmuk(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1):
    '''
    Project back lamk1 given muk
    '''
    zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
    yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
    B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
    B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
    
    Nk_muk = np.array([-B0k_muk[1], B0k_muk[0]])
    Nk1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
    dotTestMin = np.dot( yk1 - zk, Nk_muk )
    dotTestMax = np.dot( yk1 - zk, Nk1_lamk1 )
    #print("       dotTestMin: ", dotTestMin, "  dotTestMax: ", dotTestMax)
    # Test if lamk < lamMin
    if( dotTestMin < 0):
        # Means that we need to find lambdaMin
        tMin = lambda lam: t1(lam, x0, xk1, B0k1, Bk1, zk, B0k_muk)
        # Find the root of tMin
        rootMin = root_scalar(tMin, bracket=[0, 1])
        lamk1 = rootMin.root
        #print("       lambda < lambdaMin")
    if( dotTestMax < 0):
        # Means that we need to find lambdaMax
        tMax = lambda lam: t2(lam, x0, xk1, B0k1, Bk1, zk)
        rootMax = root_scalar(tMax, bracket=[0, 1])
        lamk1 = rootMax.root
        #print("       lambda > lambdaMax")
    lamk1 = project_box(lamk1) # Such that 0<=lamk1 <=1
    return lamk1

def project_mukGivenlamk1(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1):
     '''
     Project back muk given lamk1
     '''
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
     # Compute the normals
     N0k1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
     N0k_muk = np.array([-B0k_muk[1], B0k_muk[1]])
     dotTestMin =  np.dot(N0k1_lamk1, yk1 - zk)
     dotTestMax = np.dot(N0k_muk, yk1 - zk)
     if(dotTestMin<0):
          print("  failed dotTestMin")
          print("  zk: ", zk, "  yk1: ", yk1)
          tMin = lambda mu: t4(mu, x0, xk, B0k, Bk, yk1, B0k1_lamk1)
          rootMin = root_scalar(tMin, bracket = [0,1])
          muk = rootMin.root
     if(dotTestMax<0):
          print("  failed dotTestMax")
          print("  zk: ", zk, "  yk1: ", yk1)
          tMax = lambda mu: t3(mu, x0, xk, B0k, Bk, yk1)
          rootMax = root_scalar(tMax, bracket = [0,1])
          muk = rootMax.root
     muk = project_box(muk)
     return muk

def project_blockBackwards(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1):
     '''
     Project a block [muk, lamk1] such that it is feasible
     project back muk given lamk1
     '''
     lamk1 = project_box(lamk1)
     # Fiven lamk1 project back muk
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
     # Compute the normals
     N0k1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
     N0k_muk = np.array([-B0k_muk[1], B0k_muk[1]])
     dotTestMin =  np.dot(N0k1_lamk1, yk1 - zk)
     dotTestMax = np.dot(N0k_muk, yk1 - zk)
     if(dotTestMin<0):
          print("  failed dotTestMin")
          print("  zk: ", zk, "  yk1: ", yk1)
          tMin = lambda mu: t4(mu, x0, xk, B0k, Bk, yk1, B0k1_lamk1)
          rootMin = root_scalar(tMin, bracket = [0,1])
          muk = rootMin.root
     if(dotTestMax<0):
          print("  failed dotTestMax")
          print("  zk: ", zk, "  yk1: ", yk1)
          tMax = lambda mu: t3(mu, x0, xk, B0k, Bk, yk1)
          rootMax = root_scalar(tMax, bracket = [0,1])
          muk = rootMax.root
     muk = project_box(muk)
     return muk, lamk1
     
def project_ontoLine(d):
     '''
     Projection of a vector d along the line x=y
     '''
     projector = 0.5*np.array([ [1,1], [1,1]])
     return projector@d

def close_to_identity(lamk, muk):
     '''
     Computes the distance of the vector [lamk, muk]
     to the line lamk = muk
     '''
     p = np.array([lamk, muk])
     proj = project_ontoLine(p)
     return norm(p - proj)

def backTrClose_firstBlock(alpha0, dmu1, dlam2, dmu2, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [mu1, lam2, mu2]
     it considers two directions: steepest descent
                                  steepest descent projected onto the line lam2 = mu2
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     alpha = alpha0
     d_middle = np.array([dlamk2, dmu2]) # "normal" steespest descent
     d_middle_proj = project_ontoLine(d_middle) # steepest descent projected onto the line lam2 = mu2
     params_test = np.copy(params)
     params_test[0] = params[0] - alpha*dmu1
     params_test[1:3] = params[1:3] - alpha*d_middle
     params_test_proj = np.copy(params)
     params_test_proj[0] = params[0] - alpha*dmu1
     params_test_proj[1:3] = params[1:3] - alpha*d_middle_proj
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_test and f_before <= f_test_proj and i < 25):
          alpha = alpha*0.75
          params_test[0] = params[0] - alpha*dmu1
          params_test[1:3] = params[1:3] - alpha*d_middle
          params_test_proj[0] = params[0] - alpha*dmu1
          params_test_proj[1:3] = params[1:3] - alpha*d_middle_proj
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test and f_before <= f_test_proj ):
          paramsOpt = params
     elif( f_test < f_test_proj):
          paramsOpt = f_test
     elif( f_test_proj <= f_test ):
          paramsOpt = f_test_proj
     return paramsOpt

def backTr_firstBlock(alpha0, dmu1, dlam2, dmu2, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [mu1, lam2, mu2]
     it considers one direction:  steepest descent
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     alpha = alpha0
     d_middle = np.array([dlamk2, dmu2]) # "normal" steespest descent
     params_test = np.copy(params)
     params_test[0] = params[0] - alpha*dmu1
     params_test[1:3] = params[1:3] - alpha*d_middle
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_test and i < 25):
          alpha = alpha*0.75
          params_test[0] = params[0] - alpha*dmu1
          params_test[1:3] = params[1:3] - alpha*d_middle
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test ):
          paramsOpt = params
     else:
          paramsOpt = f_test
     return paramsOpt


def backTrClose_block(alpha0, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [lamk, muk]
     it considers two directions: steepest descent
                                  steepest descent projected onto the line lamk = muk
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     alpha = alpha0
     d_middle = np.array([dlamk, dmuk]) # "normal" steespest descent
     d_middle_proj = project_ontoLine(d_middle) # steepest descent projected onto the line lamk = muk
     params_test = np.copy(params)
     params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
     params_test_proj = np.copy(params)
     params_test_proj[k:(k+2)] = params[k:(k+2)] - alpha*d_middle_proj
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_test and f_before <= f_test_proj and i < 25):
          alpha = alpha*0.75
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          params_test_proj[k:(k+2)] = params[k:(k+2)] - alpha*d_middle_proj
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test and f_before <= f_test_proj ):
          paramsOpt = params
     elif( f_test < f_test_proj):
          paramsOpt = f_test
     elif( f_test_proj <= f_test ):
          paramsOpt = f_test_proj
     return paramsOpt

def backTr_block(alpha0, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [lamk, muk]
     it considers one direction:  steepest descent
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     alpha = alpha0
     d_middle = np.array([dlamk, dmuk]) # "normal" steespest descent
     params_test = np.copy(params)
     params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_test and f_before <= f_test_proj and i < 25):
          alpha = alpha*0.75
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test ):
          paramsOpt = params
     else:
          paramsOpt = f_test
     return paramsOpt

def backtrClose_lastBlock(alpha0, nRegions, dlamn, dmun, dlamn1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [lamn, mun, lamn1]
     it considers two directions: steepest descent
                                  steepest descent projected onto the line lam2 = mu2
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     alpha = alpha0
     d_middle = np.array([dlamn, dmun]) # "normal" steepest descent
     d_middle_proj = project_ontoLine(d_middle) # steepest descent projected onto the line lamn = mun
     params_test = np.copy(params)
     params_test[(2*nRegions - 3):(2*nRegions - 1)] = params[(2*nRegions - 3):(2*nRegions - 1)] - alpha*d_middle
     params_test[2*nRegions -1] = params[2*nRegions - 1] - alpha*dlamn1
     params_test_proj = np.copy(params)
     params_test_proj[(2*nRegions - 3):(2*nRegions - 1)] = params[(2*nRegions - 3):(2*nRegions - 1)] - alpha*d_middle_proj
     params_test_proj[2*nRegions -1] = params[2*nRegions - 1] - alpha*dlamn1
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_test and f_before <= f_test_proj and i < 25):
          alpha = alpha*0.75
          params_test[(2*nRegions - 3):(2*nRegions - 1)] = params[(2*nRegions - 3):(2*nRegions - 1)] - alpha*d_middle
          params_test[2*nRegions -1] = params[2*nRegions - 1] - alpha*dlamn1
          params_test_proj[(2*nRegions - 3):(2*nRegions - 1)] = params[(2*nRegions - 3):(2*nRegions - 1)] - alpha*d_middle_proj
          params_test_proj[2*nRegions -1] = params[2*nRegions - 1] - alpha*dlamn1
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test and f_before <= f_test_proj ):
          paramsOpt = params
     elif( f_test < f_test_proj):
          paramsOpt = f_test
     elif( f_test_proj <= f_test ):
          paramsOpt = f_test_proj
     return paramsOpt

def backtr_lastBlock(alpha0, nRegions, dlamn, dmun, dlamn1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [lamn, mun, lamn1]
     it considers one direction:  steepest descent
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     alpha = alpha0
     d_middle = np.array([dlamn, dmun]) # "normal" steepest descent
     params_test = np.copy(params)
     params_test[(2*nRegions - 3):(2*nRegions - 1)] = params[(2*nRegions - 3):(2*nRegions - 1)] - alpha*d_middle
     params_test[2*nRegions -1] = params[2*nRegions - 1] - alpha*dlamn1
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_test and f_before <= f_test_proj and i < 25):
          alpha = alpha*0.75
          params_test[(2*nRegions - 3):(2*nRegions - 1)] = params[(2*nRegions - 3):(2*nRegions - 1)] - alpha*d_middle
          params_test[2*nRegions -1] = params[2*nRegions - 1] - alpha*dlamn1
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test ):
          paramsOpt = params
     else:
          paramsOpt = f_test
     return paramsOpt
     

def project_block(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1):
    '''
    Project a block [mu_k, lambda_k+1] such that it is feasible
    project back lamk1 given muk
    '''
    zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
    yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
    B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
    B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
    
    Nk_muk = np.array([-B0k_muk[1], B0k_muk[0]])
    Nk1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
    dotTestMin = np.dot( yk1 - zk, Nk_muk )
    dotTestMax = np.dot( yk1 - zk, Nk1_lamk1 )
    #print("       dotTestMin: ", dotTestMin, "  dotTestMax: ", dotTestMax)
    # Test if lamk < lamMin
    if( dotTestMin < 0):
        # Means that we need to find lambdaMin
        tMin = lambda lam: t1(lam, x0, xk1, B0k1, Bk1, zk, B0k_muk)
        # Find the root of tMin
        rootMin = root_scalar(tMin, bracket=[0, 1])
        lamk1 = rootMin.root
        #print("       lambda < lambdaMin")
    if( dotTestMax < 0):
        # Means that we need to find lambdaMax
        tMax = lambda lam: t2(lam, x0, xk1, B0k1, Bk1, zk)
        rootMax = root_scalar(tMax, bracket=[0, 1])
        lamk1 = rootMax.root
        #print("       lambda > lambdaMax")
    if( muk < 0):
        muk = 0
    elif( muk > 1):
        muk = 1
    if( lamk1 < 0):
        lamk1 = 0
    elif( lamk1 > 1):
        lamk1 = 1
    return muk, lamk1


def project_middleBlock(k, mukM1, lamk, muk, lamk1, x0, listB0k, listxk, listBk):
     '''
     Projection so that when we move the block [lamk, muk] the blocks
     [mukM1, lamk], [muk, lamk1] are feasible
     '''
     B0kM1 = listB0k[k-1]
     xkM1 = listxk[k]
     BkM1 = listBk[k-1]
     B0k = listB0k[k]
     xk = listxk[k+1]
     Bk = listBk[k]
     B0k1 = listB0k[k+1]
     xk1 = listxk[k+2]
     Bk1 = listBk[k+1]
     mukM1, lamk = project_block(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
     muk, lamk1 = project_block(muk, lamk1, x0, B0k, xk, B0k, B0k1, xk1, Bk1)
     return mukM1, lamk, muk, lamk1

def project_4Block0(mu1, lam2, mu2, lam3, x0, B01, x1, B1, B02, x2, B2, B03, x3, B3):
     '''
     Project the block [mu1, lam2, mu2] back to the feasible set
     (mainly: given mu1, project back lam2, given lam3 project back mu2)
     '''
     mu1, lam2 = project_block(mu1, lam2, x0, B01, x1, B1, B02, x2, B2)
     # Then we need to project back mu2 given lam3
     y3 = itt.hermite_boundary(lam3, x0, B03, x3, B3)
     B03_lam3 = itt.gradientBoundary(lam3, x0, B03, x3, B3)
     z2 = itt.hermite_boundary(mu2, x0, B02, x2, B2)
     B02_mu2 = itt.gradientBoundary(mu2, x0, B02, x2, B2)
     # Compute the normals
     N03_lam3 = np.array([-B03_lam3[1], B03_lam3[0]])
     N02_mu2 = np.array([-B02_mu2[1], B02_mu2[0]])
     dotTestMin = np.dot(N03_lam3, y3 - z2)
     dotTestMax = np.dot(N02_mu2, y3 - z2)
     # Test if mu2 < mu2Min
     if( dotTestMin < 0):
          # Means that we need to find mu2Min
          print(" dotTestMin -- 4 block0")
          tMin = lambda mu: t4(mu, x0, x2, B02, B2, y3, B03_lam3)
          # Find the root of tMin
          rootMin = root_scalar(tMin, bracket = [0,1])
          mu2 = rootMin.root
     if( dotTestMax < 0):
          # Means that we need to find mu2Max
          print(" dotTestMax -- 4 block0")
          tMax = lambda mu: t3(mu, x0, x2, B02, B2, y3)
          rootMax = root_scalar(tMax, bracket = [0,1])
          mu2 = rootMax.root
     if( mu2 < 0):
          mu2 = 0
     elif( mu2 > 1):
          mu2 = 1
     return mu1, lam2, mu2

def project_4Block(mukM1, lamk, muk, lamk1, muk1, lamk2, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, B0k1, xk1, Bk1, B0k2, xk2, Bk2):
     '''
     Project block [lamk, muk, lamk1, muk1] back to the feasible set
     mainly:  given mukM1 project back lamk,
              given lamk2 project back muk1,
              given muk project back lamk1
     '''
     # First step: given mukM1 project back lamk
     mukM1, lamk = project_block(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
     # Second step: given lamk2 project back muk1
     yk2 = itt.hermite_boundary(lamk2, x0, B0k2, xk2, Bk2)
     B0k2_lamk2 = itt.gradientBoundary(lamk2, x0, B0k2, xk2, Bk2)
     zk1 = itt.hermite_boundary(muk1, x0, B0k1, xk1, Bk1)
     B0k1_muk1 = itt.gradientBoundary(muk1, x0, B0k1, xk1, Bk1)
     #Compute the normals
     N0k2_lamk2 = np.array([-B0k2_lamk2[1], B0k2_lamk2[0]])
     N0k1_muk1 = np.array([-B0k1_muk1[1], B0k1_muk1[0]])
     dotTestMin = np.dot(N0k2_lamk2, yk2 - zk1)
     dotTestMax = np.dot(N0k1_muk1, yk2 - zk1)
     #Test if muk1 < muk1Min
     if(dotTestMin<0):
          print(" dotTestMin -- 4 block")
          tMin = lambda mu: t4(mu, x0, xk1, B0k1, Bk1, yk2, B0k2_lamk2)
          rootMin = root_scalar(tMin, bracket = [0,1])
          muk1 = rootMin.root
     if(dotTestMax<0):
          print(" dotTestMax -- 4 block")
          tMax = lambda mu: t3(mu, x0, xk1, B0k1, Bk1, yk2)
          rootMax = root_scalar(tMax, bracket = [0,1])
          muk1 = rootMax.root
     muk1 = project_box(muk1)
     # Third step: given muk project back lamk1
     muk, lamk1 = project_block(muk, lamk1, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
     return lamk, muk, lamk1, muk1

def project_4Blockn(munM1, lamn, mun, lamn1, x0, B0nM1, xnM1, BnM1, B0n, xn, Bn, B0n1, xn1, Bn1):
     '''
     Project block [lamn, mun, lamn1] back to the feasible set
     mainly:  given munM1 project back lamn
              given lamn1 project back mun
     '''
     lamn1 = project_box(lamn1)
     # First step: given munM1 project back lamn
     munM1, lamn = project_block(munM1, lamn, x0, B0nM1, xnM1, BnM1, B0n, xn, Bn)
     # Second step: given lamn1 project back mun
     yn1 = itt.hermite_boundary(lamn1, x0, B0n1, xn1, Bn1)
     B0n1_lamn1 = itt.gradientBoundary(lamn1, x0, B0n1, xn1, Bn1)
     zn = itt.hermite_boundary(mun, x0, B0n, xn, Bn)
     B0n_mun = itt.gradientBoundary(mun, x0, B0n, xn, Bn)
     # Compute the normals
     N0n1_lamn1 = np.array([-B0n1_lamn1[1], B0n1_lamn1[0]])
     N0n_mun = np.array([-B0n_mun[1], B0n_mun[0]])
     dotTestMin = np.dot(N0n1_lamn1, yn1 - zn)
     dotTestMax = np.dot(N0n_mun, yn1 - zn)
     if(dotTestMin<0):
          print(" dotTestMin -- 4 block n")
          print("B0n: ", B0n, "  xn: ", xn, "  Bn: ", Bn)
          print("  mun: ", mun)
          print("zn: ", zn, "   yn1:", yn1)
          print(" t4 at mun: ", t4(mun, x0, xn, B0n, Bn, yn1, B0n1_lamn1))
          tMin = lambda mu: t4(mu, x0, xn, B0n, Bn, yn1, B0n1_lamn1)
          rootMin = root_scalar(tMin, bracket = [0,1])
          mun = rootMin.root
     if(dotTestMax<0):
          print(" dotTestMax -- 4 block n")
          print("zn: ", zn, "   yn1:", yn1)
          print("B0n: ", B0n, "  xn: ", xn, "  Bn: ", Bn)
          print("  mun: ", mun)
          print(" t3 at mun: ", t3(mun, x0, xn, B0n, Bn, yn1))
          tMax = lambda mu: t3(mu, x0, xn, B0n, Bn, yn1)
          rootMax = root_scalar(tMax, bracket = [0,1])
          mun = rootMax.root
     mun = project_box(mun)
     return lamn, mun, lamn1
     
     

def backTr_block(alpha0, k, dmuk, dlamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the step size for block coordinate gradint descent
     '''
     i = 0
     alpha = alpha0
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     #print("     params before backTr: ", params)
     params_test = np.copy(params)
     params_test[k] = params[k] - alpha*dmuk
     params_test[k+1] = params[k+1] - alpha*dlamk1
     #print("     test params backTr: ", params_test)
     f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     #print("         Params init:         ", params_test)
     while( f_before <= f_after and i < 25):
         alpha = alpha*0.5
         params_test[k] = params[k] - alpha*dmuk
         params_test[k+1] = params[k+1] - alpha*dlamk1
         #print("     test params backTr: ", params_test)
         f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         i += 1
     
     if( f_before <= f_after):
         alpha = 0
     return alpha

def backTr_middleBock(alpha0, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the step size for the middle block coordinate gradient descent
     '''
     i = 0
     alpha = alpha0
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     params_test = np.copy(params)
     params_test[k] = params[k] - alpha*dlamk
     params_test[k+1] = params[k+1] - alpha*dmuk
     f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_after and i < 25):
          alpha = alpha*0.5
          params_test[k] = params[k] - alpha*dlamk
          params_test[k+1] = params[k+1] - alpha*dmuk
          f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_after):
         alpha = 0
     return alpha

def project_box(param):
     if(param < 0):
          param = 0
     elif(param > 1):
          param = 1
     return param

def backTr_coord(alpha0, k, d, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Coordinate backtracking to find the step size for coordinate gradient descent
     '''
     i = 0
     alpha = alpha0
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     params_test = np.copy(params)
     params_test[k] = params[k] - alpha*d
     f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_after and i < 25):
          alpha = alpha*0.5
          params_test[k] = params[k] - alpha*d
          f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_after):
          alpha = 0
     return alpha

def backTr_4Block0(alpha0, dmu1, dlam2, dmu2, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking for the 4 pass update for the first block [mu1, lam2, mu2]
     '''
     i = 0
     alpha = alpha0
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     params_test = np.copy(params)
     params_test[0] = params[0] - alpha*dmu1
     params_test[1] = params[1] - alpha*dlam2
     params_test[2] = params[2] - alpha*dmu2
     f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_before <= f_after and i < 25):
          alpha = alpha*0.5
          params_test[0] = params[0] - alpha*dmu1
          params_test[1] = params[1] - alpha*dlam2
          params_test[2] = params[2] - alpha*dmu2
          f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if (f_before <= f_after):
          alpha = 0
     return alpha


def backTr_4Block(alpha0, k, dlamk, dmuk, dlamk1, dmuk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking for the 4 pass update in the general block [lamk, muk, lamk1, muk1]
     '''
     i = 0
     alpha = alpha0
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     params_test = np.copy(params)
     params_test[k] = params[k] - alpha*dlamk
     params_test[k+1] = params[k+1] - alpha*dmuk
     params_test[k+2] = params[k+2] - alpha*dlamk1
     params_test[k+3] = params[k+3] - alpha*dmuk1
     f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while(f_before <= f_after and i < 25):
          alpha = alpha*0.5
          params_test[k] = params[k] - alpha*dlamk
          params_test[k+1] = params[k+1] - alpha*dmuk
          params_test[k+2] = params[k+2] - alpha*dlamk1
          params_test[k+3] = params[k+3] - alpha*dmuk1
          f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if(f_before <= f_after):
          alpha = 0
     return alpha

def backTr_4Blockn(alpha0, n, dlamn, dmun, dlamn1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking for the 4 pass update for the last block [lamn, mun, lamn1]
     '''
     i = 0
     alpha = alpha0
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     params_test = np.copy(params)
     params_test[2*n - 3] = params[2*n - 3] - alpha*dlamn
     params_test[2*n - 2] = params[2*n - 2] - alpha*dmun
     params_test[2*n - 1] = params[2*n - 1] - alpha*dlamn1
     f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while(f_before <= f_after and i < 25):
          alpha = alpha*0.5
          params_test[2*n - 3] = params[2*n - 3] - alpha*dlamn
          params_test[2*n - 2] = params[2*n - 2] - alpha*dmun
          params_test[2*n - 1] = params[2*n - 1] - alpha*dlamn1
          f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if(f_before <= f_after):
          alpha = 0
     return alpha
     
          

def get_sk(muk, lamk):
     if(muk > lamk):
          sk = 1
     else:
          sk = -1
     return sk


def gradient_TY(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Gradient of the objective function TY
     '''
     n = len(listxk) - 2
     gradient = np.zeros((2*n + 1))
     muk = params[0]
     lamk1 = params[1]
     muk1 = params[2]
     B01_mu = itt.gradientBoundary(muk, x0, listB0k[0], listxk[1], listBk[0])
     yk1 = itt.hermite_boundary(lamk1, x0, listB0k[1], listxk[2], listBk[1])
     zk = itt.hermite_boundary(muk, x0, listB0k[0], listxk[1], listBk[0])
     gradient[0] = partial_fObj_mu1(params[0], x0, T0, grad0, x1, T1, grad1, B01_mu, yk1, zk)
     gradient[1] = partial_fObj_lambdak1(muk, muk1, lamk1, x0, listB0k[0], listxk[1], listBk[0], listB0k[1], listxk[2], listBk[1], listIndices[1], listIndices[0])
     for j in range(1,n):
          k = 2*j
          etakM1 = listIndices[j-1]
          etak = listIndices[j]
          etak1 = listIndices[j+1]
          lamk = params[k-1]
          muk = params[k]
          lamk1 = params[k+1]
          gradient[k] = partial_fObj_muk(muk, lamk, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], listB0k[j+1], listxk[j+2], listBk[j+1], etak, etakM1)
          gradient[k+1] = partial_fObj_lambdak1(muk, muk1, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], listB0k[j+1], listxk[j+2], listBk[j+1], etak1, etak)
     gradient[2*n] = 0
     return gradient
          

def forwardPassUpdate(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
    '''
    Updates blocks [muk, lamk1]
    Forward pass for the block coordinate gradient descent goes from k = 0 to k = 2n-1
    '''
    # First coordinate in the block: muk
    params = np.copy(params0)
    muk = params[0]
    lamk1 = params[1]
    B0k_muk = itt.gradientBoundary(muk, x0, listB0k[0], listxk[1], listBk[0])
    yk1 = itt.hermite_boundary(lamk1, x0, listB0k[1], listxk[2], listBk[1])
    zk = itt.hermite_boundary(muk, x0, listB0k[0], listxk[1], listBk[0])
    # Compute direction for muk
    partial_muk = partial_fObj_mu1(muk, x0, T0, grad0, x1, T1, grad1, B0k_muk, yk1, zk)
    #print("  partial muk: ", partial_muk)
    alpha = backTr_coord(1, 0, partial_muk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    muk = muk - alpha*partial_muk
    muk = project_box(muk)
    params[0] = muk
    # Now coordinate descent for lamk1
    muk1 = params[2]
    partial_lamk1 = partial_fObj_lambdak1(muk, muk1, lamk1, x0, listB0k[0], listxk[1], listBk[0], listB0k[1], listxk[2], listBk[1], listIndices[1], listIndices[0])
    #print("  partial lamk1: ", partial_lamk1)
    alpha = backTr_coord(1, 1, partial_lamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    lamk1 = lamk1 - alpha*partial_lamk1
    #print("   initial before projecting: ", muk, " ,  ", lamk1)
    # Project back so that it is feasible
    muk, lamk1 = project_block(muk, lamk1, x0, listB0k[0], listxk[1], listBk[0], listB0k[1], listxk[2], listBk[1])
    params[0] = muk
    params[1] = lamk1
    #print("  initial", params)

    n = len(listxk) - 2
    for j in range(1, n):
        k = 2*j
        params_test = np.copy(params)
        etakM1 = listIndices[j-1]
        etak = listIndices[j]
        etak1 = listIndices[j+1]
        lamk = params[k-1]
        muk = params[k]
        lamk1 = params[k+1]
        f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        # Compute direction for muk
        partial_muk = partial_fObj_muk(muk, lamk, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], listB0k[j+1], listxk[j+2], listBk[j+1], etak, etakM1)
        #print("  partial muk: ", partial_muk)
        alpha = backTr_coord(1, k, partial_muk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        muk = muk - alpha*partial_muk
        muk = project_box(muk)
        params[k] = muk
        # Now coordinate descent for lamk1
        muk1 = params[k+1]
        partial_lamk1 = partial_fObj_lambdak1(muk, muk1, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], listB0k[j+1], listxk[j+2], listBk[j+1], etak1, etak)
        #print("  partial lamk1: ", partial_lamk1)
        alpha = backTr_coord(1, k+1, partial_lamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        lamk1 = lamk1 - alpha*partial_lamk1
        # Project back so that it is feasible
        #print("   before projecting: ", muk, " ,  ", lamk1)
        muk, lamk1 = project_block(muk, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], listB0k[j+1], listxk[j+2], listBk[j+1])
        params_test[k] = muk
        params_test[k+1] = lamk1
        f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        if(f_test < f_before):
             params[k] = muk
             params[k+1] = lamk1
        #print("   coordinate descent: ", params)
        # print("  At block iteration ", j, "  muk = ", params[k], "   lamk1 = ", params[k+1])
        # print("  At block iteration ", j, " params:", params)
        #print("  Objective function value: ", fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk) )
        #itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params)
    return params

def forwardMiddlePassUpdate(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Updates blocks [lamk, muk]
     Forward middle pass for the block coordinate gradient descent goes from k = 1 to k = 2n-1
     Doesnt update mu1 or lam(n+1)
     '''
     params = np.copy(params0)
     n = len(listxk) - 2
     for j in range(1,(n-1)):
          params_test = np.copy(params)
          f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          #print("fObj before middle update: ", f_before)
          k = 2*j - 1
          etakM1 = listIndices[j-1]
          etak = listIndices[j]
          etak1 = listIndices[j+1]
          lamk1 = params[k+2]
          lamk = params[k]
          muk = params[k+1]
          lamkM1 = params[k-2]
          mukM1 = params[k-1]
          # Compute directions for lamk1, muk1
          partial_lamk = partial_fObj_lambdak1(mukM1, muk, lamk, x0, listB0k[j-1], listxk[j], listBk[j-1], listB0k[j], listxk[j+1], listBk[j], etak1, etak)
          print("  partial lamk: ", partial_lamk)
          partial_muk = partial_fObj_muk(muk, lamk, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], listB0k[j+1], listxk[j+2], listBk[j+1], etak, etakM1)
          print("  partial muk: ", partial_muk)
          # Compute step size
          alpha = backTr_middleBock(1, k, partial_lamk, partial_muk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          # Compute update
          lamk = lamk - alpha*partial_lamk
          muk = muk - alpha*partial_muk
          #print(" Proposed lamk ", lamk, "  proposed muk ", muk)
          # Project back all the necessary parts involved in this update
          mukM1, lamk, muk, lamk1 = project_middleBlock(k, mukM1, lamk, muk, lamk1, x0, listB0k, listxk, listBk)
          params_test[k-1] = mukM1
          params_test[k] = lamk
          params_test[k+1] = muk
          params_test[k+2] = lamk1
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          #print("fObj after middle update: ", f_test)
          if (f_test < f_before):
               params[k-1] = mukM1
               params[k] = lamk
               params[k+1] = muk
               params[k+2] = lamk1
     return params

def forward4PassUpdate(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Updates blocks [lamk, muk, lamk1, muk1]
     for the first one it just updates [mu1, lam2, mu2]
     for the final it just updates [mun, lamn, lam(n+1)]
     '''
     params = np.copy(params0)
     n = len(listxk) - 2
     # FOR THE FIRST BLOCK [mu1, lam2, mu2]
     mu1 = params[0]
     lam2 = params[1]
     mu2 = params[2]
     B01_mu1 = itt.gradientBoundary(mu1, x0, listB0k[0], listxk[1], listBk[0])
     z1 = itt.hermite_boundary(mu1, x0, listB0k[0], listxk[1], listBk[0])
     y2 = itt.hermite_boundary(lam2, x0, listB0k[1], listxk[2], listBk[1])
     z2 = itt.hermite_boundary(mu2, x0, listB0k[1], listxk[2], listBk[1])
     # Compute directions
     partial_mu1 = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B01_mu1, y2, z1)
     partial_lam2 = partial_fObj_lambdak1(mu1, mu2, lam2, x0, listB0k[0], listxk[1], listBk[0], listB0k[1], listxk[2], listBk[1], listIndices[1], listIndices[0])
     if(n > 1):
          partial_mu2 = partial_fObj_muk(mu2, lam2, params[3], x0, listB0k[2], listxk[3], listBk[2], listB0k[3], listxk[4], listBk[3], listIndices[1], listIndices[0])
     else:
          partial_mu2 = 0 # this is when we only have one section mu2 = 1 always
     # Backtracking
     alpha = backTr_4Block0(1, partial_mu1, partial_lam2, partial_mu2, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     # Update
     mu1 = mu1 - alpha*partial_mu1
     lam2 = lam2 - alpha*partial_lam2
     mu2 = mu2 - alpha*partial_mu2
     # Project back if necessary
     if(n>1):
          lam3 = params[3]
          mu1, lam2, mu2 = project_4Block0(mu1, lam2, mu2, lam3, x0, listB0k[0], listxk[1], listBk[0], listB0k[1], listxk[2], listBk[1], listB0k[2], listxk[3], listBk[2])
     else:
          mu1, lam2 = project_block(mu1, lam1, x0, listB0k[0], listxk[1], listBk[0], listB0k[1], listxk[2], listBk[1])
     params[0] = mu1
     params[1] = lam2
     params[2] = mu2

     # NOW WE CYCLE
     for j in range(1, n-1):
          k = 2*j - 1 # From 1 to 2n - 5
          mukM1 = params[k-1]
          lamk = params[k]
          muk = params[k+1]
          lamk1 = params[k+2]
          muk1 = params[k+3]
          lamk2 = params[k+4]
          etakM1 = listIndices[j-1]
          etak = listIndices[j]
          etak1 = listIndices[j+1]
          B0kM1 = listB0k[j-1]
          xkM1 = listxk[j]
          BkM1 = listBk[j-1]
          B0k = listB0k[j]
          xk = listxk[j+1]
          Bk = listBk[j]
          B0k1 = listB0k[j+1]
          xk1 = listxk[j+2]
          Bk1 = listBk[j+1]
          B0k2 = listB0k[j+2]
          xk2 = listxk[j+3]
          Bk2 = listBk[j+2]
          # Compute directions
          partial_lamk = partial_fObj_lambdak1(mukM1, muk, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, etakM1, etak)
          partial_muk = partial_fObj_muk(muk, lamk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etakM1)
          partial_lamk1 = partial_fObj_lambdak1(muk, muk1, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etak1)
          partial_muk1 = partial_fObj_muk(muk1, lamk1, lamk2, x0, B0k1, xk1, Bk1, B0k2, xk2, Bk2, etak1, etak)
          alpha = backTr_4Block(1, k, partial_lamk, partial_muk, partial_lamk1, partial_muk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          lamk = lamk - alpha*partial_lamk
          muk = muk - alpha*partial_muk
          lamk1 = lamk1 - alpha*partial_lamk1
          muk1 = muk1 - alpha*partial_muk1
          # Project back
          lamk, muk, lamk1, muk1 = project_4Block(mukM1, lamk, muk, lamk1, muk1, lamk2, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, B0k1, xk1, Bk1, B0k2, xk2, Bk2)
     if(n>1):
     # FOR THE LAST BLOCK [lamn, mun, lamn1]
          munM1 = params[2*n -2]
          lamn = params[2*n-3]
          mun = params[2*n-2]
          lamn1 = params[2*n-1]
          B0nM1 = listB0k[n-2]
          xnM1 = listxk[n-1]
          BnM1 = listBk[n-2]
          B0n = listB0k[n-1]
          xn = listxk[n]
          Bn = listBk[n-1]
          B0n1 = listB0k[n]
          xn1 = listxk[n+1]
          Bn1 = listBk[n]
          etanM1 = listIndices[n-2]
          etan = listIndices[n-1]
          etan1 = listIndices[n]
          partial_lamn = partial_fObj_lambdak1(munM1, mun, lamn, x0, B0nM1, xnM1, BnM1, B0n, xn, Bn, etanM1, etan)
          partial_mun = partial_fObj_muk(mun, lamn, lamn1, x0, B0n, xn, Bn, B0n1, xn1, Bn1, etan, etanM1)
          partial_lamn1 = partial_fObj_lambdak1(mun, 1, lamn1, x0, B0n, xn, Bn, B0n1, xn1, Bn1, etan, etan1)
          alpha = backTr_4Blockn(1, n, partial_lamn, partial_mun, partial_lamn1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          lamn = lamn - alpha*partial_lamn
          mun = mun - alpha*partial_mun
          lamn1 = lamn1 - alpha*partial_lamn1
          lamn, mun, lamn1 = project_4Blockn(munM1, lamn, mun, lamn1, x0, B0nM1, xnM1, BnM1, B0n, xn, Bn, B0n1, xn1, Bn1)
          params[2*n-3] = lamn
          params[2*n-2] = mun
          params[2*n-1] = lamn1
     return params
          
     

def backwardPassUpdate(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
    '''
    Backward pass for the block coordinate gradient descent goes from k = 2n-1 to k = 0
    '''
    params = np.copy(params0)
    n = len(listxk) -2
    for j in range(n, 1, -1):
        # Going backwards through the blocks
        k = 2*j - 1 # from 2n-1 to 3
        params_test = np.copy(params)
        f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        # First coordinate in the block: lamk1
        muk1 = params[k+1]
        lamk1 = params[k]
        muk = params[k-1]
        lamk = params[k-2]
        etak1 = listIndices[j]
        etak = listIndices[j-1]
        etakM1 = listIndices[j-2]
        B0k1 = listB0k[j]
        xk1 = listxk[j+1]
        Bk1 = listBk[j]
        B0k = listB0k[j-1]
        xk = listxk[j]
        Bk = listBk[j-1]
        # Compute directions
        partial_lamk1 = partial_fObj_lambdak1(muk, muk1, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etak1)
        # Compute step size
        alpha_lamk1 = backTr_coord(1, k, partial_lamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        lamk1 = lamk1 - alpha_lamk1*partial_lamk1
        lamk1 = project_box(lamk1)
        params[k] = lamk1
        # Second coordinate in the block: muk
        # Compute direction for muk
        partial_muk = partial_fObj_muk(muk, lamk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etakM1)
        alpha_muk = backTr_coord(1, k-1, partial_muk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        muk = muk - alpha_muk*partial_muk
        # Project back so that it is feasible
        #print("   before projecting: ", muk, " ,  ", lamk1)
        muk, lamk1 = project_blockBackwards(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
        params_test[k] = muk
        params_test[k+1] = lamk1
        f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
        if(f_test < f_before):
             params[k] = muk
             params[k+1] = lamk1
        #print("   coordinate descent: ", params)
        #itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params)
    # Now we update mu1 and lam2
    # First we update lam2
    f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    params_test = np.copy(params)
    lamk1 = params[1]
    muk1 = params[2]
    muk = params[0]
    etak1 = listIndices[1]
    etak = listIndices[0]
    B0k1 = listB0k[1]
    xk1 = listB0k[2]
    Bk = listBk[1]
    B0k = listB0k[0]
    xk = listxk[1]
    Bk = listBk[0]
    partial_lamk1 = partial_fObj_lambdak1(muk, muk1, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etak1)
    # Compute direction
    alpha_lamk1 = backTr_coord(1, 1, partial_lamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    # Compute step size
    alpha_lamk1 = backTr_coord(1, 1, partial_lamk1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    lamk1 = lamk1 - alpha_lamk1*partial_lamk1
    lamk1 = project_box(lamk1)
    # For mu1
    y2 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
    z1 = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
    B01_mu = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
    partial_muk = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B01_mu, y2, z1)
    alpha_muk = backTr_coord(1, 0, partial_muk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    muk = muk - alpha_muk*partial_muk
    # Project back so that it is feasible
    muk, lamk1 = project_blockBackwards(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
    params_test[0] = muk
    params_test[1] = lamk1
    f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    if(f_test < f_before):
         params[0] = muk
         params[1] = lamk1
    #print("  initial", params)
    ##itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params)
    return params



def blockCoordinateGradient(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol):
    '''
    Block coordinate projected gradient descent for an update in a triangle fan without tops
    '''
    params = np.copy(params0)
    params = np.append(params, [1])
    print(params)
    listObjVals = []
    listGrads = []
    # Compute the initial gradient
    gradk = gradient_TY(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
    iter = 0
    while(norm(gradk)>tol and iter<maxIter):
        # Start with the forward pass update
         print("\nStarting iteration: ", iter)
         paramsUp = forwardPassUpdate(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         fVal = fObj_noTops(paramsUp, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listObjVals.append(fVal)
         gradk = gradient_TY(paramsUp, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listGrads.append(norm(gradk))
         print("Objective function value: ", fVal )
         #itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *paramsUp[:-1])
         paramsMiddle = forwardMiddlePassUpdate(paramsUp, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         fVal = fObj_noTops(paramsMiddle, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listObjVals.append(fVal)
         gradk = gradient_TY(paramsMiddle, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listGrads.append(norm(gradk))
         print("Objective function value: ", fVal )
         #itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *paramsMiddle[:-1])
         paramsDown = backwardPassUpdate(paramsMiddle, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         fVal = fObj_noTops(paramsDown, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listObjVals.append(fVal)
         gradk = gradient_TY(paramsDown, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listGrads.append(norm(gradk))
         print("Objective function value: ", fVal )
         itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *paramsDown[:-1])
         params4Block = forward4PassUpdate(paramsDown, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         fVal = fObj_noTops(params4Block, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listObjVals.append(fVal)
         gradk = gradient_TY(params4Block, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         listGrads.append(norm(gradk))
         print("Objective function value: ", fVal )
         params = params4Block
         gradk = gradient_TY(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
         iter += 1
         itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params[:-1])
    return params, listObjVals, listGrads
        
    
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

#itt.plotFan3(x0, B01, B02, B03, B0Hat, x1, B1, x2, B2, x3, B3, xHat, BHat, mu1 = mu1, lam2 = lam2, mu2 = mu2, lam3 = lam3, mu3 = mu3, lam4 = lam4)

#plt.scatter(xSource[0], xSource[1], s = 20, c = "#ff00f2")
#plt.plot([xSource[0], xHat[0]], [xSource[1], xHat[1]], linestyle = ':', c = "#ab00a3")


#f_init = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)


##### Test the projection

# lam2 = 0
# mu1 = 1
# Bk_muk = itt.gradientBoundary(mu1, x0, B01, x1, B1)
# Bk1_lamk1 = itt.gradientBoundary(lam2, x0, B02, x2, B2)
# yk1 = itt.hermite_boundary(lam2, x0, B02, x2, B2)
# zk = itt.hermite_boundary(mu1, x0, B01, x1, B1)

# #itt.plotFan3(x0, B01, B02, B03, B0Hat, x1, B1, x2, B2, x3, B3, xHat, BHat, mu1 = mu1, lam2 = lam2, mu2 = mu2, lam3 = lam3, mu3 = mu3, lam4 = lam4)

# mu1, lam2 = project_block(mu1, lam2, Bk_muk, Bk1_lamk1, yk1, zk, x0, x2, B02, B2)

# #itt.plotFan3(x0, B01, B02, B03, B0Hat, x1, B1, x2, B2, x3, B3, xHat, BHat, mu1 = mu1, lam2 = lam2, mu2 = mu2, lam3 = lam3, mu3 = mu3, lam4 = lam4)


##### Test the forward and backward updates

print("Start test foward pass update \n\n")

mu1 = 0.5
lam2 = 0.65
mu2 = 0.5
lam3 = 0.75
mu3 = 0.5
lam4 = 0.45
params = [mu1, lam2, mu2, lam3, mu3, lam4]



# Compute the projected gradient descent

maxIter = 20
tol = 1e-8

paramsOpt, listVals, listGrads = blockCoordinateGradient(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol)
#itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *params)
#itt.plotFan3(x0, *listB0k, listxk[1], listBk[0], listxk[2], listBk[1], listxk[3], listBk[2], xHat, listBk[3], *paramsOpt[:-1])

#plt.close("all")

gradient = gradient_TY(paramsOpt, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

params2 = paramsOpt - 0.001*gradient

gradient2 = gradient_TY(params2, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

print("\nGradient Final: \n", gradient)

fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
plt.loglog(range(0, len(listVals)), listVals)


fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
plt.loglog(range(0, len(listGrads)), listVals)

#paramsUp = forwardPassUpdate(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

# print(paramsUp)
# print("Objective function value:", fObj_noTops(paramsUp, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk) )

# print(" \n\n Start test backward pass update \n\n")


# paramsDown, grad2T1 = backwardPassUpdate(paramsUp, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

# print(paramsDown)

# print("Objective function value:", fObj_noTops(paramsDown, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk) )


# print("Start test2 foward pass update \n\n")

# paramsUp = forwardPassUpdate(paramsDown, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

# print(paramsUp)
# print("Objective function value:", fObj_noTops(paramsUp, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk) )

# print(" \n\n Start test2 backward pass update \n\n")


# paramsDown, grad2T2 = backwardPassUpdate(paramsUp, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)

# print(paramsDown)

# print("Objective function value:", fObj_noTops(paramsDown, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk) )




