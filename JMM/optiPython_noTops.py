
import pdb # FOR DEBUGGING

import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from scipy.optimize import root_scalar # to project back blocks [mu_k, lambda_k+1]
import colorcet as cc
import matplotlib.colors as clr
import json
from numba import njit


colormap2 = "cet_linear_worb_100_25_c53_r"
colormap2_r = "cet_linear_worb_100_25_c53"



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


def forwardPassUpdate_noTops(params0, gammas, theta_gamma, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     gammas : radius for the circle centered at [lamk, muk], if such circle intersects the line lamk = muk then do a close update, it is an array of size nRegions - 1 (notice that it can be an array of size 0)
     theta_gamma : rate of decrease in the circle centered at [lamk, muk]
     Updates blocks     [mu1]
                        [lam2, mu2]
                        [lamk, muk]
                        [lamn1]
     '''
     # First parameter to update: mu1
     params = np.copy(params0)
     mu1 = params[0]
     lam2 = params[1]
     B0k = listB0k[0]
     xk = listxk[1]
     Bk = listBk[0]
     B0k1 = listB0k[1]
     xk1 = listxk[2]
     Bk1 = listBk[1]
     B0k_muk = itt.gradientBoundary(mu1, x0, B0k, xk, Bk)
     yk1 = itt.hermite_boundary(lam2, x0, B0k1, xk1, Bk1)
     zk = itt.hermite_boundary(mu1, x0, B0k, xk, Bk)
     # Compute direction for muk
     dmuk = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B0k_muk, yk1, zk)
     alpha = backTr_coord_noTops(2, 0, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     mu1 = mu1 - alpha*dmuk
     mu1 = project_mukGivenlamk1(mu1, lam2, x0, B0k, xk, Bk, B0k1, xk1, Bk1) # project back so that it is feasible
     params[0] = mu1
     # Now we start with the blocks of size 2
     nRegions = len(listxk) - 2
     for j in range(1, nRegions):
          k = 2*j # From 2 to 2nRegions - 2
          gamma = gammas[j - 1]
          mukM1 = params[k-2]
          lamk = params[k-1]
          muk = params[k]
          lamk1 = params[k+1]
          B0kM1 = listB0k[j-1]
          xkM1 = listxk[j]
          BkM1 = listBk[j-1]
          B0k = listB0k[j]
          xk = listxk[j+1]
          Bk = listBk[j]
          B0k1 = listB0k[j+1]
          xk1 = listxk[j+2]
          Bk1 = listBk[j+1]
          etakM1 = listIndices[j-1]
          etak = listIndices[j]
          etak1 = listIndices[j+1]
          # Compute directions
          dlamk = partial_fObj_recCr(mukM1, muk, lamk, x0, B0kM1, xkM1, BkM1, x0, B0k, xk, Bk, etakM1, etak)
          dmuk = partial_fObj_shCr(muk, lamk, lamk1, x0, B0k, xk, Bk, x0, B0k1, xk1, Bk1, etak, etakM1)
          # See if we need to do a close update or not
          r = close_to_identity(lamk, muk)
          if( r <= gamma ):
               # Meaning we have to do a close update
               lamk, muk = backTrClose_block_noTops(2, k-1, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
               gamma = gamma*theta_gamma
               gammas[j - 1] = gamma # Update this gamma
          else:
               # We don't have to consider a close update
               lamk, muk = backTr_block_noTops(2, k-1, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          # Either way we need to project back to the feasible set
          lamk = project_lamk1Givenmuk(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
          muk = project_mukGivenlamk1(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
          # Update
          params[k-1] = lamk
          params[k] = muk
     # Finally update lamn1
     mun = params[2*nRegions - 2]
     lamn1 = params[2*nRegions - 1]
     mun1 = 1 # Always
     B0kM1 = listB0k[nRegions - 1]
     xkM1 = listxk[nRegions ]
     BkM1 = listBk[nRegions - 1]
     B0k = listB0k[nRegions ]
     xk = listxk[nRegions + 1]
     Bk = listBk[nRegions ]
     etakM1 = listIndices[nRegions - 2]
     etak = listIndices[nRegions - 1]
     # Compute direction
     dlamn1 = partial_fObj_recCr(mun, mun1, lamn1, x0, B0kM1, xkM1, BkM1, x0, B0k, xk, Bk, etakM1, etak)
     # Compute step size
     alpha = backTr_coord_noTops(2, (2*nRegions - 1), dlamn1, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     # Update
     lamn1 = lamn1 - alpha*dlamn1
     # Project back
     lamn1 = project_lamk1Givenmuk(mun, lamn1, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
     # Update
     params[2*nRegions - 1] = lamn1
     return params, gammas
     


@njit
def project_box(param):
     if(param < 0):
          param = 0
     elif(param > 1):
          param = 1
     else:
          param = param
     return param

def backTr_coord_noTops(alpha0, k, d, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Coordinate backtracking to find the step size for coordinate gradient descent
     '''
     i = 0
     alpha = alpha0*0.2/(max(abs(d), 1))
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     params_test = np.copy(params)
     params_test[k] = params[k] - alpha*d
     f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_after < f_before and i < 8 ):
          alpha = alpha*1.3 # increase step size
          params_test[k] = params[k] - alpha*d
          f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     i = 0
     while( f_before <= f_after and i < 25):
          alpha = alpha*0.2
          params_test[k] = params[k] - alpha*d
          f_after = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_after):
          alpha = 0
     return alpha

@njit
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
     gradient[1] = partial_fObj_recCr1(muk, muk1, lamk1, x0, listB0k[0], listxk[1], listBk[0], listB0k[1], listxk[2], listBk[1], listIndices[1], listIndices[0])
     for j in range(1,n):
          k = 2*j
          etakM1 = listIndices[j-1]
          etak = listIndices[j]
          etak1 = listIndices[j+1]
          lamk = params[k-1]
          muk = params[k]
          lamk1 = params[k+1]
          muk1 = params[k+2]
          gradient[k] = partial_fObj_shCr(muk, lamk, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], x0, listB0k[j+1], listxk[j+2], listBk[j+1], etak, etakM1)
          gradient[k+1] = partial_fObj_recCr1(muk, muk1, lamk1, x0, listB0k[j], listxk[j+1], listBk[j], listB0k[j+1], listxk[j+2], listBk[j+1], etak1, etak)
     gradient[2*n] = 0
     return gradient



def blockCoordinateGradient(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, maxIter, tol, theta_gamma = 1, plotSteps = True, saveIterates = False):
     '''
     Block coordinate subgradient descent (modified) for an update in a triangle fan without
     tops. Inspired by gradient sampling but in this case we know where our (geometric)
     singularities are and thus we don't need to calculate explicitly the convex hull
     '''
     params = np.copy(params0)
     params = np.append(params, [1])
     #print(" Initial params: \n", params)
     # Initialize the useful things
     listObjVals = []
     listGrads = []
     listChangefObj = []
     listIterates = []
     listGradIterations = []
     nRegions = len(listxk) -2
     gammas = 0.05*np.ones((nRegions - 1)) # Might be an array of length 0, it's fine
     gradk = gradient_TY(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     norm_gradk = norm(gradk)
     fVal = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     listGrads.append(norm_gradk)
     listObjVals.append(fVal)
     iter = 0
     change_fVal = 1
     while( abs(change_fVal) > tol and iter < maxIter):
          # Start with a forward pass
          params, gammas = forwardPassUpdate_noTops(params, gammas, theta_gamma, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          gradk = gradient_TY(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          norm_gradk = norm(gradk)
          fVal_prev = fVal
          fVal = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          change_fVal = fVal_prev - fVal
          listGrads.append(norm_gradk)
          listObjVals.append(fVal)
          listChangefObj.append( change_fVal )
          iter += 1
          if plotSteps:
               itt.plotFann(x0, listB0k, listxk, listBk, params = params)
          if saveIterates:
               listIterates.append( params )
               listGradIterations.append( gradk)
     if saveIterates:
          return params, listObjVals, listGrads, listChangefObj, listIterates, listGradIterations
     else:
          return params, listObjVals, listGrads, listChangefObj
     



