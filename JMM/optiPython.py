


# Before going into C with the optimization rutine
# we test the projected coordinate gradient descent here
# to make sure it works or it makes sense to try this approach




import pdb


import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt
from scipy.optimize import root_scalar # to project back blocks [mu_k, lambda_k+1]
import colorcet as cc
import matplotlib.colors as clr


# colormap2  = clr.LinearSegmentedColormap.from_list('Retro',
#                                                    [(0,    '#000000'),
#                                                     (0.1, '#2c3454'),
#                                                     (0.25, '#0033ff'),
#                                                     (0.60, '#00f3ff'),
#                                                     (1,    '#e800ff')], N=256)

colormap2 = "cet_linear_worb_100_25_c53_r"
colormap2_r = "cet_linear_worb_100_25_c53"

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

def findRtan(r, xkM1, xk, BkM1Bk_0, BkM1Bk_1, pointFrom):
     '''
     Find a_tan
     '''
     a_tan = itt.hermite_boundary(r, xkM1, BkM1Bk_0, xk, BkM1Bk_1)
     BkM1Bk_tan = itt.gradientBoundary(r, xkM1, BkM1Bk_0, xk, BkM1Bk_1)
     N_tan = np.array([-BkM1Bk_tan[1], BkM1Bk_tan[0]])
     return np.dot(N_tan, a_tan - pointFrom)


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
     if( norm(y2 - z1) < 1e-8 ):
          return der_hermite_inter
     else:
          return der_hermite_inter - np.dot(B01_mu, y2 - z1)/norm(y2 - z1)


def partial_fObj_recCr(shFrom, shTo, rec, x0From, B0From, x1From, B1From, x0To, B0To, x1To, B1To, etaInside, etaOutside):
     '''
     Partial of the objective function with respect to a receiver that creeps to a shooter
     (originally the receiver is lamk, the shooter where it comes from is mukM1 and the shooter
     to where it creeps is muk)
     '''
     shooterFrom = itt.hermite_boundary(shFrom, x0From, B0From, x1From, B1From)
     receiver = itt.hermite_boundary(rec, x0To, B0To, x1To, B1To)
     B_atShooterTo = itt.gradientBoundary(shTo, x0To, B0To, x1To, B1To)
     B_atReceiver = itt.gradientBoundary(rec, x0To, B0To, x1To, B1To)
     secondDer_B_atReceiver = itt.secondDer_Boundary(rec, x0To, B0To, x1To, B1To)
     B_halves = itt.gradientBoundary( (shTo + rec)/2, x0To, B0To, x1To, B1To)
     secondDer_Bhalves_atReceiver = itt.secondDer_Boundary( (shTo + rec)/2, x0To, B0To, x1To, B1To)
     perL_receiver = partial_L_lamk(shTo, rec, B_atShooterTo, B_halves, secondDer_Bhalves_atReceiver, B_atReceiver, secondDer_B_atReceiver)
     etaMin = min(etaInside, etaOutside)
     if( shFrom == rec ):
          return etaMin*perL_receiver
     else:
          return etaInside*np.dot(B_atReceiver, receiver - shooterFrom) + etaMin*perL_receiver
     


def partial_fObj_shCr(sh, recFrom, recTo, x0From, B0From, x1From, B1From, x0To, B0To, x1To, B1To, etaInside, etaOutside):
     '''
     Partial of the objective function with respect to a shooter that comes from a creeping ray from a receiver
     '''
     shooter = itt.hermite_boundary(sh, x0From, B0From, x1From, B1From)
     receiverTo = itt.hermite_boundary(recTo, x0To, B0To, x1To, B1To)
     B_atShooter = itt.gradientBoundary(sh, x0From, B0From, x1From, B1From)
     secondDer_B_atShooter = itt.secondDer_Boundary(sh, x0From, B0From, x1From, B1From)
     B_halves = itt.gradientBoundary( (sh + recFrom)/2, x0From, B0From, x1From, B1From)
     secondDer_Bhalves_atShooter = itt.secondDer_Boundary( (sh + recFrom)/2, x0From, B0From, x1From, B1From)
     B_atReceiver = itt.gradientBoundary(recFrom, x0From, B0From, x1From, B1From)
     parL_shooter = partial_L_muk(sh, recFrom, B_atShooter, secondDer_B_atShooter, B_halves, secondDer_Bhalves_atShooter, B_atReceiver)
     etaMin = min(etaInside, etaOutside)
     if( sh == recTo ):
          return etaMin*parL_shooter
     else:
          return etaInside*np.dot(-B_atShooter, receiverTo - shooter)/norm(receiverTo - shooter) + etaMin*parL_shooter



def partial_fObj_recCr1(muk, muk1, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, etak, etak1):
     '''
     Partial of the objective function of the "next" receiver that creeps to a shooter
     '''
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_muk1 = itt.gradientBoundary(muk1, x0, B0k1, xk1, Bk1)
     secondDer_B0k1_lamk1 = itt.secondDer_Boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_halves = itt.gradientBoundary( (muk1 + lamk1)/2, x0, B0k1, xk1, Bk1)
     secondDer_B0k1halves_lamk1 = itt.secondDer_Boundary((muk1 + lamk1)/2, x0, B0k1, xk1, Bk1)
     B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
     perL_lamk1 = partial_L_lamk(muk1, lamk1, B0k1_muk1, B0k1_halves, secondDer_B0k1halves_lamk1, B0k1_lamk1, secondDer_B0k1_lamk1)
     etaMin = min(etak, etak1)
     if( muk == lamk1 ):
          return etaMin*perL_lamk1
     else:
          return etak*np.dot(B0k1_lamk1, yk1 - zk)/norm(yk1 - zk) + etaMin*perL_lamk1


def partial_fObj_recSt(shFrom, shTo, rec, x0From, B0From, x1From, B1From, x0To, B0To, x1To, B1To, etaInside, etaOutside):
     '''
     Partial of the objective function (generalized) with respect to a RECEIVER that shoots directly (with
     a straight line) to a shooter (thinking about ak and bk)
     '''
     shooterFrom = itt.hermite_boundary(shFrom, x0From, B0From, x1From, B1From)
     receiver = itt.hermite_boundary(rec, x0To, B0To, x1To, B1To)
     shooterTo = itt.hermite_boundary(shTo, x0To, B0To, x1To, B1To)
     B_atReceiver = itt.gradientBoundary(rec, x0To, B0To, x1To, B1To)
     if( rec == 0 and shTo == 0 and shFrom == 0):
          return 0
     elif( np.all(shooterTo == receiver) and np.any(receiver != shooterFrom) ):
          return etaOutside*np.dot(B_atReceiver, receiver - shooterFrom)/norm( receiver - shooterFrom)
     elif( np.all(receiver == shooterFrom)  ):
          return etaInside*np.dot(- B_atReceiver, shooterTo - receiver)/norm(shooterTo - receiver)
     else:
          return etaOutside*np.dot(B_atReceiver, receiver - shooterFrom)/norm( receiver - shooterFrom) + etaInside*np.dot(- B_atReceiver, shooterTo - receiver)/norm(shooterTo - receiver)


def partial_fObj_shSt(sh, recFrom, recTo, x0From, B0From, x1From, B1From, x0To, B0To, x1To, B1To, etaInside, etaOutside):
     '''
     Partial of the objective function (generalized) with respect to a SHOOTER that comes from a straight ray
     from a receiver and shoots in a straight line to another receiver
     '''
     receiverFrom = itt.hermite_boundary(recFrom, x0From, B0From, x1From, B1From)
     shooter = itt.hermite_boundary(sh, x0From, B0From, x1From, B1From)
     receiverTo = itt.hermite_boundary(recTo, x0To, B0To, x1To, B1To)
     B_atShooter = itt.gradientBoundary(sh, x0From, B0From, x1From, B1From)
     if( sh == 0 and recFrom == 0 and recTo == 0 ):
          return 0
     elif( np.all(receiverTo == shooter) and np.any(shooter != receiverFrom) ):
          return etaOutside*np.dot( B_atShooter, shooter - receiverFrom)/norm(shooter - receiverFrom)
     elif( np.all(shooter == receiverFrom) ):
          return etaInside*np.dot( -B_atShooter, receiverTo - shooter)/norm(receiverTo - shooter)
     else:
          return etaOutside*np.dot( B_atShooter, shooter - receiverFrom)/norm(shooter - receiverFrom) + etaInside*np.dot( -B_atShooter, receiverTo - shooter)/norm(receiverTo - shooter)


def partial_fObj_collapsedShooter(shFrom, sh, recTo, x0From, B0From, x1From, B1From, x0This, B0This, x1This, B1This, x0To, B0To, x1To, B1To, etaPrev, etaNext):
     '''
     Partial of the objective function (generalized) with respect to a "COLLAPSED SHOOTER" (one where both
     the shooter and a receiver on an edge have collapsed to the same point). This points comes from a
     straight ray from a shooter on another edge and shoots in a straight line to another receiver
     in another edge (we are dealing with 3 edges here)
     '''
     shooterFrom = itt.hermite_boundary(shFrom, x0From, B0From, x1From, B1From)
     collapsedShooter = itt.hermite_boundary(sh, x0This, B0This, x1This, B1This)
     B_collapsedShooter = itt.gradientBoundary(sh, x0This, B0This, x1This, B1This)
     receiverTo = itt.hermite_boundary(recTo, x0To, B0To, x1To, B1To)
     if( np.all(shooterFrom == collapsedShooter) and  np.all( receiverTo == collapsedShooter) ):
          return 0
     elif( np.all(shooterFrom == collapsedShooter) and np.any(receiverTo != collapsedShooter) ):
          return etaNext*np.dot( -B_collapsedShooter, receiverTo - collapsedShooter)/norm( receiverTo - collapsedShooter)
     elif( np.all( receiverTo == collapsedShooter) and np.any( shooterFrom != collapsedShooter) ):
          return etaPrev*np.dot( B_collapsedShooter, collapsedShooter - shooterFrom)/norm(collapsedShooter - shooterFrom)
     else:
          return etaPrev*np.dot( B_collapsedShooter, collapsedShooter - shooterFrom)/norm(collapsedShooter - shooterFrom) + etaNext*np.dot( -B_collapsedShooter, receiverTo - collapsedShooter)/norm( receiverTo - collapsedShooter)


def project_lamk1Givenmuk(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1):
    '''
    Project back lamk1 given muk
    '''
    lamk1 = project_box(lamk1)
    zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
    yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
    B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
    B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
    
    Nk_muk = np.array([-B0k_muk[1], B0k_muk[0]])
    Nk1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
    dotTestMin = np.dot( yk1 - zk, Nk_muk )
    dotTestMax = np.dot( yk1 - zk, Nk1_lamk1 )

    # See if this orientation is correct
    if( np.dot( Nk_muk, xk1 - x0) < 0):
         Nk_muk = -Nk_muk
    if( np.dot( Nk1_lamk1, x0 - xk) < 0):
         Nk1_lamk1 = -Nk1_lamk1
    
    #print("       dotTestMin: ", dotTestMin, "  dotTestMax: ", dotTestMax)
    # Test if lamk < lamMin
    if( dotTestMin < 0):
        # print("  failed dotTestMin project lamk1 given muk")
        # print("  zk: ", zk, "  yk1: ", yk1)
        # print("  muk: ", muk, " lamk1: ", lamk1)
        # print("  B0k: ", B0k, "  xk: ", xk, "  Bk: ", Bk)
        # print("  B0k1: ", B0k1, "  xk1: ", xk1, "  Bk1: ", Bk1)
        # Means that we need to find lambdaMin
        tMin = lambda lam: t1(lam, x0, xk1, B0k1, Bk1, zk, B0k_muk)
        # Find the root of tMin
        rootMin = root_scalar(tMin, method = "secant", x0 = 0.4, x1 = 0.5)
        lamk1 = rootMin.root
        yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
        B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
        Nk1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
        if( np.dot( Nk1_lamk1, x0 - xk) < 0):
             Nk1_lamk1 = -Nk1_lamk1
        dotTestMax = np.dot( yk1 - zk, Nk1_lamk1 )
        #print("       lambda < lambdaMin")
    if( dotTestMax < 0):
        # print("  failed dotTestMax project lamk1 given muk")
        # print("  zk: ", zk, "  yk1: ", yk1)
        # print("  muk: ", muk, " lamk1: ", lamk1)
        # print("  B0k: ", B0k, "  xk: ", xk, "  Bk: ", Bk)
        # print("  B0k1: ", B0k1, "  xk1: ", xk1, "  Bk1: ", Bk1)
        # Means that we need to find lambdaMax
        tMax = lambda lam: t2(lam, x0, xk1, B0k1, Bk1, zk)
        rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
        lamk1 = rootMax.root
        #print("       lambda > lambdaMax")
    lamk1 = project_box(lamk1) # Such that 0<=lamk1 <=1
    return lamk1

def project_mukGivenlamk1(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1):
     '''
     Project back muk given lamk1
     '''
     muk = project_box(muk)
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
     
     # Compute the normals
     N0k1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
     N0k_muk = np.array([-B0k_muk[1], B0k_muk[0]])

     # See if this orientation is correct
     if( np.dot( N0k_muk, xk1 - x0) < 0):
          N0k_muk = -N0k_muk
     if( np.dot( N0k1_lamk1, x0 - xk) < 0):
          N0k1_lamk1 = -N0k1_lamk1
     
     # Compute the tests
     dotTestMin =  np.dot(N0k1_lamk1, yk1 - zk)
     dotTestMax = np.dot(N0k_muk, yk1 - zk)
     if(dotTestMin<0 ):
          # print("  failed dotTestMin project muk given lamk1")
          # print("  zk: ", zk, "  yk1: ", yk1)
          # print("  muk: ", muk, " lamk1: ", lamk1)
          # print("  B0k: ", B0k, "  xk: ", xk, "  Bk: ", Bk)
          # print("  B0k1: ", B0k1, "  xk1: ", xk1, "  Bk1: ", Bk1)
          tMin = lambda mu: t4(mu, x0, xk, B0k, Bk, yk1, B0k1_lamk1)
          rootMin = root_scalar(tMin, method = "secant", x0 = 0.4, x1 = 0.5)
          muk = rootMin.root
          zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
          B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
          N0k_muk = np.array([-B0k_muk[1], B0k_muk[0]])
          if( np.dot( N0k_muk, xk1 - x0) < 0):
               N0k_muk = -N0k_muk
          dotTestMax = np.dot(N0k_muk, yk1 - zk)
     if(dotTestMax<0 ):
          # print("  failed dotTestMax project muk given lamk1")
          # print("  zk: ", zk, "  yk1: ", yk1)
          # print("  muk: ", muk, " lamk1: ", lamk1)
          # print("  B0k: ", B0k, "  xk: ", xk, "  Bk: ", Bk)
          # print("  B0k1: ", B0k1, "  xk1: ", xk1, "  Bk1: ", Bk1)
          tMax = lambda mu: t3(mu, x0, xk, B0k, Bk, yk1)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          muk = rootMax.root
     muk = project_box(muk)
     return muk

     
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



def backTrClose_block_noTops(alpha0, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [lamk, muk]
     it considers two directions: steepest descent
                                  steepest descent projected onto the line lamk = muk
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     d_middle = np.array([dlamk, dmuk]) # "normal" steespest descent
     params_test = np.copy(params)
     alpha = alpha0*0.2*1/(max(norm(d_middle), 1))
     params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
     params_test_proj = np.copy(params)
     params_test_proj[k:(k+2)] = project_ontoLine(params_test[k:(k+2)])
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( (f_test < f_before or f_test_proj < f_before) and i < 8):
          alpha = alpha*1.3 # Increase step size
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          params_test_proj[k:(k+2)] = project_ontoLine(params_test[k:(k+2)])
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     i = 0
     while( f_before <= f_test and f_before <= f_test_proj and i < 25 ):
          alpha = alpha*0.2
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          params_test_proj[k:(k+2)] = project_ontoLine(params_test[k:(k+2)])
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          f_test_proj = fObj_noTops(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test and f_before <= f_test_proj ):
          return params[k], params[k+1]
     elif( f_test < f_test_proj):
          return params_test[k], params_test[k+1]
     elif( f_test_proj <= f_test ):
          return params_test_proj[k], params_test_proj[k+1]

def backTr_block_noTops(alpha0, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk):
     '''
     Backtracking to find the next block [lamk, muk]
     it considers one direction:  steepest descent
     '''
     f_before = fObj_noTops(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     i = 0
     d_middle = np.array([dlamk, dmuk]) # "normal" steespest descent
     params_test = np.copy(params)
     alpha = alpha0*0.2*1/(max(norm(d_middle), 1))
     params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
     # Compare the function value
     f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
     while( f_test < f_before and i < 8):
          alpha = alpha*1.3 # Increase step size
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     i = 0
     while( f_before <= f_test  and i < 25):
          alpha = alpha*0.2
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          f_test = fObj_noTops(params_test, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk)
          i += 1
     if( f_before <= f_test ):
          return params[k], params[k+1]
     else:
          return params_test[k], params_test[k+1]


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
     


def plotResults(x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listB0k, listxk,
                listBk, params0, paramsOpt, listObjVals, listGradNorms, listChangefObj,
                listChangeParams = None, trueSol = None, contours = True,
                listBkBk1 = None, indCrTop = None, paramsCrTop0 = None, 
                indStTop = None, paramsStTop0 = None, 
                paramsCrTop = None, paramsStTop = None):
     '''
     Plots results from blockCoordinateGradient
     '''
     # First plot the original parameters
     itt.plotFann(x0, listB0k, listxk, listBk, params = params0,
                  listBkBk1 = listBkBk1, indCrTop = indCrTop,
                  paramsCrTop = paramsCrTop0, indStTop = indStTop,
                  paramsStTop = paramsStTop0)
     plt.title("Initial parameters")
     # Then plot the parameters found by this method
     itt.plotFann(x0, listB0k, listxk, listBk, params = paramsOpt,
                  listBkBk1 = listBkBk1, indCrTop = indCrTop,
                  paramsCrTop = paramsCrTop, indStTop = indStTop,
                  paramsStTop = paramsStTop)
     plt.title("Optimal parameters found")
     # Now plot the function value at each iteration
     fig = plt.figure(figsize=(800/96, 800/96), dpi=96) 
     plt.semilogy( range(0, len(listChangefObj)), listChangefObj, c = "#394664", linewidth = 0.8)
     plt.xlabel("Iteration")
     plt.ylabel("Change of function value")
     plt.title("Change of function value at each iteration")
     # Now plot the function value at each iteration
     fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
     plt.semilogy( range(0, len(listObjVals)), listObjVals, c = "#396064", linewidth = 0.8)
     plt.xlabel("Iteration")
     plt.ylabel("Function value")
     plt.title("Function value at each iteration")
     # Now plot the norm of the gradient at each iteration
     fig = plt.figure(figsize=(800/96, 800/96), dpi=96) 
     plt.semilogy( range(0, len(listGradNorms)), listGradNorms, c = "#4d3964", linewidth = 0.8)
     plt.xlabel("Iteration")
     plt.ylabel("Norm of (sub)gradient")
     plt.title("Norm of (sub)gradient at each iteration")
     # Plot the change in the parameters
     fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
     plt.semilogy( range(0, len(listChangeParams)), listChangeParams, c = "#394664", linewidth = 0.8)
     plt.xlabel("Iteration")
     plt.ylabel("Norm of change in parameters")
     plt.title("Norm of change in parameters at each iteration")
     # For pairs of parameters MU AND LAMBDA
     nRegions = len(listxk) - 2
     fObjMesh = np.empty((200,200))
     for k in range(2*nRegions - 1):
          # We plot level sets by changins sequential parameters (2 at a time)
          param1, param2 = np.meshgrid( np.linspace(0,1,200), np.linspace(0,1,200) )
          # Compute the solution, compute this level set
          for i in range(200):
               for j in range(200):
                    p1 = param1[i,j]
                    p2 = param2[i,j]
                    paramsMesh = np.copy(paramsOpt)
                    paramsMesh[k] = p1
                    paramsMesh[k+1] = p2
                    fObjMesh[i,j] = fObj_generalized(paramsMesh, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                                     listxk, listB0k, listBk, listBkBk1,
                                                     indCrTop = indCrTop, paramsCrTop = paramsCrTop,
                                                     indStTop = indStTop, paramsStTop = paramsStTop)
          # Plot it
          fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
          im = plt.imshow(fObjMesh, cmap = colormap2, extent = [0,1,0,1], origin = "lower")
          plt.scatter(paramsOpt[k], paramsOpt[k+1], c = "white", marker = "*", label = "optimum found")
          if(contours):
               plt.contour(param1[0, :], param1[0, :], fObjMesh, colors = ["white"], extent = [0,1,0,1], origin = "lower", levels = 30, linewidths = 0.5)
          plt.title("Level set of objective function")
          plt.xlabel("Parameter " + str(k))
          plt.ylabel("Parameter " + str(k+1))
          plt.legend()
          plt.colorbar(im)
     # For pairs of parameters R AND S for a creeping ray update
     if( indCrTop is not None):
          fObjMesh = np.empty((200,200))
          for kCrTop in range(len(indCrTop)):
               k = 2*kCrTop
               # We plot level sets by changins sequential parameters (2 at a time)
               rr, ss = np.meshgrid( np.linspace(0,1,200), np.linspace(0,1,200) )
               # Compute the solution, compute this level set
               for i in range(200):
                    for j in range(200):
                         rk = rr[i,j]
                         sk = ss[i,j]
                         paramsMesh = np.copy(paramsCrTop)
                         paramsMesh[k] = rk
                         paramsMesh[k+1] = sk
                         fObjMesh[i,j] = fObj_generalized(paramsOpt, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                                          listxk, listB0k, listBk, listBkBk1,
                                                          indCrTop = indCrTop, paramsCrTop = paramsMesh,
                                                          indStTop = indStTop, paramsStTop = paramsStTop)
               # Plot it
               fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
               im = plt.imshow(fObjMesh, cmap = colormap2, extent = [0,1,0,1], origin = "lower")
               plt.scatter(paramsCrTop[k], paramsCrTop[k+1], c = "white", marker = "*", label = "optimum found")
               if(contours):
                    plt.contour(rr[0, :], rr[0, :], fObjMesh, colors = ["white"], extent = [0,1,0,1], origin = "lower", levels = 30, linewidths = 0.5)
               plt.title("Level set of objective function, parameters for creeping rays on side edges changed, r" + str(kCrTop+1) + " s" + str(kCrTop+1))
               plt.xlabel("r" + str(kCrTop+1))
               plt.ylabel("s" + str(kCrTop+1))
               plt.legend()
               plt.colorbar(im)
     # For pairs of parameters R AND S for a creeping ray update
     if( indStTop is not None):
          fObjMesh = np.empty((200,200))
          for kStTop in range(len(indStTop)):
               k = 2*kStTop
               # We plot level sets by changins sequential parameters (2 at a time)
               rr, ss = np.meshgrid( np.linspace(0,1,200), np.linspace(0,1,200) )
               # Compute the solution, compute this level set
               for i in range(200):
                    for j in range(200):
                         rk = rr[i,j]
                         sk = ss[i,j]
                         paramsMesh = np.copy(paramsStTop)
                         paramsMesh[k] = rk
                         paramsMesh[k+1] = sk
                         fObjMesh[i,j] = fObj_generalized(paramsOpt, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                                          listxk, listB0k, listBk, listBkBk1,
                                                          indCrTop = indCrTop, paramsCrTop = paramsCrTop,
                                                          indStTop = indStTop, paramsStTop = paramsMesh)
               # Plot it
               fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
               im = plt.imshow(fObjMesh, cmap = colormap2, extent = [0,1,0,1], origin = "lower")
               plt.scatter(paramsStTop[k], paramsStTop[k+1], c = "white", marker = "*", label = "optimum found")
               if(contours):
                    plt.contour(rr[0, :], ss[0, :], fObjMesh, colors = ["white"], extent = [0,1,0,1], origin = "lower", levels = 30, linewidths = 0.5)
               plt.title("Level set of objective function, parameters for straight rays on side edges changed, r" + str(kStTop) + " s" + str(kStTop))
               plt.xlabel("r" + str(k))
               plt.ylabel("s" + str(k+1))
               plt.legend()
               plt.colorbar(im)
     # If we know the true solution for this triangle fan, plot the decrease in the error
     if trueSol is not None:
          errRel = np.array(listObjVals)
          errRel = abs(errRel-trueSol)/trueSol
          fig = plt.figure(figsize=(800/96, 800/96), dpi=96)
          plt.semilogy( range(0, len(listObjVals)), errRel, c = "#00011f", linewidth = 0.8)
          plt.xlabel("Iteration")
          plt.ylabel("Relative error")
          plt.title("Relative error at each iteration")

######################################################
######################################################
######################################################
######################################################
######################################################
# Optimization for a generalized triangle fan i.e. the tops of the triangle fan are curved parametric curves

def fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices, listxk, listB0k, listBk, listBkBk1,
                     indCrTop = None, paramsCrTop = None, indStTop = None, paramsStTop = None):
     '''
     Generalized objective function for when the tops on the triangle fan are also parametric curves.
     '''
     # Set indStTop if indCrTop is given
     if(paramsStTop is None or indStTop is None):
          indStTop = [-1]
          paramsStTop = [0,0]
     # Set indCrTop if indStTop is given
     if(paramsCrTop is None or indCrTop is None):
          indCrTop = [-1]
          paramsCrTop = [0,0]
     currentCrTop = 0
     currentStTop = 0
     n = len(listxk) - 2
     muk = params[0]
     etak = listIndices[0]
     Bk = listBk[0]
     B0k = listB0k[0]
     zk = itt.hermite_boundary(muk, x0, B0k, x1, Bk)
     sum = hermite_interpolationT(muk, x0, T0, grad0, x1, T1, grad1)
     for j in range(1, n+1):
          # j goes from 1 to n
          # We have to circle around all the regions
          k = 2*j - 1  # Starts in k = 1, all the way to k = 2n - 1
          nTop = j # Number of boundary on the top that we are considering
          mukM1 = params[k-1]
          lamk = params[k]
          muk = params[k+1]
          B0k = listB0k[j]
          xkM1 = listxk[j]
          xk = listxk[j+1]
          Bk = listBk[j]
          BkBk1_0 = listBkBk1[k-1] # grad of hkhk1 at xk
          BkBk1_1 = listBkBk1[k] # grad of hkhk1 at xk1
          etakPrev = etak
          etak = listIndices[j]
          etaMin = min(etakPrev, etak)
          # Compute the points
          zkPrev = zk
          zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
          yk = itt.hermite_boundary(lamk, x0, B0k, xk, Bk)
          # Then we need to know if we have points on the triangle top
          # See if there are points on the triangle top and if they are associated with a creeping ray
          if( nTop == indCrTop[currentCrTop] ):
               # This means that there is creeping along this triangle top
               # Hence from zkPrev the path goes to ak and creeps to bk
               # which then shoots to yk and then creeps to zk
               etaRegionOutside = listIndices[n + j]
               etaMinCr = min(etaRegionOutside, etakPrev)
               rk = paramsCrTop[2*currentCrTop]
               sk = paramsCrTop[2*currentCrTop + 1]
               ak = itt.hermite_boundary(rk, xkM1, BkBk1_0, xk, BkBk1_1)
               bk = itt.hermite_boundary(sk, xkM1, BkBk1_0, xk, BkBk1_1)
               sum += etakPrev*norm( ak - zkPrev ) # shoots from zkPrev to ak
               sum += etaMinCr*arclengthSimpson(rk, sk, xkM1, BkBk1_0, xk, BkBk1_1) # creeps from ak to bk
               sum += etakPrev*norm( yk - bk ) # shoots from bk to yk
               sum += etaMin*arclengthSimpson(muk, lamk, x0, B0k, xk, Bk) # creeps from yk to zk
               # Update the current index of the creeping updates
               if (currentCrTop  < len(indCrTop) - 1):
                    currentCrTop += 1
          elif( nTop == indStTop[currentStTop]):
               # This means that there is no creeping along this triangle top, it goes
               # straight through that ouside region
               # from zkPrev the path goes to ak, from ak it goes to bk, from bk to yk
               # from yk it creeps to zk
               etaRegionOutside = listIndices[n + j]
               rk = paramsStTop[2*currentStTop]
               sk = paramsStTop[2*currentStTop + 1]
               ak = itt.hermite_boundary(rk, xkM1, BkBk1_0, xk, BkBk1_1)
               bk = itt.hermite_boundary(sk, xkM1, BkBk1_0, xk, BkBk1_1)
               sum += etakPrev*norm( ak - zkPrev ) # shoots from zkPrev to ak
               sum += etaRegionOutside*norm( bk - ak )  # shoots from ak to bk
               sum += etakPrev*norm( yk - bk ) # shoots from bk to yk
               sum += etaMin*arclengthSimpson(muk, lamk, x0, B0k, xk, Bk) # creeps from yk to zk
               # Update the current index of the creeping updates
               if (currentCrTop  < len(indCrTop) - 1):
                    currentCrTop += 1
          else:
               # This means that there are no points along this triangle top, we proceed as "usual"
               sum += etakPrev*norm(yk - zkPrev) + etaMin*arclengthSimpson(muk, lamk, x0, B0k, xk, Bk)
     return sum
          
################################################################################################
################################################################################################
################################################################################################
################################################################################################
#### Knowing when an update is feasible in a generalized triangle


# Project back lamk given mukM1 when there is no creeping or shooting through the side edge

def project_lamkGivenmuk1_noCr(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, BkM1Bk_0, BkM1Bk_1):
     '''
     Project back lamk given mukM1 when there is no creeping or shooting though the side edge xkxk1
     In this case we assume that the edge xkxk1 curves towards the inside of the curvy triangle
     we also assume that lamk in the previous iteration, lamkPrev is such that lamk > lamkPrev
     (otherwise we don't need to project back like this, we would just need a box projection)
     '''
     lamk = project_box(lamk) # Very first step in every projection method in here
     zkM1 = itt.hermite_boundary(mukM1, x0, B0kM1, xkM1, BkM1)
     yk = itt.hermite_boundary(lamk, x0, B0k, xk, Bk)
     B0kM1_mukM1 = itt.gradientBoundary(mukM1, x0, B0kM1, xkM1, BkM1)
     B0k_lamk = itt.gradientBoundary(lamk, x0, B0k, xk, Bk)
     # We need to find a_tan = hkM1hk(r_tan)
     rPass = lambda r: findRtan(r, xkM1, xk, BkM1Bk_0, BkM1Bk_1, zkM1)
     rootTan = root_scalar(rPass, method = "secant", x0 = 0.4, x1 = 0.5)
     r_tan = rootTan.root
     a_tan = itt.hermite_boundary(r_tan, xkM1, BkM1Bk_0, xk, BkM1Bk_1)
     BkM1Bk_tan = itt.gradientBoundary(r_tan, xkM1, BkM1Bk_0, xk, BkM1Bk_1)

     # The normals - depends on the position of x0 and xk  with respect to xkM1
     NkM1_mukM1 = np.array([-B0kM1_mukM1[1], B0kM1_mukM1[0]])
     Nk_lamk = np.array([-B0k_lamk[1], B0k_lamk[0]])
     N_tan = np.array([-BkM1Bk_tan[1], BkM1Bk_tan[0]])
     NkM1Nk_0 = np.array([-BkM1Bk_0[1], BkM1Bk_0[0]])
     if( np.dot(NkM1_mukM1, xk - x0) < 0):
          NkM1_mukM1 = -NkM1_mukM1
     if( np.dot(Nk_lamk, x0 - xkM1) < 0):
          Nk_lamk = -Nk_lamk
     if( np.dot(NkM1Nk_0, x0 - xk) < 0):
          N_tan = -N_tan

     # Tests
     dotTestMin_fromh0kM1 = np.dot( yk - zkM1, NkM1_mukM1 ) # should be positive
     dotTestMax_fromh0kM1 = np.dot( yk - zkM1, Nk_lamk ) # should be positive
     dotTestMax_fromhkM1hk = np.dot( yk - a_tan, N_tan ) # should be positive

     # Test if lamk < lamMin
     if( dotTestMin_fromh0kM1 < 0 ):
          #print(" dotTestMin_fromh0kM1 failed", yk, ",   ", zkM1)
          tMin = lambda lam: t1(lam, x0, xk, B0k, Bk, zkM1, B0kM1_mukM1)
          rootMin = root_scalar(tMin, method = "secant", x0 = 0.4, x1 = 0.5)
          lamk = rootMin.root
          yk = itt.hermite_boundary(lamk, x0, B0k, xk, Bk)
          B0k_lamk = itt.gradientBoundary(lamk, x0, B0k, xk, Bk)
          Nk_lamk = np.array([-B0k_lamk[1], B0k_lamk[0]])
          if( np.dot(Nk_lamk, x0 - xkM1) < 0):
               Nk_lamk = -Nk_lamk
          dotTestMax_fromh0kM1 = np.dot( yk - zkM1, Nk_lamk )
     # Test if lamk > lamkMax (from h0hkM1)
     if( dotTestMax_fromh0kM1 < 0 ):
          tMax = lambda lam: t2(lam, x0, xk, B0k, Bk, zkM1)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          lamk = rootMax.root
          # Update this lambda to test for max from hkM1hk
          rPass = lambda r: findRtan(r, xkM1, xk, BkM1Bk_0, BkM1Bk_1, zkM1)
          rootTan = root_scalar(rPass, method = "secant", x0 = 0.4, x1 = 0.5)
          r_tan = rootTan.root
          a_tan = itt.hermite_boundary(r_tan, xkM1, BkM1Bk_0, xk, BkM1Bk_1)
          BkM1Bk_tan = itt.gradientBoundary(r_tan, xkM1, BkM1Bk_0, xk, BkM1Bk_1)
          N_tan = np.array([-BkM1Bk_tan[1], BkM1Bk_tan[0]])
          if( np.dot(NkM1Nk_0, x0 - xk) < 0):
               N_tan = -N_tan
          dotTestMax_fromhkM1hk = np.dot( yk - a_tan, N_tan )
     # Test if lamk > lamkMax (from hkM1hk)
     if( dotTestMax_fromhkM1hk < 0 ):
          tMax = lambda lam: t1(lam, x0, xk, B0k, Bk, a_tan, BkM1Bk_tan)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          lamk = rootMax.root
     lamk = project_box(lamk)
     return lamk

# Project back muk given lamk1 when there is no creeping or shooting through the side edge

def project_mukGivenlamk1_noCr(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1):
     '''
     Project back muk given lamk1 when there is no creeping or shooting through the side
     edge xkxk1. In this case we assume that the edge xkxk1 curves towards the inside of the
     curvy trianlge. We also assume that muk in the previous iteration, mukPrev is such
     that muk>mukPrev (otherwise we don't need to project back like this, we can just use a box
     projection)
     '''
     muk = project_box(muk)
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
     # We need to find a_tan = hkhk1(r_tan)
     rPass = lambda r: findRtan(r, xk, xk1, BkBk1_0, BkBk1_1, yk1)
     rootTan = root_scalar(rPass, method = "secant", x0 = 0.4, x1 = 0.5)
     r_tan = rootTan.root
     a_tan = itt.hermite_boundary(r_tan, xk, BkBk1_0, xk1, BkBk1_1)
     BkBk1_tan = itt.gradientBoundary(r_tan, xk, BkBk1_0, xk1, BkBk1_1)

     # Compute the normals
     N0k1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
     N0k_muk = np.array([-B0k_muk[1], B0k_muk[0]])
     N_tan = np.array([-BkBk1_tan[1], BkBk1_tan[0]])
     NkNk1_0 = np.array([-BkBk1_0[1], BkBk1_1[0]])
     # Check if this is the correct orientation - check the position of x0, xk, xk1
     if( np.dot( N0k1_lamk1, x0 - xk) < 0):
          N0k1_lamk1 = -N0k1_lamk1
     if( np.dot( N0k_muk, xk1 - x0) < 0):
          N0k_muk = -N0k_muk
     if( np.dot(NkNk1_0, x0 - xk) < 0):
          N_tan = -N_tan

     # Tests
     dotTestMin_fromh0hk =  np.dot(N0k1_lamk1, yk1 - zk) # Should be positive
     dotTestMax_fromh0hk = np.dot(N0k_muk, yk1 - zk) # Should be positive
     dotTestMax_fromhkhk1 = np.dot(N_tan, zk - a_tan) # Should be positive
     
     # Test if muk < mukMin
     if(dotTestMin_fromh0hk < 0 ):
          tMin = lambda mu: t4(mu, x0, xk, B0k, Bk, yk1, B0k1_lamk1)
          rootMin = root_scalar(tMin, bracket = [0,1])
          muk = rootMin.root
          zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
          B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
          N0k_muk = np.array([-B0k_muk[1], B0k_muk[0]])
          if( np.dot( N0k_muk, xk1 - x0) < 0):
               N0k_muk = -N0k_muk
          dotTestMax = np.dot(N0k_muk, yk1 - zk)
     # Test if muk > mukMax (from h0hk)
     if(dotTestMax_fromh0hk < 0 ):
          tMax = lambda mu: t3(mu, x0, xk, B0k, Bk, yk1)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          muk = rootMax.root
          # Update this mu to test for max from hkhk1
          rPass = lambda r: findRtan(r, xk, xk1, BkBk1_0, BkBk1_1, yk1)
          rootTan = root_scalar(rPass, method = "secant", x0 = 0.4, x1 = 0.5)
          r_tan = rootTan.root
          a_tan = itt.hermite_boundary(r_tan, xk, BkBk1_0, xk1, BkBk1_1)
          BkBk1_tan = itt.gradientBoundary(r_tan, xk, BkBk1_0, xk1, BkBk1_1)
          N_tan = np.array([-BkBk1_tan[1], BkBk1_tan[0]])
          if( np.dot(NkNk1_0, x0 - xk) < 0):
               N_tan = -N_tan
          dotTestMax_fromhkhk1 = np.dot(N_tan, zk - a_tan)
     if(dotTestMax_fromhkhk1 < 0):
          tMax = lambda mu: t4(mu, x0, xk, B0k, Bk, a_tan, BkBk1_tan)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          #print( tMax(muk) )
          muk = rootMax.root
          #print( tMax(muk))
     muk = project_box(muk)
     return muk


def project_rkGivenmuk(rk, muk, x0, B0k, xk, Bk, xk1, Bk1, BkBk1_0, BkBk1_1):
     '''
     Project back rk on the side boundary hkhk1 given muk on the bottom
     boundary h0hk.
     '''
     rk = project_box(rk)
     ak = itt.hermite_boundary(rk, xk, BkBk1_0, xk1, BkBk1_1)
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
     BkBk1_rk = itt.gradientBoundary(rk, xk, BkBk1_0, xk1, BkBk1_1)

     # Compute the normals
     N0k_muk = np.array([-B0k_muk[1], B0k_muk[0]])
     NkNk1_rk = np.array([-BkBk1_rk[1], BkBk1_rk[0]])
     NkNk1_0 = np.array([-BkBk1_0[1], BkBk1_1[0]])
     # See if this orientation is correct - depends on the position of x0 and xk with respect to xk1
     if( np.dot(N0k_muk, xk1 - x0) < 0):
          N0k_muk = -N0k_muk
     if( np.dot(NkNk1_0, x0 - xk) < 0):
          NkNk1_rk = -NkNk1_rk

     # Compute the tests
     dotTestMin = np.dot( N0k_muk, ak - zk)
     dotTestMax = np.dot( NkNk1_rk, zk - ak)

     # Test if rk < rMin
     if( dotTestMin < 0):
          tMin = lambda r: t1(r, xk, xk1, BkBk1_0, BkBk1_1, zk, B0k_muk)
          rootMin = root_scalar(tMin, method = "secant", x0 = 0.4, x1 = 0.5)
          rk = rootMin.root
          ak = itt.hermite_boundary(rk, xk, BkBk1_0, xk1, BkBk1_1)
          BkBk1_rk = itt.gradientBoundary(rk, xk, BkBk1_0, xk1, BkBk1_1)
          NkNk1_rk = np.array([-BkBk1_rk[1], BkBk1_rk[0]])
          if( np.dot(NkNk1_0, x0 - xk) < 0):
               NkNk1_rk = -NkNk1_rk
          dotTestMax = np.dot( NkNk1_rk, ak - zk)
     if( dotTestMax < 0):
          tMax = lambda r: t2(r, xk, xk1, BkBk1_0, BkBk1_1, zk)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          rk = rootMax.root
     rk = project_box(rk)
     return rk

def project_skGivenlamk1(sk, lamk1, x0, B0k1, xk1, Bk1, xk, BkBk1_0, BkBk1_1):
     '''
     Project back sk on the side boundary hkhk1 given lamk1 on the top
     boundary h0hk1
     '''
     sk = project_box(sk)
     bk = itt.hermite_boundary(sk, xk, BkBk1_0, xk1, BkBk1_1)
     yk1 = itt.hermite_boundary(lamk1, x0, B0k1, xk1, Bk1)
     B0k1_lamk1 = itt.gradientBoundary(lamk1, x0, B0k1, xk1, Bk1)
     BkBk1_sk = itt.gradientBoundary(sk, xk, BkBk1_0, xk1, BkBk1_1)

     # Compute the normals
     Nk1_lamk1 = np.array([-B0k1_lamk1[1], B0k1_lamk1[0]])
     NkNk1_sk = np.array([-BkBk1_sk[1], BkBk1_sk[0]])
     NkNk1_0 = np.array([-BkBk1_0[1], BkBk1_1[0]])
     # See if this orientation is correct
     if( np.dot( Nk1_lamk1,  x0 - xk) < 0):
          Nk1_lamk1 = -Nk1_lamk1
     if( np.dot( NkNk1_0, x0 - xk) < 0  ):
          NkNk1_sk = -NkNk1_sk

     # Compute the tests
     dotTestMin = np.dot( Nk1_lamk1, yk1 - bk ) # Should be positive
     dotTestMax = np.dot( NkNk1_sk , yk1 - bk ) #Should be positive

     # If  sk < sMin
     if( dotTestMin < 0):
          tMin = lambda s: t4(s, xk, xk1, BkBk1_0, BkBk1_1, yk1, B0k1_lamk1)
          rootMin = root_scalar(tMin, method = "secant", x0 = 0.4, x1 = 0.5)
          sk = rootMin.root
          bk = itt.hermite_boundary(sk, xk, BkBk1_0, xk1, BkBk1_1)
          BkBk1_sk = itt.gradientBoundary(sk, xk, BkBk1_0, xk1, BkBk1_1)
          NkNk1_sk = np.array([-BkBk1_sk[1], BkBk1_sk[0]])
          if( np.dot( NkNk1_0, x0 - xk) < 0  ):
               NkNk1_sk = -NkNk1_sk
          dotTestMax = np.dot( NkNk1_sk , yk1 - bk )
     if( dotTestMax < 0):
          tMax = lambda s: t2(s, xk, xk1, BkBk1_0, BkBk1_1, yk1)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          sk = rootMax.root
     sk = project_box(sk)
     return sk


def project_mukGivenrk(muk, rk, x0, B0k, xk, Bk, BkBk1_0, xk1, BkBk1_1):
     '''
     Project back muk given rk for an update that includes points on the side
     boundary hkhk1
     '''
     muk = project_box(muk)
     ak = itt.hermite_boundary(rk, xk, BkBk1_0, xk1, BkBk1_1)
     zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
     BkBk1_rk = itt.gradientBoundary(rk, xk, BkBk1_0, xk1, BkBk1_1)
     B0k_muk = itt.gradientBoundary(muk, x0, B0k, xk, Bk)
     NkNk1_0 = np.array([-BkBk1_0[1], BkBk1_1[0]])
     # Compute the normals
     NkNk1_rk = np.array([-BkBk1_rk[1], BkBk1_rk[0]])
     N0k_muk = np.array([-B0k_muk[1], B0k_muk[0]])
     # See if this orientation is correct
     if( np.dot( N0k_muk, xk1 - x0 ) < 0):
          N0k_muk = -N0k_muk
     if( np.dot( NkNk1_0, x0 - xk) < 0):
          NkNk1_rk = -NkNk1_rk

     # Compute the tests
     dotTestMin =  np.dot(N0k_muk, ak - zk)
     dotTestMax = np.dot(NkNk1_rk, zk - ak)
     if(dotTestMin<0 ):
          tMin = lambda mu: t2(mu, x0, xk, B0k, Bk, ak)
          rootMin = root_scalar(tMin, method = "secant", x0 = 0.4, x1 = 0.5)
          muk = rootMin.root
          zk = itt.hermite_boundary(muk, x0, B0k, xk, Bk)
          dotTestMax = np.dot(NkNk1_rk, zk - ak)
     if(dotTestMax < 0 ):
          tMax = lambda mu: t1(mu, x0, xk, B0k, Bk, ak, BkBk1_rk)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          muk = rootMax.root
     muk = project_box(muk)
     return muk

def project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1):
     '''
     Project back lamk given sk for an update that includes points on the side
     boundary hkM1hk
     '''
     lamk = project_box(lamk)
     bkM1 = itt.hermite_boundary(skM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1)
     yk = itt.hermite_boundary(lamk, x0, B0k, xk, Bk)
     BkM1Bk_skM1 = itt.gradientBoundary(skM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1)
     B0k_lamk = itt.gradientBoundary(lamk, x0, B0k, xk, Bk)

     # Compute the normals
     NkM1Nk_skM1 = np.array([-BkM1Bk_skM1[1], BkM1Bk_skM1[0]])
     N0k_lamk = np.array([-B0k_lamk[1], B0k_lamk[0]])
     NkM1Nk_0 = np.array([-BkM1Bk_0[1], BkM1Bk_0[0]])
     # See if this orientation is correct
     if( np.dot( N0k_lamk, x0 - xkM1) < 0):
          N0k_lamk = -N0k_lamk
     if( np.dot( NkM1Nk_0, x0 - xkM1) < 0):
          NkM1Nk_skM1 = -NkM1Nk_skM1

     # Compute the tests
     testMin = np.dot(N0k_lamk, yk - bkM1 )
     testMax = np.dot(NkM1Nk_skM1, yk - bkM1)
     # Test if lamk<lamMin
     if(testMin < 0):
          tMin = lambda lam: t2(lam, x0, xk, B0k, Bk, bkM1)
          rootMin = root_scalar(tMin, method = "secant", x0 = 0.4, x1 = 0.5)
          lamk = rootMin.root
          yk = itt.hermite_boundary(lamk, x0, B0k, xk, Bk)
          B0k_lamk = itt.gradientBoundary(lamk, x0, B0k, xk, Bk)
          N0k_lamk = np.array([-B0k_lamk[1], B0k_lamk[0]])
          if( np.dot( N0k_lamk, x0 - xkM1) < 0):
               N0k_lamk = -N0k_lamk
          testMax = np.dot(N0k_lamk, bkM1 - yk)
     if(testMax < 0):
          tMax = lambda lam: t1(lam, x0, xk, B0k, Bk, bkM1, BkM1Bk_skM1)
          rootMax = root_scalar(tMax, method = "secant", x0 = 0.4, x1 = 0.5)
          lamk = rootMax.root
     lamk = project_box(lamk)
     return lamk




############ TWO BY TWO PROJECTIONS


def projections_muk_lamk1(muk_candidate, lamk1_candidate, muk_free, lamk1_free, k, params, x0,
                          T0, grad0, x1, T1, grad1, xHat, listIndices,
                          listxk, listB0k, listBk, listBkBk1,
                          indCrTop, paramsCrTop, indStTop, paramsStTop, cl = 1e-1):
     '''
     Two options depending on which one is better:     - project muk given lamk1
                                                       - project lamk1 given muk
     '''
     if( muk_candidate == muk_free and lamk1_candidate == lamk1_free):
          return muk_candidate, lamk1_free
     params_muk_projected = np.copy(params)
     params_muk_projected[k] = muk_candidate
     if( abs( params[k-1] - params[k-2] ) < cl):
          params_muk_projected[k-1] = muk_candidate
     params_muk_projected[k+1] = lamk1_free
     if( abs( params[k+1] - params[k]) < cl and k < len(params) - 1):
          params_muk_projected[k+2] = lamk1_free
     f_muk_projected = fObj_generalized(params_muk_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     params_lamk1_projected = np.copy(params)
     params_lamk1_projected[k] = muk_free
     if( abs( params[k-1] - params[k-2]) < cl):
          params_lamk1_projected[k-1] = muk_free
     params_lamk1_projected[k+1] = lamk1_candidate
     if( abs(params[k+1] - params[k]) < cl and k < len(params) -1):
          params_lamk1_projected[k] = lamk1_candidate
     f_lamk1_projected = fObj_generalized(params_lamk1_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop)
     if( f_muk_projected < f_lamk1_projected):
          return muk_candidate, lamk1_free
     else:
          return muk_free, lamk1_candidate


def projections_muk_rkCr(muk_candidate, rk_candidate, muk_free, rk_free, k, kCrTop, params, x0,
                          T0, grad0, x1, T1, grad1, xHat, listIndices,
                          listxk, listB0k, listBk, listBkBk1,
                          indCrTop, paramsCrTop, indStTop, paramsStTop, cl = 1e-1):
     '''
     FOR CREEPING ALONG THE k-th side
     Two options depending on which one is better:     - project muk given rk
                                                       - project rk given muk
     '''
     if( muk_candidate == muk_free and rk_candidate == rk_free):
          return muk_candidate, rk_free
     params_muk_projected = np.copy(params)
     params_muk_projected[k] = muk_candidate
     if( abs( params[k-1] - params[k]) < cl ):
          params_muk_projected[k-1] = muk_candidate
     paramsCrTop_muk_projected = np.copy(paramsCrTop)
     paramsCrTop_muk_projected[kCrTop] = rk_free
     if( abs( paramsCrTop[kCrTop + 1] -  paramsCrTop[kCrTop ]) < cl):
          paramsCrTop_muk_projected[kCrTop + 1] = rk_free
     f_muk_projected = fObj_generalized(params_muk_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop_muk_projected, indStTop, paramsStTop)
     params_rk_projected = np.copy(params)
     params_rk_projected[k] = muk_free
     if( abs( params[k-1] - params[k]) < cl):
          params_rk_projected[k-1] = muk_free
     paramsCrTop_rk_projected = np.copy(paramsCrTop)
     paramsCrTop_rk_projected[kCrTop] = rk_candidate
     if( abs( paramsCrTop[kCrTop + 1] - paramsCrTop[kCrTop] ) < cl ):
          paramsCrTop_rk_projected[kCrTop + 1] = rk_candidate
     f_rk_projected = fObj_generalized(params_rk_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop_rk_projected, indStTop, paramsStTop)
     if( f_muk_projected < f_rk_projected ):
          return muk_candidate, rk_free
     else:
          return muk_free, rk_candidate


def projections_muk_rkSt(muk_candidate, rk_candidate, muk_free, rk_free, k, kStTop, params, x0,
                          T0, grad0, x1, T1, grad1, xHat, listIndices,
                          listxk, listB0k, listBk, listBkBk1,
                          indCrTop, paramsCrTop, indStTop, paramsStTop, cl = 1e-1):
     '''
     FOR STRAIGHT LINE THROUGH THE k-th side
     Two options depending on which one is better:     - project muk given rk
                                                       - project rk given muk
     '''
     if( muk_candidate == muk_free and rk_candidate == rk_free):
          return muk_candidate, rk_free
     params_muk_projected = np.copy(params)
     params_muk_projected[k] = muk_candidate
     if( abs( params[k-1] - params[k]) < cl ):
          params_muk_projected[k-1] = muk_candidate
     paramsStTop_muk_projected = np.copy(paramsStTop)
     paramsStTop_muk_projected[kStTop] = rk_free
     if( abs( paramsStTop[kStTop + 1] - paramsStTop[kStTop]) < cl):
          paramsStTop[kStTop + 1] = rk_free
     f_muk_projected = fObj_generalized(params_muk_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop, indStTop, paramsStTop_muk_projected)
     params_rk_projected = np.copy(params)
     params_rk_projected[k] = muk_free
     if( abs( params[k-1] - params[k]) < cl):
          params_rk_projected[k-1] = muk_free
     paramsStTop_rk_projected = np.copy(paramsStTop)
     paramsStTop_rk_projected[kStTop] = rk_candidate
     if( abs( paramsStTop[kStTop + 1] - paramsStTop[kStTop]) < cl ):
          paramsStTop_rk_projected[kStTop + 1] = rk_candidate
     f_rk_projected = fObj_generalized(params_rk_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop, indStTop, paramsStTop_rk_projected)
     if( f_muk_projected < f_rk_projected ):
          return muk_candidate, rk_free
     else:
          return muk_free, rk_candidate


def projections_skCr_lamk1(sk_candidate, lamk1_candidate, sk_free, lamk1_free, k, kCrTop, params, x0,
                          T0, grad0, x1, T1, grad1, xHat, listIndices,
                          listxk, listB0k, listBk, listBkBk1,
                           indCrTop, paramsCrTop, indStTop, paramsStTop, cl = 1e-1):
     '''
     FOR CREEPING ALONG THE k-th side
     Two options depending on which is better:   - project sk given lamk1
                                                 - project lamk1 given sk
     '''
     if( sk_candidate == sk_free and lamk1_candidate == lamk1_free):
          return sk_candidate, lamk1_free
     paramsCrTop_sk_projected = np.copy(paramsCrTop)
     paramsCrTop_sk_projected[kCrTop + 1] = sk_candidate
     if( abs( paramsCrTop[kCrTop + 1] - paramsCrTop[kCrTop] ) < cl):
          # If rk = sk numerically
          paramsCrTop_sk_projected[kCrTop] = sk_candidate
     params_sk_projected = np.copy(params)
     params_sk_projected[k] = lamk1_free
     if( abs( params[k] - params[k+1]) < cl and k < len(params) - 1):
          # If lamk1 == muk1 numerically
          params_sk_projected[k+1] = lamk1_free
     f_sk_projected = fObj_generalized(params_sk_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop_sk_projected, indStTop, paramsStTop)
     paramsCrTop_lamk1_projected = np.copy(paramsCrTop)
     paramsCrTop_lamk1_projected[kCrTop + 1] = sk_free
     if( abs(paramsCrTop[kCrTop] - paramsCrTop[kCrTop + 1]) < cl):
          # If rk = sk numerically
          paramsCrTop_lamk1_projected[kCrTop] = sk_free
     params_lamk1_projected = np.copy(params)
     params_lamk1_projected[k] = lamk1_candidate
     if( abs(params[k+1] - params[k] ) < cl and k < len(params) - 1):
          # Then lamk1 == muk1 numerically
          params_lamk1_projected[k+1] = lamk1_candidate
     f_lamk1_projected = fObj_generalized(params_lamk1_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop_lamk1_projected, indStTop, paramsStTop)
     if( f_sk_projected < f_lamk1_projected):
          return sk_candidate, lamk1_free
     else:
          return sk_free, lamk1_candidate


def projections_skSt_lamk1(sk_candidate, lamk1_candidate, sk_free, lamk1_free, k, kStTop, params, x0,
                          T0, grad0, x1, T1, grad1, xHat, listIndices,
                          listxk, listB0k, listBk, listBkBk1,
                           indCrTop, paramsCrTop, indStTop, paramsStTop, cl = 1e-1):
     '''
     FOR STRAIGHT LINE THROUGH THE k-th side
     Two options depending on which is better:   - project sk given lamk1
                                                 - project lamk1 given sk
     '''
     if( sk_candidate == sk_free and lamk1_candidate == lamk1_free):
          return sk_candidate, lamk1_free
     paramsStTop_sk_projected = np.copy(paramsStTop)
     paramsStTop_sk_projected[kStTop + 1] = sk_candidate
     if( abs(paramsStTop[kStTop + 1] - paramsStTop[kStTop] ) < cl):
          # If rk = sk numerically
          paramsStTop_sk_projected[kStTop] = sk_candidate
     params_sk_projected = np.copy(params)
     params_sk_projected[k] = lamk1_free
     if( abs( params[k] - params[k+1]) < cl and k < len(params) - 1):
          # If lamk1 == muk1 numerically
          params_sk_projected[k+1] = lamk1_free
     f_sk_projected = fObj_generalized(params_sk_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop, indStTop, paramsStTop_sk_projected)
     paramsStTop_lamk1_projected = np.copy(params)
     paramsStTop_lamk1_projected[kStTop + 1] = sk_free
     if( abs(paramsStTop[kStTop] - paramsStTop[kStTop + 1]) < cl):
          # If rk = sk numerically
          paramsStTop[kStTop] = sk_free
     params_lamk1_projected = np.copy(params)
     params_lamk1_projected[k] = lamk1_candidate
     if( abs(params[k+1] - params[k] ) < cl and k < len(params) - 1):
          # Then lamk1 == muk1 numerically
          params_lamk1_projected[k+1] = lamk1_candidate
     f_lamk1_projected = fObj_generalized(params_lamk1_projected, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop_lamk1_projected)
     if( f_sk_projected < f_lamk1_projected):
          return sk_candidate, lamk1_free
     else:
          return sk_free, lamk1_candidate


################################################################################################
################################################################################################
################################################################################################
################################################################################################
#### Optimization for the generalized objective function


###################################
# General backtracking (mainly just for mu1)


def backTr_coord(alpha0, k, d, params, x0, T0, grad0, x1, T1, grad1, xHat,
                 listIndices, listxk, listB0k, listBk, listBkBk1,
                 indCrTop, paramsCrTop, indStTop, paramsStTop):
     '''
     Backtracking for just one coordinate (the k-th coordinate) in params (not for
     the points on the side of the triangle fan (tops)
     '''
     f_before = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     i = 0
     params_test = np.copy(params)
     alpha = alpha0*1/(max(abs(d), 1))
     params_test[k] = params[k] - alpha*d
     f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop)
     # If there is a decrease in the function, try increasing alpha, the step size
     while( (f_test < f_before) and i < 8 ):
          alpha = alpha*1.3
          params_test[k] = params[k] - alpha*d
          f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          i += 1
     i = 0
     # If there is no decrease in the function, try decreasing alpha, the step size
     while( (f_before <= f_test) and i < 25 ):
          alpha = alpha*0.2
          params_test[k] = params[k] - alpha*d
          f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          i += 1
     if( f_before <= f_test ):
          return 0
     else:
          return alpha


###################################
# Backtracking for updates close to the identity


def backTrClose_block0k(alpha0, k, dlamk, dmuk, dCollapsed,
                        params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop):
     '''
     Backtracking to find the next block [lamk, muk] in the generalized objective function (considers points on the sides
     of the triangle fan)
     it considers three directions: steepest descent
                                    steepest descent projected onto the line lamk = muk
                                    steepest descent for the collapsed point (i.e. muk =∼ lamk)
     '''
     f_before = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     i = 0
     d_middle = np.array([dlamk, dmuk])
     params_test = np.copy(params)
     alpha = alpha0*1/(max(norm(d_middle), 1))
     params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
     params_test_proj = np.copy(params)
     params_test_proj[k:(k+2)] = project_ontoLine(params_test[k:(k+2)])
     params_collapsed = np.copy(params)
     params_collapsed[k:(k+2)] = params[k] - alpha*dCollapsed
     # Compare the function value
     f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop)
     f_test_proj = fObj_generalized(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
     f_test_collapsed = fObj_generalized(params_collapsed, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
     # If there is a decrease in the function, try increasing alpha, the step size
     while( ( (f_test < f_before) or (f_test_proj < f_before) or (f_test_collapsed < f_before) ) and i < 8 ):
          alpha = alpha*1.3 # increase the step size
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          params_test_proj[k:(k+2)] = project_ontoLine(params_test[k:(k+2)])
          params_collapsed[k:(k+2)] = params[k] - alpha*dCollapsed
          f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_test_proj = fObj_generalized(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_test_collapsed = fObj_generalized(params_collapsed, x0, T0, grad0, x1, T1, grad1, xHat,
                                              listIndices, listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          i += 1
     i = 0
     # Then if there is no decrease try decreasing alpha, the step size
     while( ( (f_before < f_test) and (f_before < f_test_proj) and (f_before < f_test_collapsed) ) and i < 25 ):
          alpha = alpha*0.2 # increase the step size
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          params_test_proj[k:(k+2)] = project_ontoLine(params_test[k:(k+2)])
          params_collapsed[k:(k+2)] = params[k] - alpha*dCollapsed
          #breakpoint()
          f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_test_proj = fObj_generalized(params_test_proj, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_test_collapsed = fObj_generalized(params_collapsed, x0, T0, grad0, x1, T1, grad1, xHat,
                                              listIndices, listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          i += 1
     # Now we should have a decrease or set alpha to 0
     if( f_before <= f_test and f_before <= f_test_proj and f_before <= f_test_collapsed):
          return params[k], params[k+1]
     elif( f_test < f_test_proj and f_test < f_test_collapsed ):
          return params_test[k], params_test[k+1]
     elif( f_test_proj < f_test_collapsed ):
          return params_test_proj[k], params_test_proj[k+1]
     else:
          return params_collapsed[k], params_collapsed[k+1]

def backTrClose_blockCrTop(alpha0, kCrTop, drk, dsk, dCollapsed,
                           params, x0, T0, grad0, x1, T1, grad1, xHat,
                           listIndices, listxk, listB0k, listBk, listBkBk1,
                           indCrTop, paramsCrTop, indStTop, paramsStTop):
     '''
     Backtracking to find the next block [rk, sk] in the generalized objective functions (considers points
     on the sides of the triangle fan)
     it considers three directions: steepest descent
                                    steepest descent projected onto the line rk = sk
                                    steepest descent for the collapsed point (i.e. rk =∼ sk)
     '''
     f_before = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     i = 0
     d_middle = np.array([drk, dsk])
     #breakpoint()
     paramsCrTop_test = np.copy(paramsCrTop)
     paramsCrTop_test_proj = np.copy(paramsCrTop)
     paramsCrTop_collapsed = np.copy(paramsCrTop)
     alpha = alpha0*1/(max(norm(d_middle), 1))
     paramsCrTop_test[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop:(kCrTop+2)] - alpha*d_middle
     paramsCrTop_test_proj[kCrTop:(kCrTop + 2)] = project_ontoLine(paramsCrTop_test[kCrTop:(kCrTop + 2)])
     paramsCrTop_collapsed[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop] - alpha*dCollapsed
     # Compare the function value
     f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop_test, indStTop, paramsStTop)
     f_test_proj = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop_test_proj, indStTop, paramsStTop)
     f_test_collapsed = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop_collapsed, indStTop, paramsStTop)
     # If there is a decrease in the function, try increasing alpha, the step size
     while( ( (f_test < f_before) or (f_test_proj < f_before) or (f_test_collapsed < f_before ) ) and i < 8 ):
          alpha = alpha*1.3
          paramsCrTop_test[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop:(kCrTop+2)] - alpha*d_middle
          paramsCrTop_test_proj[kCrTop:(kCrTop + 2)] = project_ontoLine(paramsCrTop_test[kCrTop:(kCrTop + 2)])
          paramsCrTop_collapsed[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop] - alpha*dCollapsed
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop_test, indStTop, paramsStTop)
          f_test_proj = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop_test_proj, indStTop, paramsStTop)
          f_test_collapsed = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                              listIndices, listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop_collapsed, indStTop, paramsStTop)
          i += 1
     i = 0
     # If there is no decrease in the function value, try decreasing alpha, the step size
     while( ((f_before <= f_test) and (f_before <= f_test_proj) and (f_before < f_test_collapsed) ) and i < 25):
          alpha = alpha*0.2
          paramsCrTop_test[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop:(kCrTop+2)] - alpha*d_middle
          paramsCrTop_test_proj[kCrTop:(kCrTop + 2)] = project_ontoLine(paramsCrTop_test[kCrTop:(kCrTop + 2)])
          paramsCrTop_collapsed[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop] - alpha*dCollapsed
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop_test, indStTop, paramsStTop)
          f_test_proj = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop_test_proj, indStTop, paramsStTop)
          f_test_collapsed = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                              listIndices, listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop_collapsed, indStTop, paramsStTop)
          i += 1
     # Now we should have a decrease or set alpha to 0
     if( f_before <= f_test and f_before <= f_test_proj and f_before <= f_test_collapsed):
          return paramsCrTop[kCrTop], paramsCrTop[kCrTop+1]
     elif( f_test < f_test_proj and f_test < f_test_collapsed ):
          return paramsCrTop[kCrTop], paramsCrTop[kCrTop+1]
     elif( f_test_proj < f_test_collapsed):
          return paramsCrTop_test_proj[kCrTop], paramsCrTop_test_proj[kCrTop+1]
     else:
          return paramsCrTop_collapsed[kCrTop], paramsCrTop_collapsed[kCrTop + 1]


def backTrClose_blockStTop(alpha0, kStTop, drk, dsk, dCollapsed,
                           params, x0, T0, grad0, x1, T1, grad1, xHat,
                           listIndices, listxk, listB0k, listBk, listBkBk1,
                           indCrTop, paramsCrTop, indStTop, paramsStTop):
     '''
     Backtracking to find the next block [rk, sk] in the generalized objective functions (considers points
     on the sides of the triangle fan)
     it considers three directions: steepest descent
                                    steepest descent projected onto the line rk = sk
                                    steepest descent for the collapsed point (i.e. rk =∼ sk)
     '''
     f_before = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     i = 0
     d_middle = np.array([drk, dsk])
     paramsStTop_test = np.copy(paramsStTop)
     paramsStTop_test_proj = np.copy(paramsStTop)
     paramsStTop_collapsed = np.copy(paramsStTop)
     alpha = alpha0*1/(max(norm(d_middle), 1))
     paramsStTop_test[kStTop:(kStTop + 2)] = paramsStTop[kStTop:(kStTop+2)] - alpha*d_middle
     paramsStTop_test_proj[kStTop:(kStTop + 2)] = project_ontoLine(paramsStTop_test[kStTop:(kStTop + 2)])
     paramsStTop_collapsed[kStTop:(kStTop + 2)] = paramsStTop[kStTop] - alpha*dCollapsed
     # Compare the function value
     f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop_test)
     f_test_proj = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop_test_proj)
     f_test_collaped = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                        listIndices, listxk, listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop, indStTop, paramsStTop_collapsed)
     # If there is a decrease in the function, try increasing alpha, the step size
     while( ( (f_test < f_before) or (f_test_proj < f_before) or (f_test_collapsed < f_before )) and i < 8 ):
          alpha = alpha*1.3
          paramsStTop_test[kStTop:(kStTop + 2)] = paramsStTop[kStTop:(kStTop+2)] - alpha*d_middle
          paramsStTop_test_proj[kStTop:(kStTop + 2)] = project_ontoLine(paramsStTop_test[kStTop:(kStTop + 2)])
          paramsStTop_collapsed[kStTop:(kStTop + 2)] = paramsStTop[kStTop] - alpha*dCollapsed
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop_test)
          f_test_proj = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop_test_proj)
          f_test_collaped = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                             listIndices, listxk, listB0k, listBk, listBkBk1,
                                             indCrTop, paramsCrTop, indStTop, paramsStTop_collapsed)
          i += 1
     i = 0
     # If there is no decrease in the function value, try decreasing alpha, the step size
     while( ( (f_before <= f_test) and (f_before <= f_test_proj) and (f_before <= f_test_collapsed ) ) and i < 25):
          alpha = alpha*0.2
          paramsStTop_test[kStTop:(kStTop + 2)] = paramsStTop[kStTop:(kStTop+2)] - alpha*d_middle
          paramsStTop_test_proj[kStTop:(kStTop + 2)] = project_ontoLine(paramsStTop_test[kStTop:(kStTop + 2)])
          paramsStTop_collapsed[kStTop:(kStTop + 2)] = paramsStTop[kStTop] - alpha*dCollapsed
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop_test)
          f_test_proj = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop_test_proj)
          f_test_collaped = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                             listIndices, listxk, listB0k, listBk, listBkBk1,
                                             indCrTop, paramsCrTop, indStTop, paramsStTop_collapsed)
          i += 1
     # Now we should have a decrease or set alpha to 0
     if( f_before <= f_test and f_before <= f_test_proj and f_before <= f_test_collapsed):
          return paramsStTop[kStTop], paramsStTop[kStTop+1]
     elif( f_test < f_test_proj and f_test < f_test_collapsed ):
          return paramsStTop_test[kStTop], paramsStTop_test[kStTop+1]
     elif( f_test_proj < f_test_collapsed):
          return paramsStTop_test_proj[kStTop], paramsStTop_test_proj[kStTop+1]
     else:
          return paramsStTop_collapsed[kStTop], paramsStTop_collapsed[kStTop + 1]

###################################
# Backtracking for updates far from the identity


def backTr_block0k(alpha0, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                   listIndices, listxk, listB0k, listBk, listBkBk1,
                   indCrTop, paramsCrTop, indStTop, paramsStTop):
     '''
     Backtracking to find the next block [lamk, muk] in the generalized objective function (considers points on the sides
     of the triangle fan)
     it considers one direction: steepest descent
     '''
     f_before = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     i = 0
     d_middle = np.array([dlamk, dmuk])
     params_test = np.copy(params)
     params_lamk = np.copy(params)
     params_muk = np.copy(params)
     alpha = alpha0*1/(max(norm(d_middle), 1))
     alpha_lamk = 1
     alpha_muk = 1
     params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
     params_lamk[k] = params[k] - alpha*dlamk
     params_muk[k+1] = params[k+1] - alpha*dmuk
     # Compare the function value
     f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop)
     f_lamk = fObj_generalized(params_lamk, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop)
     f_muk = fObj_generalized(params_muk, x0, T0, grad0, x1, T1, grad1, xHat,
                              listIndices, listxk, listB0k, listBk, listBkBk1,
                              indCrTop, paramsCrTop, indStTop, paramsStTop)
     # If there is a decrease in the function, try increasing alpha, the step size
     while( (f_test < f_before or f_lamk < f_before or f_muk < f_before) and i < 8 ):
          alpha = alpha*1.3 # increase the step size
          alpha_lamk = alpha_lamk*1.3
          alpha_muk = alpha_muk*1.3
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          params_lamk[k] = params[k] - alpha*dlamk
          params_muk[k+1] = params[k+1] - alpha*dmuk
          f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_lamk = fObj_generalized(params_lamk, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_muk = fObj_generalized(params_muk, x0, T0, grad0, x1, T1, grad1, xHat,
                                   listIndices, listxk, listB0k, listBk, listBkBk1,
                                   indCrTop, paramsCrTop, indStTop, paramsStTop)
          i += 1
     i = 0
     # Then if there is no decrease try decreasing alpha, the step size
     while( (f_before < f_test and f_before < f_lamk and f_before < f_muk)  and i < 25 ):
          alpha = alpha*0.2 # decrease the step size
          alpha_lamk = alpha_lamk*0.2
          alpha_muk = alpha_muk*0.2
          params_test[k:(k+2)] = params[k:(k+2)] - alpha*d_middle
          params_lamk[k] = params[k] - alpha*dlamk
          params_muk[k+1] = params[k+1] - alpha*dmuk
          f_test = fObj_generalized(params_test, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_lamk = fObj_generalized(params_lamk, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop)
          f_muk = fObj_generalized(params_muk, x0, T0, grad0, x1, T1, grad1, xHat,
                                   listIndices, listxk, listB0k, listBk, listBkBk1,
                                   indCrTop, paramsCrTop, indStTop, paramsStTop)
          i += 1
     # Now we should have a decrease or set alpha to 0
     if( f_lamk < f_before and f_lamk < f_test and f_lamk < f_muk):
          return params_lamk[k], params_lamk[k+1]
     elif( f_muk < f_before and f_muk < f_test ):
          return params_muk[k], params_muk[k+1]
     elif( f_before <= f_test ):
          return params[k], params[k+1]
     else:
          return params_test[k], params_test[k+1]

def backTr_blockCrTop(alpha0, kCrTop, drk, dsk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                      listIndices, listxk, listB0k, listBk, listBkBk1,
                      indCrTop, paramsCrTop, indStTop, paramsStTop):
     '''
     Backtracking to find the next block [rk, sk] in the generalized objective functions (considers points
     on the sides of the triangle fan)
     it considers one direction: steepest descent
     '''
     f_before = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     i = 0
     d_middle = np.array([drk, dsk])
     paramsCrTop_test = np.copy(paramsCrTop)
     paramsCrTop_rk = np.copy(paramsCrTop)
     paramsCrTop_sk = np.copy(paramsCrTop)
     alpha = alpha0*1/(max(norm(d_middle), 1))
     alpha_rk = 1
     alpha_sk = 1
     paramsCrTop_test[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop:(kCrTop+2)] - alpha*d_middle
     paramsCrTop_rk[kCrTop] = paramsCrTop[kCrTop] - alpha*drk
     paramsCrTop_sk[kCrTop + 1] = paramsCrTop[kCrTop + 1] - alpha*dsk
     # Compare the function value
     f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop_test, indStTop, paramsStTop)
     f_rk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                             listIndices, listxk, listB0k, listBk, listBkBk1,
                             indCrTop, paramsCrTop_rk, indStTop, paramsStTop)
     f_sk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                             listIndices, listxk, listB0k, listBk, listBkBk1,
                             indCrTop, paramsCrTop_sk, indStTop, paramsStTop)
     # If there is a decrease in the function, try increasing alpha, the step size
     while( (f_test < f_before or f_rk < f_before or f_sk < f_before) and i < 8 ):
          alpha = alpha*1.3
          alpha_rk = alpha_rk*1.3
          alpha_sk = alpha_sk*1.3
          paramsCrTop_test[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop:(kCrTop+2)] - alpha*d_middle
          paramsCrTop_rk[kCrTop] = paramsCrTop[kCrTop] - alpha*drk
          paramsCrTop_sk[kCrTop + 1] = paramsCrTop[kCrTop + 1] - alpha*dsk
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop_test, indStTop, paramsStTop)
          f_rk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop_rk, indStTop, paramsStTop)
          f_sk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop_sk, indStTop, paramsStTop)
          i += 1
     i = 0
     # If there is no decrease in the function value, try decreasing alpha, the step size
     while( (f_before <= f_test and f_before <= f_rk and f_before <= f_sk) and i < 25):
          alpha = alpha*0.2
          alpha_rk = alpha_rk*0.2
          alpha_sk = alpha_sk*0.2
          paramsCrTop_test[kCrTop:(kCrTop + 2)] = paramsCrTop[kCrTop:(kCrTop+2)] - alpha*d_middle
          paramsCrTop_rk[kCrTop] = paramsCrTop[kCrTop] - alpha*drk
          paramsCrTop_sk[kCrTop + 1] = paramsCrTop[kCrTop + 1] - alpha*dsk
          #breakpoint()
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop_test, indStTop, paramsStTop)
          f_rk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop_rk, indStTop, paramsStTop)
          f_sk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop_sk, indStTop, paramsStTop)
          i += 1
     # Now we should have a decrease or set alpha to 0
     if( f_before <= f_test and f_before <= f_rk and f_before <= f_sk):
          return paramsCrTop[kCrTop], paramsCrTop[kCrTop+1]
     elif( f_rk < f_before and f_rk < f_test and f_rk < f_sk):
          return paramsCrTop_rk[kCrTop], paramsCrTop_rk[kCrTop + 1]
     elif( f_sk < f_before and f_sk < f_test):
          return paramsCrTop_sk[kCrTop], paramsCrTop_sk[kCrTop + 1]
     else:
          return paramsCrTop_test[kCrTop], paramsCrTop_test[kCrTop+1]


def backTr_blockStTop(alpha0, kStTop, drk, dsk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                      listIndices, listxk, listB0k, listBk, listBkBk1,
                      indCrTop, paramsCrTop, indStTop, paramsStTop):
     '''
     Backtracking to find the next block [rk, sk] in the generalized objective functions (considers points
     on the sides of the triangle fan)
     it considers one direction: steepest descent
     '''
     f_before = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                        listIndices, listxk, listB0k, listBk, listBkBk1,
                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     i = 0
     d_middle = np.array([drk, dsk])
     paramsStTop_test = np.copy(paramsStTop)
     paramsStTop_rk = np.copy(paramsStTop)
     paramsStTop_sk = np.copy(paramsStTop)
     alpha = alpha0*1/(max(norm(d_middle), 1))
     alpha_rk = 1
     alpha_sk = 1
     paramsStTop_test[kStTop:(kStTop + 2)] = paramsStTop[kStTop:(kStTop+2)] - alpha*d_middle
     paramsStTop_rk[kStTop] = paramsStTop[kStTop] - alpha*drk
     paramsStTop_sk[kStTop + 1] = paramsStTop[kStTop + 1] - alpha*dsk
     # Compare the function value
     f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop_test)
     f_rk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                             listIndices, listxk, listB0k, listBk, listBkBk1,
                             indCrTop, paramsCrTop, indStTop, paramsStTop_rk)
     f_sk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                             listIndices, listxk, listB0k, listBk, listBkBk1,
                             indCrTop, paramsCrTop, indStTop, paramsStTop_sk)
     # If there is a decrease in the function, try increasing alpha, the step size
     while( (f_test < f_before or f_rk < f_before or f_sk < f_before) and i < 8 ):
          alpha = alpha*1.3
          alpha_rk = alpha_rk*1.3
          alpha_sk = alpha_sk*1.3
          paramsStTop_test[kStTop:(kStTop + 2)] = paramsStTop[kStTop:(kStTop+2)] - alpha*d_middle
          paramsStTop_rk[kStTop] = paramsStTop[kStTop] - alpha*drk
          paramsStTop_sk[kStTop + 1] = paramsStTop[kStTop + 1] - alpha*dsk
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop_test)
          f_rk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop, indStTop, paramsStTop_rk)
          f_sk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop, indStTop, paramsStTop_sk)
          i += 1
     i = 0
     # If there is no decrease in the function value, try decreasing alpha, the step size
     while( (f_before <= f_test and f_before <= f_rk and f_before <= f_sk) and i < 25):
          alpha = alpha*0.2
          alpha_rk = alpha_rk*0.2
          alpha_sk = alpha_sk*0.2
          paramsStTop_test[kStTop:(kStTop + 2)] = paramsStTop[kStTop:(kStTop+2)] - alpha*d_middle
          paramsStTop_rk[kStTop] = paramsStTop[kStTop] - alpha*drk
          paramsStTop_sk[kStTop + 1] = paramsStTop[kStTop + 1] - alpha*dsk
          f_test = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                    listIndices, listxk, listB0k, listBk, listBkBk1,
                                    indCrTop, paramsCrTop, indStTop, paramsStTop_test)
          f_rk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop, indStTop, paramsStTop_rk)
          f_sk = fObj_generalized(params, x0, T0, grad0, x1, T1, grad1, xHat,
                                  listIndices, listxk, listB0k, listBk, listBkBk1,
                                  indCrTop, paramsCrTop, indStTop, paramsStTop_sk)
          i += 1
     # Now we should have a decrease or set alpha to 0
     if( f_before <= f_test and f_before <= f_rk and f_before <= f_sk):
          return paramsStTop[kStTop], paramsStTop[kStTop+1]
     elif( f_rk < f_before and f_rk < f_test and f_rk < f_sk):
          return paramsStTop_rk[kStTop], paramsStTop_rk[kStTop + 1]
     elif( f_sk < f_before and f_sk < f_test):
          return paramsStTop_sk[kStTop], paramsStTop_sk[kStTop + 1]
     else:
          return paramsStTop_test[kStTop], paramsStTop_test[kStTop+1]

###################################
# Auxiliary functions for the forward pass update



def updateFromCrTop(n, j, currentCrTop, currentStTop, params, gammas, theta_gamma, x0, T0,
                    grad0, x1, T1, grad1, xHat, listIndices, listxk,
                    listB0k, listBk, listBkBk1, indCrTop, paramsCrTop,
                    indStTop, paramsStTop, listCurvingInwards, gradParams, gradCrTop):
     '''
     Update [rkM1, skM1] and [lamk, muk] where [rkM1, skM1] define
     a creeping ray from the top of the triangle fan (or side edge,
     however you want to look at it)
     '''
     k = 2*j - 1
     nTop = j
     gamma = gammas[j-1]
     mukM1 = params[k-1]
     lamk = params[k]
     muk = params[k+1]
     B0kM1 = listB0k[j-1]
     B0k = listB0k[j]
     xkM1 = listxk[j]
     xk = listxk[j+1]
     if( j < n ):
          xk1 = listxk[j+2] # Because otherwise this is not defined
          B0k1 = listB0k[j+1]
          Bk = listBk[j+1]
     BkM1 = listBk[j-1]
     Bk = listBk[j]
     BkM1Bk_0 = listBkBk1[k-1] # grad of hkhk1 at xk
     BkM1Bk_1 = listBkBk1[k] # grad of hkhk1 at xk1
     etakPrev = listIndices[j-1] # index of refraction from previous triangle
     etak = listIndices[j] # index of refraction in current triangle
     etaRegionOutside = listIndices[n+j] # index of refraction ouside from the side of the previous triangle
     # This is where the update begins
     etaMinCr = min(etaRegionOutside, etakPrev)
     kCrTop = 2*currentCrTop
     rkM1 = paramsCrTop[kCrTop]
     skM1 = paramsCrTop[kCrTop + 1]
     akM1 = itt.hermite_boundary(rkM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1) # receiver from mukM1
     bkM1 = itt.hermite_boundary(skM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1) # shooter to lamk
     # Compute directions
     drkM1 = partial_fObj_recCr(mukM1, skM1, rkM1, x0, B0kM1, xkM1, BkM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1, etakPrev, etaRegionOutside)
     gradCrTop[kCrTop] = drkM1
     dskM1 = partial_fObj_shCr(skM1, rkM1, lamk, xkM1, BkM1Bk_0, xk, BkM1Bk_1, x0, B0k, xk, Bk, etakPrev, etaRegionOutside)
     gradCrTop[kCrTop + 1] = dskM1     # Decide if we need close or far backtracking
     r = close_to_identity(rkM1, skM1)
     if( r <= gamma ):
          # Meaning we have to do a close update
          # Compute the "collapsed direction"
          drkM1_collapsed = partial_fObj_collapsedShooter(mukM1, rkM1, lamk, x0, B0kM1, xkM1, BkM1,
                                                          xkM1, BkM1Bk_0, xk, BkM1Bk_1,
                                                          x0, B0k, xk, Bk, etakPrev, etak)
          rkM1, skM1 = backTrClose_blockCrTop(1, kCrTop, drkM1, dskM1, drkM1_collapsed,
                                              params, x0, T0, grad0, x1, T1, grad1, xHat,
                                              listIndices, listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
     else:
          # Meaning we dont have to do a close update
          rkM1, skM1 = backTr_blockCrTop(1, kCrTop, drkM1, dskM1, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
     # Project back
     rkM1_projected = project_rkGivenmuk(rkM1, mukM1, x0, B0kM1, xkM1, BkM1, xk, Bk, BkM1Bk_0, BkM1Bk_1)
     mukM1_projected = project_mukGivenrk(mukM1, rkM1, x0, B0kM1, xkM1, BkM1, BkM1Bk_0, xk, BkM1Bk_1)
     rkM1_free = project_box(rkM1)
     mukM1_free = mukM1
     mukM1, rkM1 = projections_muk_rkCr(mukM1_projected, rkM1_projected, mukM1_free, rkM1_free, k-1,
                                        kCrTop, params, x0, T0, grad0,
                                        x1, T1, grad1, xHat, listIndices, listxk,
                                        listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     skM1_projected = project_skGivenlamk1(skM1, lamk, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
     lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
     skM1_free = project_box(skM1)
     skM1, lamk = projections_skCr_lamk1(skM1_projected, lamk_projected, skM1_free, lamk,
                                         k, kCrTop, params, x0, T0, grad0,
                                         x1, T1, grad1, xHat, listIndices, listxk,
                                         listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
     # Update
     params[k-1] = mukM1
     paramsCrTop[kCrTop] = rkM1
     paramsCrTop[kCrTop + 1] = skM1
     params[k] = lamk
     if( currentCrTop < len(indCrTop) - 1):
          currentCrTop += 1
          ### UPDATE LAMK, MUK
          # Now go with the block [lamk, muk] - we have to do this here since we have already updated rkM1 and skM1
     dlamk = partial_fObj_recCr(skM1, muk, lamk, xkM1, BkM1Bk_0, xk, BkM1Bk_1, x0, B0k, xk, Bk, etakPrev, etak)
     gradParams[k] = dlamk
     # We need to know if the next receiver is on h0k1 or on hkk1
     if( nTop + 1 == indCrTop[currentCrTop] or nTop + 1 == indStTop[currentStTop] ):
          # The next receiver is on hkhk1 NOTE THAT THE PROCESS NEVER ENTERS HERE IF j = n (number of regions)
          BkBk1_0 = listBkBk1[k+1]
          BkBk1_1 = listBkBk1[k+2]
          if( nTop + 1 == indCrTop[currentCrTop]):
               rk = paramsCrTop[2*currentCrTop]
          else:
               rk = paramsStTop[2*currentStTop]
          dmuk = partial_fObj_shCr(muk, lamk, rk, x0, B0k, xk, Bk, xk, BkBk1_0, xk1, BkBk1_1, etak, listIndices[n+j+1])
          gradParams[k+1] = dmuk
          r = close_to_identity(lamk, muk)
          # Decide which type of backtracking we have to do: close or far
          if ( r<= gamma):
               # Meaning we have to do a close update
               dmuk_collapsed = partial_fObj_collapsedShooter(skM1, muk, rk, xkM1, BkM1Bk_0, xk, BkM1Bk_1,
                                                              x0, B0k, xk, Bk,
                                                              xk, BkBk1_0, xk1, BkBk1_1,
                                                              etakPrev, etak) 
               lamk, muk = backTrClose_block0k(1, k, dlamk, dmuk, dmuk_collapsed,
                                               params, x0, T0, grad0, x1, T1, grad1, xHat,
                                               listIndices, listxk, listB0k, listBk, listBkBk1,
                                               indCrTop, paramsCrTop, indStTop, paramsStTop)
          else:
               # We have a far update
               lamk, muk = backTr_block0k(1, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Now we project, we use rk to project back muk and skM1 to project back lamk
          lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
          skM1_projected = project_skGivenlamk1(skM1, lamk, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
          lamk_free = project_box(lamk)
          skM1, lamk = projections_skCr_lamk1(skM1_projected, lamk_projected, skM1, lamk_free,
                                              k, kCrTop, params, x0, T0, grad0,
                                              x1, T1, grad1, xHat, listIndices, listxk,
                                              listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          muk_projected = project_mukGivenrk(muk, rk, x0, B0k, xk, Bk, BkBk1_0, xk1, BkBk1_1)
          rk_projected = project_rkGivenmuk(rk, muk, x0, B0k, xk, Bk, xk1, Bk1, BkBk1_0, BkBk1_1)
          muk_free = project_box(muk)
          muk, rk = projections_muk_rkCr(muk_projected, rk_projected, muk_free, rk_free,
                                         k, 2*currentCrTop,  params, x0,
                                         T0, grad0, x1, T1, grad1, xHat, listIndices,
                                         listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Update
          params[k] = lamk
          params[k+1] = muk
          paramsCrTop[kCrTop + 1] = skM1
          params[2*currentCrTop] = rk
     elif(j < n):
          # The next receiver is on h0hk1
          lamk1 = params[k+2]
          B0k1 = listB0k[j+1]
          Bk1 = listBk[j+1]
          etak1 = listIndices[j+1]
          dmuk = partial_fObj_shCr(muk, lamk, lamk1, x0, B0k, xk, Bk, x0, B0k1, xk1, Bk1, etak, etak1)
          gradParams[k+1] = dmuk
          r = close_to_identity(lamk, muk)
          # Decide which type of backtracking we have to do: close or far
          if (r <= gamma):
               # Meaning we have to do a close update
               dmuk_collapsed = partial_fObj_collapsedShooter(skM1, muk, lamk1, xkM1, BkM1Bk_0, xk, BkM1Bk_1,
                                                              x0, B0k, xk, Bk,
                                                              x0, B0k1, xk1, Bk1, etak, etak1) 
               lamk, muk = backTrClose_block0k(1, k, dlamk, dmuk, dmuk_collapsed,
                                               params, x0, T0, grad0, x1, T1, grad1, xHat,
                                               listIndices, listxk, listB0k, listBk, listBkBk1,
                                               indCrTop, paramsCrTop, indStTop, paramsStTop)
               gamma = gamma*theta_gamma
               gammas[j-1] = gamma
          else:
               # Meaning we have to do a far update
               lamk, muk = backTr_block0k(1, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Now we project back, we use lamk1 to project back muk and skM1 to project back lamk
          lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
          lamk_free = project_box(lamk)
          skM1_projected = project_skGivenlamk1(skM1, lamk_free, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
          skM1, lamk = projections_skCr_lamk1(skM1_projected, lamk_projected, skM1, lamk_free,
                                              k, kCrTop, params, x0, T0, grad0,
                                              x1, T1, grad1, xHat, listIndices, listxk,
                                              listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          paramsCrTop[kCrTop + 1] = skM1
          # Depending on dmuk and the curvature of hkhk1 we project back muk differently
          muk_free = project_box(muk)
          if( dmuk >= 0 or listCurvingInwards[j] != 1):
               # Means that we don't have a point on the side edge and we don't have to worry about that edge
               muk_projected = project_mukGivenlamk1(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
               lamk1_projected = project_lamk1Givenmuk(muk_free, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
          else:
               # We have to worry about the side edge
               muk_projected = project_mukGivenlamk1_noCr(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
               lamk1_projected = project_lamkGivenmuk1_noCr(muk_free, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
          muk, lamk1 = projections_muk_lamk1(muk_projected, lamk_projected, muk_free, lamk1,
                                             k+1, params, x0, T0, grad0, x1, T1,
                                             grad1, xHat, listIndices, listxk, listB0k, listBk,
                                             listBkBk1, indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Update
          params[k] = lamk
          params[k+1] = muk
          params[k+2] = lamk1
     else:
          # We just have to update lamk because we are on lamn1
          alpha = backTr_coord(1, k, dlamk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop)
          lamk = lamk - alpha*dlamk
          lamk_free = project_box(lamk)
          lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
          skM1_projected = project_skGivenlamk1(skM1, lamk, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
          skM1, lamk = projections_skCr_lamk1(skM1_projected, lamk_projected, skM1, lamk_free,
                                              k, kCrTop, params, x0, T0, grad0,
                                              x1, T1, grad1, xHat, listIndices, listxk,
                                              listB0k, listBk, listBkBk1, indCrTop,
                                              paramsCrTop, indStTop, paramsStTop)
          params[k] = lamk
          paramsCrTop[kCrTop +1] = skM1
     #Now we just return whatever we did
     return currentCrTop, params, paramsCrTop, gradParams, gradCrTop


def updateFromStTop(n, j, currentCrTop, currentStTop, params, gammas, theta_gamma, x0, T0,
                    grad0, x1, T1, grad1, xHat, listIndices, listxk,
                    listB0k, listBk, listBkBk1, indCrTop, paramsCrTop,
                    indStTop, paramsStTop, listCurvingInwards, gradParams, gradStTop):
     '''
     Update [rkM1, skM1] and [lamk, muk] where [rkM1, skM1] define
     a straight ray which peeks into the outter side region
     from the top of the triangle fan (or side edge,
     however you want to look at it)
     '''
     k = 2*j - 1
     nTop = j
     gamma = gammas[j-1]
     mukM1 = params[k-1]
     lamk = params[k]
     muk = params[k+1]
     B0kM1 = listB0k[j-1]
     B0k = listB0k[j]
     xkM1 = listxk[j]
     xk = listxk[j+1]
     if( j < n ):
          xk1 = listxk[j+2] # Because otherwise this is not defined
          B0k1 = listB0k[j+1]
          Bk1 = listBk[j+1]
     BkM1 = listBk[j-1]
     Bk = listBk[j]
     BkM1Bk_0 = listBkBk1[k-1] # grad of hkhk1 at xk
     BkM1Bk_1 = listBkBk1[k] # grad of hkhk1 at xk1
     etakPrev = listIndices[j-1] # index of refraction from previous triangle
     etak = listIndices[j] # index of refraction in current triangle
     etaRegionOutside = listIndices[n+j] # index of refraction ouside from the side of the previous triangle
     # The update starts here
     etaMinSt = min(etaRegionOutside, etakPrev)
     kStTop = 2*currentStTop
     rkM1 = paramsStTop[kStTop]
     skM1 = paramsStTop[kStTop + 1]
     akM1 = itt.hermite_boundary(rkM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1) # receiver from mukM1
     bkM1 = itt.hermite_boundary(skM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1) # shooter to lamk
     # Compute directions
     drkM1 = partial_fObj_recSt(mukM1, skM1, rkM1, x0, B0kM1, xkM1, BkM1, xkM1, BkM1Bk_0, xk, BkM1Bk_1, etakPrev, etaRegionOutside)
     gradStTop[kStTop] = drkM1
     dskM1 = partial_fObj_shSt(skM1, rkM1, lamk, xkM1, BkM1Bk_0, xk, BkM1Bk_1, x0, B0k, xk, Bk, etakPrev, etaRegionOutside)
     gradStTop[kStTop + 1] = dskM1
     # Decide if we need close or far backtracking
     r = close_to_identity(rkM1, skM1)
     if( r <= gamma ):
          # Meaning we have to do a close update
          drkM1_collapsed = partial_fObj_collapsedShooter(mukM1, rkM1, lamk, x0, B0kM1, xkM1, BkM1,
                                                          xkM1, BkM1Bk_0, xk, BkM1Bk_1,
                                                          x0, B0k, xk, Bk, etakPrev, etak)
          rkM1, skM1 = backTrClose_blockStTop(1, kStTop, drkM1, dskM1, drkM1_collapsed,
                                              params, x0, T0, grad0, x1, T1, grad1, xHat,
                                              listIndices, listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
     else:
          # Meaning we dont have to do a close update
          rkM1, skM1 = backTr_blockStTop(1, kStTop, drkM1, dskM1, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                         listIndices, listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
     # Project back
     rkM1_projected = project_rkGivenmuk(rkM1, mukM1, x0, B0kM1, xkM1, BkM1, xk, Bk, BkM1Bk_0, BkM1Bk_1)
     mukM1_projected = project_mukGivenrk(mukM1, rkM1, x0, B0kM1, xkM1, BkM1, BkM1Bk_0, xk, BkM1Bk_1)
     rkM1_free = project_box(rkM1)
     mukM1_free = mukM1
     mukM1, rkM1 = projections_muk_rkSt(mukM1_projected, rkM1_projected, mukM1_free, rkM1_free, k-1,
                                        kStTop, params, x0, T0, grad0,
                                        x1, T1, grad1, xHat, listIndices, listxk,
                                        listB0k, listBk, listBkBk1,
                                        indCrTop, paramsCrTop, indStTop, paramsStTop)
     skM1_projected = project_skGivenlamk1(skM1, lamk, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
     lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
     skM1_free = project_box(skM1)
     skM1, lamk = projections_skSt_lamk1(skM1_projected, lamk_projected, skM1_free, lamk,
                                         k, kStTop, params, x0, T0, grad0,
                                         x1, T1, grad1, xHat, listIndices, listxk,
                                         listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
     # Update
     params[k-1] = mukM1
     paramsStTop[kStTop] = rkM1
     paramsStTop[kStTop + 1] = skM1
     params[k] = lamk
     if( currentStTop < len(indStTop) - 1):
          currentStTop += 1
     # Now we can update lamk muk - we have to do this here since we have already updated rkM1 and skM1
     dlamk = partial_fObj_recCr(skM1, muk, lamk, xkM1, BkM1Bk_0, xk, BkM1Bk_1, x0, B0k, xk, Bk, etakPrev, etak)
     gradParams[k] = dlamk
     # We need to know if the next receiver is on h0k1 or on hkk1
     if( nTop + 1 == indStTop[currentStTop] or nTop + 1 == indCrTop[currentCrTop] ):
          # The next receiver is on hkhk1 NOTE THAT THE PROCESS NEVER ENTERS HERE IF j = n (number of regions)
          BkBk1_0 = listBkBk1[k+1]
          BkBk1_1 = listBkBk1[k+2]
          if( nTop + 1 == indStTop[currentStTop]):
               rk = paramsStTop[2*currentStTop]
          else:
               rk = paramsCrTop[2*currentCrTop]
          dmuk = partial_fObj_shCr(muk, lamk, rk, x0, B0k, xk, Bk, xk, BkBk1_0, xk1, BkBk1_1, etak, listIndices[n+j+1])
          gradParams[k+1] = dmuk
          r = close_to_identity(lamk, muk)
          # Decide which type of backtracking we have to do: close or far
          if ( r<= gamma):
               # Meaning we have to do a close update
               dmuk_collapsed = partial_fObj_collapsedShooter(skM1, muk, rk, xkM1, BkM1Bk_0, xk, BkM1Bk_1,
                                                              x0, B0k, xk, Bk,
                                                              xk, BkBk1_0, xk1, BkBk1_1,
                                                              etakPrev, etak) 
               lamk, muk = backTrClose_block0k(1, k, dlamk, dmuk, dmuk_collapsed,
                                               params, x0, T0, grad0, x1, T1, grad1, xHat,
                                               listIndices, listxk, listB0k, listBk, listBkBk1,
                                               indCrTop, paramsCrTop, indStTop, paramsStTop)
          else:
               # We have a far update
               lamk, muk = backTr_block0k(1, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Now we project, we use rk to project back muk and skM1 to project back lamk
          lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
          skM1_projected = project_skGivenlamk1(skM1, lamk, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
          lamk_free = project_box(lamk)
          skM1, lamk = projections_skSt_lamk1(skM1_projected, lamk_projected, skM1, lamk_free,
                                              k, kStTop, params, x0, T0, grad0,
                                              x1, T1, grad1, xHat, listIndices, listxk,
                                              listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          muk_projected = project_mukGivenrk(muk, rk, x0, B0k, xk, Bk, BkBk1_0, xk1, BkBk1_1)
          rk_projected = project_rkGivenmuk(rk, muk, x0, B0k, xk, Bk, xk1, Bk1, BkBk1_0, BkBk1_1)
          muk_free = project_box(muk)
          muk, rk = projections_muk_rkSt(muk_projected, rk_projected, muk_free, rk_free,
                                         k, 2*currentStTop,  params, x0,
                                         T0, grad0, x1, T1, grad1, xHat, listIndices,
                                         listxk, listB0k, listBk, listBkBk1,
                                         indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Update
          params[k] = lamk
          params[k+1] = muk
          paramsStTop[kStTop + 1] = skM1
          params[2*currentStTop] = rk
     elif( j < n):
          # The next receiver is on h0hk1
          lamk1 = params[k+2]
          B0k1 = listB0k[j+1]
          Bk1 = listBk[j+1]
          etak1 = listIndices[j+1]
          dmuk = partial_fObj_shCr(muk, lamk, lamk1, x0, B0k, xk, Bk, x0, B0k1, xk1, Bk1, etak, etak1)
          gradParams[k+1] = dmuk
          r = close_to_identity(lamk, muk)
          # Decide which type of backtracking we have to do: close or far
          if (r <= gamma):
               # Meaning we have to do a close update
               dmuk_collapsed = partial_fObj_collapsedShooter(skM1, muk, lamk1, xkM1, BkM1Bk_0, xk, BkM1Bk_1,
                                                              x0, B0k, xk, Bk,
                                                              x0, B0k1, xk1, Bk1, etak, etak1) 
               lamk, muk = backTrClose_block0k(1, k, dlamk, dmuk, dmuk_collapsed,
                                               x0, T0, grad0, x1, T1, grad1, xHat,
                                               listIndices, listxk, listB0k, listBk, listBkBk1,
                                               indCrTop, paramsCrTop, indStTop, paramsStTop)
               gamma = gamma*theta_gamma
               gammas[j-1] = gamma
          else:
               # Meaning we have to do a far update
               lamk, muk = backTr_block0k(1, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Now we project back, we use lamk1 to project back muk and skM1 to project back lamk
          lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
          lamk_free = project_box(lamk)
          skM1_projected = project_skGivenlamk1(skM1, lamk_free, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
          skM1, lamk = projections_skSt_lamk1(skM1_projected, lamk_projected, skM1, lamk_free,
                                              k, kStTop, params, x0, T0, grad0,
                                              x1, T1, grad1, xHat, listIndices, listxk,
                                              listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          paramsStTop[kStTop + 1] = skM1          # Depending on dmuk and the curvature of hkhk1 we project back muk differently
          muk_free = project_box(muk)
          if( dmuk >= 0 or listCurvingInwards[j] != 1):
               # Means that we don't have a point on the side edge and we don't have to worry about that edge
               muk_projected = project_mukGivenlamk1(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
               lamk1_projected = project_lamk1Givenmuk(muk_free, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
          else:
               # We have to worry about the side edge
               muk_projected = project_mukGivenlamk1_noCr(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
               lamk1_projected = project_lamkGivenmuk1_noCr(muk_free, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
          muk, lamk1 = projections_muk_lamk1(muk_projected, lamk_projected, muk_free, lamk1,
                                             k+1, params, x0, T0, grad0, x1, T1,
                                             grad1, xHat, listIndices, listxk, listB0k, listBk,
                                             listBkBk1, indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Update
          params[k] = lamk
          params[k+1] = muk
          params[k+2] = lamk1
     else:
          # We just have to update lamk because we are on lamn1
          alpha = backTr_coord(1, k, dlamk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop)
          lamk = lamk - alpha*dlamk
          lamk_free = project_box(lamk)
          lamk_projected = project_lamkGivenskM1(lamk, skM1, x0, B0k, xk, Bk, BkM1Bk_0, xkM1, BkM1Bk_1)
          skM1_projected = project_skGivenlamk1(skM1, lamk, x0, B0k, xk, Bk, xkM1, BkM1Bk_0, BkM1Bk_1)
          skM1, lamk = projections_skSt_lamk1(skM1_projected, lamk_projected, skM1, lamk_free,
                                              k, kStTop, params, x0, T0, grad0,
                                              x1, T1, grad1, xHat, listIndices, listxk,
                                              listB0k, listBk, listBkBk1, indCrTop,
                                              paramsCrTop, indStTop, paramsStTop)
          params[k] = lamk
          paramsStTop[kStTop +1] = skM1
     #Now we just return whatever we did
     return currentStTop, params, paramsStTop, gradParams, gradStTop


def udapteFromh0kM1(n, j, currentCrTop, currentStTop, params, gammas, theta_gamma, x0, T0,
                    grad0, x1, T1, grad1, xHat, listIndices, listxk,
                    listB0k, listBk, listBkBk1, indCrTop, paramsCrTop,
                    indStTop, paramsStTop, listCurvingInwards, gradParams):
     '''
     Update just [lamk, muk]. In this case the previous shooter comes from h0kM1
     (i.e. the previous shooter from lamk is mukM1)
     '''
     k = 2*j - 1
     nTop = j
     gamma = gammas[j-1]
     mukM1 = params[k-1]
     lamk = params[k]
     muk = params[k+1]
     B0kM1 = listB0k[j-1]
     BkM1Bk_0 = listBkBk1[k-1] # grad of hkhk1 at xk
     BkM1Bk_1 = listBkBk1[k] # grad of hkhk1 at xk1
     B0k = listB0k[j]
     xkM1 = listxk[j]
     xk = listxk[j+1]
     if( j < n ):
          xk1 = listxk[j+2] # Because otherwise this is not defined
          B0k1 = listB0k[j+1]
          Bk = listBk[j+1]
     BkM1 = listBk[j-1]
     Bk = listBk[j]
     etakPrev = listIndices[j-1] # index of refraction from previous triangle
     etak = listIndices[j] # index of refraction in current triangle
     etaRegionOutside = listIndices[n+j] # index of refraction ouside from the side of the previous triangle
     # The update starts here
     dlamk = partial_fObj_recCr(mukM1, muk, lamk, x0, B0kM1, xkM1, BkM1, x0, B0k, xk, Bk, etakPrev, etak)
     #breakpoint()
     gradParams[k] = dlamk
     # We need to see if the next receiver is on h0k1 or on hkk1
     if( nTop + 1 == indStTop[currentStTop] or nTop + 1 == indCrTop[currentStTop] ):
          # This means that the next receiver is on hkk1 NOTE THAT THE PROCESS NEVER GOES HERE IF j = n, if we are on lamn1
          BkBk1_0 = listBkBk1[k+1]
          BkBk1_1 = listBkBk1[k+2]
          if( nTop + 1 == indStTop[currentStTop]):
               rk = paramsStTop[2*currentStTop]
          else:
               rk = paramsCrTop[2*currentCrTop]
          xk1 = listxk[j+2]
          dmuk = partial_fObj_shCr(muk, lamk, rk, x0, B0k, xk, Bk, xk, BkBk1_0, xk1, BkBk1_1, etak, listIndices[n+j+1])
          gradParams[k+1] = dmuk
          # Decide which type of backtracking we have to do: close or far
          r = close_to_identity(lamk, muk)
          if( r <= gamma ):
               # Meaning we have a close update
               dmuk_collapsed = partial_fObj_collapsedShooter(mukM1, muk, rk, x0, BkM1, xkM1, BkM1,
                                                              x0, B0k, xk, Bk,
                                                              xk, BkBk1_0, xk1, BkBk1_1, etakPrev, etak) 
               lamk, muk = backTrClose_block0k(1, k, dlamk, dmuk, dmuk_collapsed,
                                               x0, T0, grad0, x1, T1, grad1, xHat,
                                               listIndices, listxk, listB0k, listBk, listBkBk1,
                                               indCrTop, paramsCrTop, indStTop, paramsStTop)
               gamma = gamma*theta_gamma
               gammas[j-1] = gamma
          else:
               # Meaning we have to do a far update
               lamk, muk = backTr_block0k(1, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Now we project back, we use rk to project back muk and mukM1 to project back lamk
          # See if we need to use hkM1k to project back lamk
          if( dlamk > 0 or listCurvingInwards[j-1] != 1):
               # There is no problem with hkM1k
               lamk_projected = project_lamk1Givenmuk(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
               mukM1_projected = project_mukGivenlamk1(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
          else:
               # There is a problem with hkM1k
               lamk_projected = project_lamkGivenmuk1_noCr(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, BkM1Bk_0, BkM1Bk_1)
               mukM1_projected = project_mukGivenlamk1_noCr(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, BkM1Bk_0, BkM1Bk_1)
          lamk_free = project_box(lamk)
          mukM1, lamk = projections_muk_lamk1(mukM1_projected, lamk_projected, mukM1, lamk_free, k-1,
                                              params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                              listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          # Project back muk using rk
          muk_projected = project_mukGivenrk(muk, rk, x0, B0k, xk, Bk, BkBk1_0, xk1, BkBk1_1)
          muk_free = project_box(muk)
          rk_projected = project_rkGivenmuk(rk, muk, x0, B0k, xk, Bk, xk1, Bk1, BkBk1_0, BkBk1_1)
          if( nTop + 1 == indStTop[currentStTop]):
               muk, rk = projections_muk_rkSt(muk_projected, rk_projected, muk_free, rk, k-1, 2*currentStTop,
                                              params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                              listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
               paramsStTop[2*currentStTop] = rk
          else:
               muk, rk = projections_muk_rkCr(muk_projected, rk_projected, muk_free, rk, k-1, 2*currentCrTop,
                                              params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                              listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
               paramsCrTop[2*currentCrTop] = rk
          # Update
          params[k-1] = mukM1
          params[k] = lamk
          params[k+1] = muk
     elif(j <  n):
          # This means that the next receiver is on h0k1
          breakpoint()
          dmuk = partial_fObj_shCr(muk, lamk, lamk1, x0, B0k, xk, Bk, x0, B0k1, xk1, Bk1, etak, etakM1)
          gradParams[k+1] = dmuk
          # Decide which type of backtracking we have to do: close or far
          r = close_to_identity(lamk, muk)
          if( r<= gamma):
               # Meaning we have a close update
               dmuk_collapsed = partial_fObj_collapsedShooter(mukM1, muk, lamk1, x0, BkM1, xkM1, BkM1,
                                                              x0, B0k, xk, Bk,
                                                              x0, B0k1, xk1, Bk1, etakPrev, etak) 
               lamk, muk = backTrClose_block0k(1, k, dlamk, dmuk, dmuk_collapsed,
                                               params, x0, T0, grad0, x1, T1, grad1, xHat,
                                               listIndices, listxk, listB0k, listBk, listBkBk1,
                                               indCrTop, paramsCrTop, indStTop, paramsStTop)
               gamma = gamma*theta_gamma
               gammas[j-1] = gamma
          else:
               # Meaning we have to do a far update
               lamk, muk = backTr_block0k(1, k, dlamk, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                                          listIndices, listxk, listB0k, listBk, listBkBk1,
                                          indCrTop, paramsCrTop, indStTop, paramsStTop)
               # Now we project back, we use rk to project back muk and mukM1 to project back lamk
               # See if we need to use hkM1k to project back lamk
               if( dlamk > 0 or listCurvingInwards[j-1] != 1):
                    # There is no problem with hkM1k
                    lamk_projected = project_lamk1Givenmuk(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
                    mukM1_projected = project_mukGivenlamk1(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
               else:
                    # There is a problem with hkM1k
                    lamk_projected = project_lamkGivenmuk1_noCr(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, BkM1Bk_0, BkM1Bk_1)
                    mukM1_projected = project_mukGivenlamk1_noCr(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, BkM1Bk_0, BkM1Bk_1)
               lamk_free = project_box(lamk)
               mukM1, lamk = projections_muk_lamk1(mukM1_projected, lamk_projected, mukM1, lamk_free, k-1,
                                                   params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                                   listxk, listB0k, listBk, listBkBk1,
                                                   indCrTop, paramsCrTop, indStTop, paramsStTop)
               # Project back muk using lamk1
               if( dmuk > 0 or listCurvingInwards[j] != 1):
                    # There is no problem with hkk1
                    muk_projected = project_mukGivenlamk1(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
                    lamk1_projected = project_lamk1Givenmuk(muk, lamk1, x0, B0k, xk, B, B0k1, xk1, Bk1)
               else:
                    # There is a problem with hkk1
                    BkBk1_0 = listBkBk1[k+1]
                    BkBk1_1 = listBkBk1[k+2]
                    muk_projected = project_mukGivenlamk1_noCr(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
                    lamk1_projected = project_lamkGivenmuk1_noCr(muk, lamk1, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
               muk_free = project_box(muk)
               muk, lamk1 = projections_muk_lamk1(muk_projected, lamk1_projected, muk_free, lamk1, k+1,
                                                  params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                                  listxk, listB0k, listBk, listBkBk1,
                                                  indCrTop, paramsCrTop, indStTop, paramsStTop)
               params[k] = lamk
               params[k+1] = muk
               params[k+2] = lamk1
     else:
          # This means that j = n, we are on lamn1
          alpha = backTr_coord(1, k, dlamk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                               listIndices, listxk, listB0k, listBk, listBkBk1,
                               indCrTop, paramsCrTop, indStTop, paramsStTop)
          lamk = lamk - alpha*dlamk
          lamk_free = project_box(lamk)
          if( dlamk > 0 or listCurvingInwards[j-1] != 1):
               # There is no problem with hkM1k
               lamk_projected = project_lamk1Givenmuk(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
               mukM1_projected = project_mukGivenlamk1(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk)
          else:
               # There is a problem with hkM1k
               lamk_projected = project_lamkGivenmuk1_noCr(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, BkM1Bk_0, BkM1Bk_1)
               mukM1_projected = project_mukGivenlamk1_noCr(mukM1, lamk, x0, B0kM1, xkM1, BkM1, B0k, xk, Bk, BkM1Bk_0, BkM1Bk_1)
          mukM1, lamk = projections_muk_lamk1(mukM1_projected, lamk_projected, mukM1, lamk_free, k-1,
                                              params, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                              listxk, listB0k, listBk, listBkBk1,
                                              indCrTop, paramsCrTop, indStTop, paramsStTop)
          params[k-1] = mukM1
          params[k] = lamk
     return params, paramsCrTop, paramsStTop, gradParams
     


def forwardPassUpdate(params0, gammas, theta_gamma, x0, T0, grad0, x1, T1, grad1, xHat,
                      listIndices, listxk, listB0k, listBk, listBkBk1,
                      indCrTop, paramsCrTop0, indStTop, paramsStTop0, listCurvingInwards):
     '''
     gammas: radius for the circle centered at [lamk, muk], if such circle intersects the line lamk = muk
     then do a close update.
     Updates blocks     [ mu1 ]
                        [ lamk, muk ]
                        [  rk,   sk ]
                        [ lamn1 ]
     listCurvingInwards is just a list of length n such that if the k-th side edge of the triangle
     fan is curving inwards (to the triangle fan) then listCurvingInwards[k] = 1, it's 0 if it's not
     In the foor loop the strategy is to update blocks in the following order
                        [  rkM1,  skM1 ]   (if applicable)
                        [  lamk,  muk  ]     
     '''
     gradParams = np.empty((len(params0), 1))
     gradParams[-1] = 0 # Because mun1 = 1 always
     # Set indStTop if indCrTop is given
     if(paramsStTop0 is None or indStTop is None):
          indStTop = [-1]
          paramsStTop = np.array([0,0])
          gradStTop = np.empty((0,1))
     else:
          paramsStTop = np.copy(paramsStTop0)
          gradStTop = np.empty((len(paramsStTop0), 1))
     # Set indCrTop if indStTop is given
     if(paramsCrTop0 is None or indCrTop is None):
          indCrTop = [-1]
          paramsCrTop = np.array([0,0])
          gradCrTop = np.empty((0, 1))
     else:
          paramsCrTop = np.copy(paramsCrTop0)
          gradCrTop = np.empty((len(paramsCrTop), 1))
     currentCrTop = 0
     currentStTop = 0
     # First parameter to update: mu1
     params = np.copy(params0)
     mu1 = params[0]
     lam2 = params[1]
     B0k = listB0k[0]
     etak = listIndices[0]
     xk = listxk[1]
     Bk = listBk[0]
     B0k1 = listB0k[1]
     xk1 = listxk[2]
     Bk1 = listBk[1]
     B0k_muk = itt.gradientBoundary(mu1, x0, B0k, xk, Bk)
     BkBk1_0 = listBkBk1[0]
     BkBk1_1 = listBkBk1[1]
     yk1 = itt.hermite_boundary(lam2, x0, B0k1, xk1, Bk1)
     zk = itt.hermite_boundary(mu1, x0, B0k, xk, Bk)
     # Compute direction for muk
     if( indCrTop[0] != 1 and indStTop[0] != 1):
          dmuk = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B0k_muk, yk1, zk)
     else:
          rk = paramsCrTop[0]
          ak = itt.hermite_boundary( rk, xk, BkBk1_0, xk1, BkBk1_1)
          #breakpoint()      ##############################################################################
          dmuk = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B0k_muk, ak, zk)
     alpha = backTr_coord(2, 0, dmuk, params, x0, T0, grad0, x1, T1, grad1, xHat,
                          listIndices, listxk, listB0k, listBk, listBkBk1,
                          indCrTop, paramsCrTop, indStTop, paramsStTop)
     mu1 = mu1 - alpha*dmuk
     #breakpoint() ##############################################################################
     # Then we need to project it back WE ASSUME THAT PARAMS0 IS A FEASIBLE SET!
     if( indCrTop[0] != 1 and indStTop[0] != 1 and listCurvingInwards[0] != 1):
          # Means that we don't have a point on the side edge and we don't have to worry about that edge
          # Consider moving mu1, consider moving lam2
          mu1_projected = project_mukGivenlamk1(mu1, lam2, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
          mu1_free = project_box(mu1)
          lam2_projected = project_lamk1Givenmuk(mu1_free, lam2, x0, B0k, xk, Bk, B0k1, xk1, Bk1)
          mu1, lam2 = projections_muk_lamk1(mu1_projected, lam2_projected, mu1_free, lam2, 0, params, x0,
                                            T0, grad0, x1, T1, grad1, xHat, listIndices,
                                            listxk, listB0k, listBk, listBkBk1,
                                            indCrTop, paramsCrTop, indStTop, paramsStTop)
          params[0] = mu1
          params[1] = lam2
          gradParams[0] = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B0k_muk, yk1, zk)
     elif( indCrTop[0] != 1 and indStTop[0] != 1 ):
          # Means that we don't have a point on the side edge but we do need to worry about it
          mu1_projected = project_mukGivenlamk1_noCr(mu1, lam2, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
          mu1_free = project_box(mu1)
          lam2_projected = project_lamkGivenmuk1_noCr(mu1_free, lam2, x0, B0k, xk, Bk, B0k1, xk1, Bk1, BkBk1_0, BkBk1_1)
          mu1, lam2 = projections_muk_lamk1(mu1_projected, lam2_projected, mu1_free, lam2, 0, params, x0,
                                            T0, grad0, x1, T1, grad1, xHat, listIndices,
                                            listxk, listB0k, listBk, listBkBk1,
                                            indCrTop, paramsCrTop, indStTop, paramsStTop)
          params[0] = mu1
          params[1] = lam2
          gradParams[0] = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B0k_muk, yk1, zk)
     else:
          # Means that we have a point on the side edge, we can project using that point
          mu1_projected = project_mukGivenrk(mu1, rk, x0, B0k, xk, Bk, BkBk1_0, xk1, BkBk1_1)
          mu1_free = project_box(mu1)
          if( indCrTop[0] == 1):
               r1 = paramsCrTop[0]
               r1_projected = project_rkGivenmuk(r1, mu1_free, x0, B0k, xk, Bk, xk1, Bk1, BkBk1_0, BkBk1_1)
               mu1, r1 = projections_muk_rkCr(mu1_projected, r1_projected, mu1_free, r1, 0, 0, params, x0,
                                            T0, grad0, x1, T1, grad1, xHat, listIndices,
                                            listxk, listB0k, listBk, listBkBk1,
                                            indCrTop, paramsCrTop, indStTop, paramsStTop)
               paramsCrTop[0] = r1
          else:
               r1 = paramsStTop[0]
               r1_projected = project_rkGivenmuk(r1, mu1_free, x0, B0k, xk, Bk, xk1, Bk1, BkBk1_0, BkBk1_1)
               mu1, r1 = projections_muk_rkSt(mu1_projected, r1_projected, mu1_free, r1, 0, 0, params, x0,
                                            T0, grad0, x1, T1, grad1, xHat, listIndices,
                                            listxk, listB0k, listBk, listBkBk1,
                                            indCrTop, paramsCrTop, indStTop, paramsStTop)
               paramsStTop[0] = r1
          params[0] = mu1
          gradParams[0] = partial_fObj_mu1(mu1, x0, T0, grad0, x1, T1, grad1, B0k_muk, ak, zk)
     #breakpoint() ##############################################################################
     # Now we start with blocks of size 2
     n = len(listxk) - 2
     for j in range(1, n+1):
          # j goes from 1 to n
          #breakpoint() ##############################################################################
          ###### Compute the directions - first we want to update [  rkM1,  skM1 ] if applicable
          if( j == indCrTop[currentCrTop] ):
               # This means that the current muk comes from a point on the side edge of the previous triangle (comes from creeping)
               currentCrTop, params, paramsCrTop, gradParams, gradCrTop = updateFromCrTop(n, j, currentCrTop, currentStTop,
                                                                                          params, gammas, theta_gamma,
                                                                                          x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                          listIndices, listxk, listB0k, listBk,
                                                                                          listBkBk1, indCrTop, paramsCrTop, indStTop,
                                                                                          paramsStTop, listCurvingInwards,
                                                                                          gradParams, gradCrTop)
          elif( j == indStTop[currentStTop]):
               # This means that the current muk comes from a point on the side edge of the previous triangle (comes from a straight ray)
               currentStTop, params, paramsStTop, gradParams, gradStTop = updateFromStTop(n, j, currentCrTop, currentStTop,
                                                                                          params, gammas, theta_gamma,
                                                                                          x0, T0, grad0, x1, T1, grad1, xHat,
                                                                                          listIndices, listxk, listB0k, listBk,
                                                                                          listBkBk1, indCrTop, paramsCrTop, indStTop,
                                                                                          paramsStTop, listCurvingInwards,
                                                                                          gradParams, gradStTop)
          else:
               # Meaning that we dont have [rkM1, skM1] on the side of the edge, in this case we just need to update [lamk, muk]
                params, paramsCrTop, paramsStTop, gradParams = udapteFromh0kM1(n, j, currentCrTop, currentStTop, params, gammas,
                                                                               theta_gamma, x0, T0, grad0, x1, T1, grad1, xHat,
                                                                               listIndices, listxk, listB0k, listBk, listBkBk1,
                                                                               indCrTop, paramsCrTop, indStTop, paramsStTop,
                                                                               listCurvingInwards, gradParams)
     # THIS IS THE END OF THE FOR LOOP
     # We are done
     return params, paramsCrTop, paramsStTop, gradParams, gradCrTop, gradStTop



def blockCoordinateGradient_generalized(params0, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                        listxk, listB0k, listBk, listBkBk1, indCrTop, paramsCrTop0,
                                        indStTop, paramsStTop0, listCurvingInwards, theta_gamma = 1,
                                        tol = 1e-12, maxIter = 30, plotSteps = False):
     '''
     Block coordiante subgradient descent (modified) for a generalized triangle fan.
     '''
     paramsk = np.copy(params0)
     n = len(listxk) - 2
     # Add artificial mun1 if necessary
     if( len(paramsk) < 2*n + 1):
          paramsk = np.append(paramsk, [1])
     if( paramsCrTop0 == None):
          paramsCrTopk = np.array([0,0])
     else:
          paramsCrTopk = np.copy(paramsCrTop0)
     if( paramsStTop0 == None):
          paramsStTopk = np.array([0,0])
     else:
          paramsStTopk = np.copy(paramsStTop0)
     # Initialize useful things
     listObjVals = []
     listGradNorms = []
     listChangefObj = []
     listChangeParams = []
     gammas = 0.1*np.ones((n + 1))
     iter = 0
     change_fVal = 1
     normChangeP = 1
     while( abs(change_fVal) > tol and iter < maxIter and normChangeP > tol):
          paramskM1, paramsCrTopkM1, paramsStTopkM1 = paramsk, paramsCrTopk, paramsStTopk
          # Forward pass
          paramsk, paramsCrTopk, paramsStTopk, gradParamsk, gradCrTopk, gradStTopk = forwardPassUpdate(paramsk, gammas,
                                                                                                 theta_gamma, x0,
                                                                                                 T0, grad0, x1,
                                                                                                 T1, grad1, xHat,
                                                                                                 listIndices, listxk,
                                                                                                 listB0k, listBk, listBkBk1,
                                                                                                 indCrTop, paramsCrTopk,
                                                                                                 indStTop, paramsStTopk,
                                                                                                 listCurvingInwards)
          # Update lists
          fk = fObj_generalized(paramsk, x0, T0, grad0, x1, T1, grad1, xHat, listIndices,
                                listxk, listB0k, listBk, listBkBk1,
                                indCrTop = indCrTop, paramsCrTop = paramsCrTopk,
                                indStTop = indStTop, paramsStTop = paramsStTopk)
          gradk = sqrt( norm(gradParamsk)**2 + norm(gradCrTopk)**2 + norm(gradStTopk)**2)
          normChangeP = sqrt( norm( paramskM1 - paramsk)**2 + norm( paramsCrTopkM1 - paramsCrTopk)**2 + norm(paramsStTopkM1 - paramsStTopk)**2 )
          listChangeParams.append(normChangeP)
          listObjVals.append(fk)
          listGradNorms.append(gradk)
          if( iter > 0):
               change_fVal = listObjVals[-2] - fk
               listChangefObj.append(change_fVal)
          if( plotSteps):
               itt.plotFann(x0, listB0k, listxk, listBk, params = paramsk, listBkBk1 = listBkBk1, indStTop = indStTop, paramsStTop = paramsStTopk, indCrTop = indCrTop, paramsCrTop = paramsCrTopk)
               plt.title("Triangle fan after " + str(iter) + " iterations")
          iter += 1
     return paramsk, paramsCrTopk, paramsStTopk, gradParamsk, gradCrTopk, gradStTopk, listObjVals, listGradNorms, listChangefObj, listChangeParams
          
          


