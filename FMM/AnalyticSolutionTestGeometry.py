from cmath import inf
from scipy.optimize import NonlinearConstraint, minimize
import numpy as np
from numpy.linalg import norm 
from math import sqrt, asin, pi
import matplotlib.pyplot as plt
from numpy import subtract as sb

colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)
my_dpi=96
mIter = 100

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

def regA1Bool(xhat):
    '''
    This function determines if xhat can be accessed directly from x0
    '''
    ytan = (26*(7 + sqrt(14))*xhat[0])/(273 - 42*sqrt(13)) + (140*sqrt(13) + 130*sqrt(14))/(91 - 14*sqrt(13))
    if ( (xhat[1] >= ytan) | (xhat[1]<= -10) |  ((xhat[0]<=0) & (norm(xhat)>= 10) & ( xhat[1]<= 20/sqrt(13))  ) ):
        ans = True
    else:
        ans = False
    return ans

def segmentBetweenPoints(lam, a, b):
    '''
    Parametrization of the segment from a to b (using lam)
    '''
    return (1-lam)*a + lam*b

def inCtBool(x):
    '''
    Function that determines if x is on Ct (top of the snowman)
    '''
    if ( sqrt( x[0]**2 + (x[1]-5*sqrt(2))**2  )== 5*sqrt(2) & norm(x)>=10 ):
        ans = 1
    return bool(ans == 1)

def inCbBool(x):
    '''
    Function that determines if x is on Cb (bottom of the snowman AND NOT THE ARCH BETWEEN reg2 and reg3)
    '''
    if (norm(x)== 10 & x[1]< 5*sqrt(2)):
        ans = 1
    return bool(ans == 1)

def inArchBool(x):
    '''
    Function that determines if x is on the arch that separates reg2 from reg3
    '''
    if(norm(x)==10 & x[1]> 5*sqrt(2)):
        ans = 1
    return bool(ans == 1)

def parametrization(param):
    '''
    This function is going to help us in the optimization part so that the
    segments parametrized are valid for their parameters between 0 and 1
    '''
    return -(param-1)*(param)

def segmentInReg3(param, x1, x2):
    '''
    Auxiliary function for the constrain that the segment from x1 to x2 is inside reg3
    '''
    return 10-norm( param*x2 + (1-param)*x1 )

def inTopCircle(param, x1, x2):
    '''
    Auxiliary function for the constrain that the segment from x1 to x2 is in the top circle (NOT REG2)
    '''
    return 5*sqrt(2) - norm( param*x2 + (1-param)*x1 )

def notInTopCircle(param, x1, x2):
    '''
    Auxiliary function for the constrain that the segment from x1 to x2 is NOT in the top circle
    '''
    return norm( param*x2 + (1-param)*x1 ) - 5*sqrt(2)

def notInSegment3(param, x1, x2):
    '''
    Auxiliary function for the contraint that the segment from x1 to x2 is NOT inside reg3
    '''
    return norm( param*x2 + (1-param)*x1 ) - 10

def onCb(x):
    '''
    Auxiliary function for the constraint that x is on Cb (including the arch)
    '''
    return norm(x)-10

def onCt(x):
    '''
    Auxiliary function for the constraint that x is on Ct (including the bottom half of the circle)
    '''
    return sqrt(x[0]**2 + (x[1]-5*sqrt(2) )**2)

def topCb(x):
    '''
    Auxiliary function for the constraint that x[1]>= 5*sqrt(2)
    '''
    return x[1]-5*sqrt(2)

def bottomCb(x):
    '''
    Auxiliary function for the constraint that x[1]<= 5*sqrt(2)
    '''
    return 5*sqrt(2) - x[1]

def ArcLengthCb(pb1, pb2):
    '''
    Auxiliary function to calculate the arclength from pb1 to pb2, both of them on Cb
    '''
    if ((norm(sb(pb1, pb2))/20) <= 1):
        ans = inf
    else:
        if(  (pb1[1]*pb2[1])<= 0  ):
            ans = 20*pi - 20*asin( (norm(sb(pb1, pb2))/20) )
        else:
            ans = 20*asin( (norm(sb(pb1, pb2))/20) )
    return ans

def ArcLengthCt(pt1, pt2):
    '''
    Auxiliary function to calculate the arclength from pt1 to pt2, both of them on Ct
    '''
    if ((norm(sb(pt1, pt2))/(10*sqrt(2))) <= 1):
        ans = inf
    else:
        if ((pt1[1]-5*sqrt(2))*(pt2[1]-5*sqrt(2)) <= 0 ):
            ans = 10*sqrt(2)*pi - (10*sqrt(2))*asin( (norm(sb(pt1, pt2))/(10*sqrt(2))) )
        else:
            ans = (10*sqrt(2))*asin( (norm(sb(pt1, pt2))/(10*sqrt(2))) )
    return ans
    

def SnowSolution(xhat):
    '''
    This is the analytic solution to the Eikonal equation with the snowman domain starting at [-15, -10]
    '''
    x0 = np.array(  [-15, -10]  )
    p1 = np.array( [0, -10] )
    p3 = np.array( [5*sqrt(2), 5*sqrt(2)] )
    p4 = np.array( [-5*sqrt(2), 5*sqrt(2)] )
    regHat = whichRegion(xhat) # we need to know in which region xhat is
    if(regHat == 1):
        if( regA1Bool(xhat)  ):
            eik = norm( sb(x0, xhat) )
            type_route = 0
        else:
            # TYPE 1 1->3->1
            f_t3 = lambda z: norm( sb(x0, z[0:2]) ) + 5*norm( sb(z[0:2], z[2:4]) ) + norm( sb(z[2:4], xhat) )
            cons3 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
                    {'type':'eq', 'fun': lambda z: onCb(z[2:4] ) }, #pt2 on Cb
                    {'type':'ineq', 'fun': lambda z: notInSegment3(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
                    {'type':'ineq', 'fun': lambda z: notInTopCircle(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
                    {'type':'ineq', 'fun': lambda z: parametrization(z[4])  }, # 0<= lambda <= 1
                    {'type':'ineq', 'fun': lambda z: segmentInReg3(z[5], z[0:2], z[2:4]) }, #segment from pt1 to pt3 in reg3
                    {'type':'ineq', 'fun': lambda z: parametrization(z[5])  }, # 0<= mu <= 1
                    {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
                    {'type':'ineq', 'fun': lambda z: bottomCb(z[2:4])  } #pt2 not in arch
                    ]
            t3_opt = minimize(f_t3, [0, -10, 5*sqrt(2), 5*sqrt(2), 0.5, 0.5], constraints = cons3, options={'maxiter':mIter} )
            t3 = f_t3( t3_opt.x)
            # TYPE 2 1->3->2->1
            f_t4 = lambda z: norm( sb(x0, z[0:2]) ) + 5*norm( sb(z[0:2], z[2:4]) ) + 2*norm( sb(z[2:4], z[4:6]) ) +  norm( sb(z[4:6], xhat) )
            cons4 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
                    {'type':'eq', 'fun': lambda z: onCb(z[2:4] ) }, #pt2 on Cb
                    {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
                    {'type':'ineq', 'fun': lambda z: topCb(z[2:4])  }, #pt2 in arch
                    {'type':'eq', 'fun': lambda z: onCt(z[4:6])  }, #pt3 in Ct
                    {'type':'ineq', 'fun': lambda z: notInSegment3(z[6], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
                    {'type':'ineq', 'fun': lambda z: notInTopCircle(z[6], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
                    {'type':'ineq', 'fun': lambda z: parametrization(z[6])  }, # 0<= lambda <= 1
                    {'type':'ineq', 'fun': lambda z: segmentInReg3(z[7], z[0:2], z[2:4]) }, #segment from pt1 to pt3 in reg3
                    {'type':'ineq', 'fun': lambda z: parametrization(z[7])  }, # 0<= mu <= 1
                    {'type':'ineq', 'fun': lambda z: inTopCircle(z[8], z[2:4], z[4:6])  }, #segment from pt2 to pt3 in reg2
                    {'type':'ineq', 'fun': lambda z: parametrization(z[8])  }, # 0<= alpha <= 1
                    {'type':'ineq', 'fun': lambda z: notInSegment3(z[8], z[2:4], z[4:6])  } #segment from pt2 to pt3 not in reg3
                    ]
            t4_opt = minimize(f_t4, [0, -10, 5*sqrt(2), 5*sqrt(2), 0, 10*sqrt(2), 0.5, 0.5, 0.5], constraints = cons4, options={'maxiter':mIter})
            t4 = f_t4(t4_opt.x)
            # TYPE 3 1->around 3->1
            f_t5 = lambda z: norm( sb(x0, z[0:2]) ) + ArcLengthCb(z[0:2], z[2:4]) + norm( sb(z[2:4], xhat) )
            cons5 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
                    {'type':'eq', 'fun': lambda z: onCb(z[2:4] ) }, #pt2 on Cb
                    {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
                    {'type':'ineq', 'fun': lambda z: bottomCb(z[2:4])  }, #pt2 not in arch
                    {'type':'ineq', 'fun': lambda z: notInSegment3(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
                    {'type':'ineq', 'fun': lambda z: notInTopCircle(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
                    {'type':'ineq', 'fun': lambda z: parametrization(z[4])  }, # 0<= lambda <= 1
                    {'type':'ineq', 'fun': lambda z: segmentInReg3(z[5], z[0:2], z[2:4]) }, #segment from pt1 to pt2 in reg3
                    {'type':'ineq', 'fun': lambda z: parametrization(z[5])  } # 0<= mu <= 1
                    ]
            t5_opt = minimize(f_t5, [0, -10, 5*sqrt(2),-5*sqrt(2) , 0.5, 0.5], constraints = cons5, options={'maxiter':mIter} )
            t5 = f_t5(t5_opt.x)
            # # TYPE 4  1 -> 3 -> 1 -> 2 -> 1
            # f_t6 = lambda z: norm( sb(x0, z[0:2]) ) + 5*norm( sb(z[0:2], z[2:4]) ) + norm( sb(z[2:4], z[4:6]) + 2*norm( sb(z[4:6], z[6:8]) ) + norm( sb( z[6:8], xhat )  )  )
            # cons6 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
            #         {'type':'eq', 'fun': lambda z: onCb(z[2:4] ) }, #pt2 on Cb
            #         {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
            #         {'type':'ineq', 'fun': lambda z: bottomCb(z[2:4])  }, #pt2 not in arch
            #         {'type':'eq', 'fun': lambda z: onCt(z[4:6])  }, #pt3 in Ct
            #         {'type':'eq', 'fun': lambda z: onCt(z[6:8])  }, #pt4 in Ct
            #         {'type':'ineq', 'fun': lambda z: topCb(z[4:6])  }, #pt3 upper half of Ct
            #         {'type':'ineq', 'fun': lambda z: topCb(z[6:8])  }, #pt4 upper half of Ct
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[8], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[8], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[8])  }, # 0<= lambda <= 1
            #         {'type':'ineq', 'fun': lambda z: segmentInReg3(z[9], z[0:2], z[2:4]) }, #segment from pt1 to pt3 in reg3
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[9])  }, # 0<= mu <= 1
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[10], z[2:4], z[4:6])  }, #segment from pt2 to pt3 not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[10], z[2:4], z[4:6])  }, #segment from pt2 to pt3 not in the top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[10])  }, # 0<= alpha <= 1
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[11], z[4:6], z[6:8])  }, #segment from pt3 to pt4 not in reg3
            #         {'type':'ineq', 'fun': lambda z: inTopCircle(z[11], z[4:6], z[6:8])  }, #segment from pt3 to pt4 in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[11])  }, # 0<= gamma <= 1
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[12], z[2:4], xhat)  }, #segment from pt2 to xhat not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[12], z[2:4], xhat)  }, #segment from pt2 to xhat not in the top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[12])  } # 0<= beta <= 1
            #         ]
            # t6_opt = minimize(f_t6, [0, -10, 5*sqrt(2), 5*sqrt(2), 5*sqrt(2), 5*sqrt(2), 0, 10*sqrt(2), 0.5, 0.5, 0.5, 0.5, 0.5], constraints = cons6, options={'maxiter':mIter} )
            # t6 = f_t6(t6_opt.x)
            # # TYPE 5 1->around 3 -> 3 -> 1
            # f_t7 = lambda z: norm( sb(x0, z[0:2]) ) + ArcLengthCb(z[0:2], z[2:4]) + 5*norm(sb(z[2:4], z[4:6])) + norm( sb(z[4:6], xhat) )
            # cons7 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
            #         {'type':'eq', 'fun': lambda z: onCb(z[2:4] ) }, #pt2 on Cb
            #         {'type':'eq', 'fun': lambda z: onCb(z[4:6] ) }, #pt2 on Cb
            #         {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
            #         {'type':'ineq', 'fun': lambda z: bottomCb(z[2:4])  }, #pt2 not in arch
            #         {'type':'ineq', 'fun': lambda z: bottomCb(z[4:6])  }, #pt3 not in arch
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[6], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[6], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[6])  }, # 0<= lambda <= 1
            #         {'type':'ineq', 'fun': lambda z: segmentInReg3(z[7], z[0:2], z[2:4]) }, #segment from pt1 to pt2 in reg3
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[7])  }, # 0<= mu <= 1
            #         {'type':'ineq', 'fun': lambda z: segmentInReg3(z[8], z[2:4], z[4:6]) }, #segment from pt2 to pt3 in reg3
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[8])  } # 0<= alpha <= 1
            #         ]
            # t7_opt = minimize(f_t7, [0, -10, 10, 0, 5*sqrt(2), 5*sqrt(2), 0.5, 0.5, 0.5], constraints = cons7, options={'maxiter':mIter} )
            # t7 = f_t7(t7_opt.x)
            # t8_opt = minimize(f_t7, [-10, 0, 0, -10, 10, 0, 0, 0, 0], constraints = cons7, options={'maxiter':mIter} )
            # t8 = f_t7(t8_opt.x)
            # # 1 -> around 3 to p3-> 2 -> 1
            # f_t9 = lambda z: norm(x0, z[0:2]) + ArcLengthCb(z[0:2], p3) + 2*norm(sb(p3, z[2:4])) + norm(sb(z[2:4], xhat))
            # cons9 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
            #         {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[4])  }, # 0<= lambda <= 1
            #         {'type':'ineq', 'fun': lambda z: segmentInReg3(z[5], p3, z[0:2]) }, #segment from pt1 to p3 in reg3
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[5])  }, # 0<= mu <= 1
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[6], p3, z[2:4])  }, #segment from p3 to pt2 not in reg3
            #         {'type':'ineq', 'fun': lambda z: inTopCircle(z[6], p3, z[2:4])  }, #segment from p3 to pt2 in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[6])  }, # 0<= alpha <= 1
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[7], z[2:4], xhat) }, # segment from pt2 to xhat not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[7], z[2:4], xhat) }, # segment from pt2 to xhat not in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[7])  } # 0<= beta <= 1
            #         ]
            # t9_opt = minimize(f_t9, [0, -10, 0, 10*sqrt(2), 0.5, 0.5, 0.5, 0.5], constraints = cons9, options = {'maxiter': mIter})
            # t9 = f_t9(t9_opt.x)
            # # 1 -> around 3 to p4-> 2 -> 1
            # f_t10 = lambda z: norm(x0, z[0:2]) + ArcLengthCb(z[0:2], p4) + 2*norm(sb(p4, z[2:4])) + norm(sb(z[2:4], xhat))
            # cons9 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
            #         {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[4], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[4])  }, # 0<= lambda <= 1
            #         {'type':'ineq', 'fun': lambda z: segmentInReg3(z[5], p4, z[0:2]) }, #segment from pt1 to p4 in reg3
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[5])  }, # 0<= mu <= 1
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[6], p4, z[2:4])  }, #segment from p4 to pt2 not in reg3
            #         {'type':'ineq', 'fun': lambda z: inTopCircle(z[6], p4, z[2:4])  }, #segment from p4 to pt2 in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[6])  }, # 0<= alpha <= 1
            #         {'type':'ineq', 'fun': lambda z: notInSegment3(z[7], z[2:4], xhat) }, # segment from pt2 to xhat not in reg3
            #         {'type':'ineq', 'fun': lambda z: notInTopCircle(z[7], z[2:4], xhat) }, # segment from pt2 to xhat not in top circle
            #         {'type':'ineq', 'fun': lambda z: parametrization(z[7])  } # 0<= beta <= 1
            #         ]
            # t10_opt = minimize(f_t10, [0, -10, 0, 10*sqrt(2), 0.5, 0.5, 0.5, 0.5], constraints = cons9, options = {'maxiter': mIter})
            # t10 = f_t10(t10_opt.x)
            # 1->around 3-> around 2 -> 1
            f_t11 = lambda z: norm( sb(x0, z[0:2]) ) + ArcLengthCb(z[0:2], z[2:4]) + ArcLengthCt(z[2:4], z[4:6]) + norm( sb(z[4:6], xhat) )
            cons11 = [{'type':'eq', 'fun': lambda z: onCb(z[0:2] ) }, #pt1 on Cb
                    {'type':'eq', 'fun': lambda z: onCb(z[2:4] ) }, #pt2 on Cb
                    {'type':'eq', 'fun': lambda z: onCt(z[4:6] ) }, #pt3 on Ct
                    {'type':'ineq', 'fun': lambda z: bottomCb(z[0:2])  }, #pt1 not in arch
                    {'type':'ineq', 'fun': lambda z: bottomCb(z[2:4])  }, #pt2 not in arch
                    {'type':'ineq', 'fun': lambda z: topCb(z[2:4])  }, #pt2 top half of upper circle
                    {'type':'ineq', 'fun': lambda z: notInSegment3(z[6], x0, z[0:2]) }, # segment from x0 to pt1 not in reg3
                    {'type':'ineq', 'fun': lambda z: notInTopCircle(z[6], x0, z[0:2]) }, # segment from x0 to pt1 not in top circle
                    {'type':'ineq', 'fun': lambda z: parametrization(z[6])  }, # 0<= lambda <= 1
                    {'type':'ineq', 'fun': lambda z: segmentInReg3(z[7], z[0:2], z[2:4]) }, #segment from pt1 to pt2 in reg3
                    {'type':'ineq', 'fun': lambda z: parametrization(z[7])  }, # 0<= mu <= 1
                    {'type':'ineq', 'fun': lambda z: notInSegment3(z[8], z[2:4], z[4:6]) }, # segment from pt2 to pt3 not in reg3
                    {'type':'ineq', 'fun': lambda z: inTopCircle(z[8], z[2:4], z[4:6]) }, # segment from pt2 to pt3 in top circle
                    {'type':'ineq', 'fun': lambda z: parametrization(z[8])  } # 0<= omega <= 1
                    ]
            t11_opt = minimize(f_t11, [0, -10, 5*sqrt(2),5*sqrt(2) , 0, 10*sqrt(2), 0.5, 0.5, 0.5], constraints = cons11, options={'maxiter':mIter} )
            t11 = f_t11(t11_opt.x)
            possible_values = [t3, t4, t5, t11]
            eik = min(possible_values)
            type_route = possible_values.index(eik) + 1
    else:
        eik = 50
        type_route = -1
    return eik, type_route

nx = int(18*50)
ny = int(21*50) 
x, y = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
solution = np.zeros(x.shape)
solution_type = np.zeros(x.shape)
for i in range(ny):
    for j in range(nx):
        xhat = np.array([ x[i,j], y[i,j] ])
        sol, route_taken = SnowSolution(xhat)
        solution[i, j] = sol
        solution_type[i,j] = route_taken
 
        
fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im1 = plt.imshow( solution, cmap = colormap2, extent=[-18,18,-18,24]  )
plt.title("Exact solution test")
ax = plt.gca()
ax.invert_yaxis()
plt.show(block = False)
plt.colorbar(im1)


fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
plt.axis('equal')
im2 = plt.imshow( solution_type, cmap = colormap2, extent=[-18,18,-18,24]  )
plt.title("Solution type test")
ax = plt.gca()
ax.invert_yaxis()
plt.show(block = False)
plt.colorbar(im2)


fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# Plot the surface.
surf = ax.plot_surface(x, y, solution, cmap=colormap2, linewidth=0, antialiased=False)
plt.show(block = False)

plt.show()
