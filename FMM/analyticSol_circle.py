# ANALYTIC SOLUTION TO THE CIRCLE
from math import sqrt, cos, sin, pi, atan
from scipy.optimize import minimize, minimize_scalar
import numpy as np
from cmath import asin
from matplotlib.patches import Arc
import matplotlib as plt


def regA1Bool(xi, yi, zx, zy, inner_angleSource, tau3, x0, center, R,):
    '''
    This function determines if xhat can be accessed directly from x0. This is in general for a source and a circle since
    we are just looking at the angles and distances. 
    '''
    # It is in regA1 if the angle from x0 to xhat is greater than the angle from x0 to either of tha tangent points (meaning
    # that xhat lies "outside" the tangent lines) or
    # if the distance from x0 to xhat is smaller than the distance from x0 to either of the tangent points (meaning that
    # xhat lies inside the cone from x0 to the circle)
    rtan = sqrt((zx - x0[0])**2 + (zy - x0[1])**2) # Distance from the source to either of the tangent points 
    rx0 = sqrt( (center[0]-x0[0])**2 + (center[1]-x0[1])**2 ) # Distance from the center of the circle to the source x0
    uc, vc = center[0]-x0[0]/rx0, center[1]-x0[1]/rx0 # direction from the center of the circle to the center of the circle (unitary)
    dirTarget_x, dirTarget_y = (xi - x0[0])/tau3, (yi - x0[1])/tau3 # direction from the source to the xHat (unitary)
    if abs(np.arccos(dirTarget_x*uc + dirTarget_y*vc)) >= abs(inner_angleSource) or tau3 <= rtan:
        ans = True
    else:
        ans = False
    return ans


def paramCircle(theta, x0, center, R):
    '''
    Parametrization of the boundary - the circle with radius R centered at center
    '''
    return R*cos(theta) + center[0], R*sin(theta) + center[1]

def pointsTangentFromSource(xSource, ySource, xCenter, yCenter, R):
    '''
    Function that returns the coordinates of the two points on the tangent lines of the circle centered at xCenter, yCenter with 
    radius R that go through xSource, ySouce and the inner angle
    '''
    rSource = sqrt( (xSource - xCenter)**2 + (ySource - yCenter)**2  ) # distance from the source point to the center of the circle
    angle_Source = np.arctan2( ySource - yCenter, xSource - xCenter ) # signed angle from (1, 0) to the source point
    thtan = np.arccos(R/rSource)
    zx_1, zy_1 = R*np.cos( angle_Source + thtan  ), R*np.sin( angle_Source + thtan  )
    zx_2, zy_2 = R*np.cos( angle_Source - thtan  ), R*np.sin( angle_Source - thtan  )
    inner_angleSource = pi - pi/2 - thtan # Angle between ray from Source to center and the tangents (since the inner angles of a triangle must add up pi)
    return zx_1, zy_1, zx_2, zy_2, angle_Source,thtan, inner_angleSource

def insideTwoSegmentLine(xi, yi, angle_Source, thtan, x0, center, R, eta1, eta2, eps):
    '''
    Minimization problem. The target is inside the circle, we need to find a point on the boundary that is 
    directly accessible from the source. The path taken from the source to the target is made up from
    two segments of straight lines (Snell's law)
    '''
    def f_t2(theta):
        p_x, p_y = paramCircle(theta, x0, center, R)
        return eta1*sqrt( (x0[0] - p_x)**2 + (x0[1] - p_y)**2  ) + eta2*sqrt( (p_x - xi)**2 + (p_y - yi)**2  )
    bounds = (angle_Source - thtan, angle_Source + thtan)
    opt_theta = minimize_scalar(f_t2, None, bounds, method='bounded', tol=eps).x 
    tau_opt = f_t2(opt_theta)
    return tau_opt, opt_theta
    
def insideCreepingRay(xi, yi, zx_1, zy_1, zx_2, zy_2, thtan, angle_Source, x0, center, R, eta1, eta2, eps):
    '''
    In this case there is no minimization problem, we just need to calculate the tangent points from xhat to the circle
    outside the circle + reached by shed ray i.e. points outside the circle reached by rays which creep along the 
    circle before being shed tangentially into the shadow of the circle.
    '''
    bounds = (angle_Source + thtan, angle_Source - thtan + 2*pi)
    def fToOptimize1(theta):
        p_x, p_y = paramCircle(theta, x0, center, R)
        centralAngle = abs( 2*asin( (sqrt( (p_x-zx_1)**2 + (p_y-zy_1)**2 ))/(2*R) )  )
        # print('Angle from the (1,0) to one of the tangent points z1 or z2 :  ', angle_z)
        # print('Central angle from the tangent point to the optimum point p1: ', centralAngle)
        return eta1*sqrt(  (x0[0]-zx_1)**2 + (x0[1]-zy_1)**2 ) + eta1*R*centralAngle + eta2*sqrt(  (xi-p_x)**2 + (yi-p_y)**2 )
    opt_theta1 = minimize_scalar(fToOptimize1,  None, bounds, method='bounded', tol=eps).x 
    tau_opt1 = fToOptimize1(opt_theta1)
    def fToOptimize2(theta):
        p_x, p_y = paramCircle(theta, x0, center, R)
        centralAngle = abs( 2*asin( (sqrt( (p_x-zx_2)**2 + (p_y-zy_2)**2 ))/(2*R) )  )
        # print('Angle from the (1,0) to one of the tangent points z1 or z2 :  ', angle_z)
        # print('Central angle from the tangent point to the optimum point p1: ', centralAngle)
        return eta1*sqrt(  (x0[0]-zx_2)**2 + (x0[1]-zy_2)**2 ) + eta1*R*centralAngle + eta2*sqrt(  (xi-p_x)**2 + (yi-p_y)**2 )
    opt_theta2 = minimize_scalar(fToOptimize2,  None, bounds, method='bounded', tol=eps).x 
    tau_opt2 = fToOptimize2(opt_theta2)
    if(tau_opt1 < tau_opt2):
        tau_opt = tau_opt1
        opt_theta = opt_theta1
        angle_z = np.arctan(zy_1/zx_1)
    else:
        tau_opt = tau_opt2
        opt_theta = opt_theta2
        angle_z = np.arctan(zy_2/zx_2)
    return tau_opt, opt_theta, angle_z

def outsideShedRay(xi, yi, zx_1, zy_1, zx_2, zy_2, thtan, angle_Source, x0, center, R, eta1):
    '''
    In this case there is no minimization problem, we just need to calculate the tangent points from xhat to the circle
    outside the circle + reached by shed ray i.e. points outside the circle reached by rays which creep along the 
    circle before being shed tangentially into the shadow of the circle.
    '''
    n1 = sqrt( (xi - zx_1)**2 + (yi - zy_1)**2 ) # distance from p1 to the target
    n2 = sqrt( (xi - zx_2)**2 + (yi - zy_2)**2 ) # distance from p2 to the target
    # From which point its going to go around depends on which is closer (this makes sense because the speed of sound is piecewise constant)
    if ( n1 <= n2):
        zx, zy = zx_1, zy_1
    else:
        zx, zy = zx_2, zy_2
    angle_z = abs( angle_Source - thtan )
    px_1, py_1, px_2, py_2, angle_xHat,thtan_hat, inner_angleXhat = pointsTangentFromSource(xi, yi, center[0], center[1], R)
    if( sqrt((px_1 - zx)**2 + (py_1 - zy)**2) <=  sqrt((px_2 - zx)**2 + (py_2 - zy)**2)  ):
        px, py = px_1, py_1
    else:
        px, py = px_2, py_2
    centralAngle = abs( 2*asin( (sqrt( (px-zx)**2 + (py-zy)**2 ))/(2*R) )  )
    # print('Angle from the (1,0) to one of the tangent points z1 or z2:  ', angle_z)
    # print('Central angle from the tangent point to the optimum point p1: ', centralAngle)
    opt_theta = np.arctan2(py, px)
    tau5 = eta1*sqrt(  (x0[0] - zx)**2 + (x0[1] - zy)**2  ) + eta1*( centralAngle*R ) + eta1*sqrt(  (xi - px)**2 + (yi - py)**2  )
    return tau5, zx, zy, px, py, opt_theta

def outsideThroughCircle(xi, yi, angle_Source, thtan, x0, center, R, eta1, eta2):
    # the possible values for the first point on the circle are the same as before (those points on the circle that are direcly accessible from x0)
    bounds1 = (angle_Source - thtan, angle_Source + thtan )
    # The same idea is true for the second point but instead they are those points on the circle that are directly accessible from xhat
    zxHat_1, zyHat_1, zxHat_2, zyHat_2, angle_xHat, thtan_Hat, inner_anglexHat = pointsTangentFromSource(xi, yi, center[0], center[1], R)
    bounds2 = (angle_xHat - thtan_Hat, angle_xHat + thtan_Hat )
    # Set up the optimization problem
    def f_t4(theta_vec):
        theta1 = theta_vec[0]
        theta2 = theta_vec[1]
        px_1, py_1 = paramCircle(theta1, x0, center, R)
        px_2, py_2 = paramCircle(theta2, x0, center, R)
        return eta1*sqrt( (px_1 - x0[0])**2 + (py_1 - x0[1])**2 ) + eta2*sqrt( (px_1 - px_2)**2 + (py_1 - py_2)**2 ) + eta1*sqrt( (xi - px_2)**2 + (yi- py_2)**2 )
    t4_opti = minimize( f_t4, [angle_Source, angle_xHat], bounds= [ bounds1, bounds2  ] ).x
    tau4 = f_t4(t4_opti)
    theta1 = t4_opti[0]
    theta2 = t4_opti[1]
    px_1, py_1 = paramCircle(theta1, x0, center, R)
    px_2, py_2 = paramCircle(theta2, x0, center, R)
    return tau4, theta1, px_1, py_1, theta2, px_2, py_2

# def drawPathTaken(path_taken, type_path, x0, center, R,):
#     if(type_path == 1): # means the target is inside the circle and is reached by 2 segments of lines
#         plt.scatter(x0[0], x0[1], s=20, c='k', zorder=2)
#         plt.scatter(path_taken[0][2], path_taken[1][2], s=20, c='#0800ff', zorder=2)
#         plt.scatter(path_taken[0][1], path_taken[1][1], s=20, c='k', zorder=2)
#         plt.plot(path_taken[0], path_taken[1], c='k', linewidth=1, zorder=2)
#     if(type_path == 2): # means the target is inside the circle and it is reached by a creeping ray
#         plt.scatter(x0[0], x0[1], s=20, c='k', zorder=2)
#         plt.scatter(path_taken[0][2], path_taken[1][2], s=20, c='#0800ff', zorder=2)
#         plt.scatter(path_taken[0][1], path_taken[1][1], s=20, c='k', zorder=2)
#         ax = fig.gca()
#         ax.add_patch(Arc((center[0], center[1]), 2*R,  2*R, theta1 = 180*path_taken[2]/pi, theta2 = 180*path_taken[3]/pi, edgecolor="#000536", lw=1.5))


def trueSolution(xi, yi, x0, center, R, eta1, eta2, eps = np.finfo(np.float64).resolution, path = False):
    '''
    Analytic solution
    '''
    tau = np.inf
    opt_theta = np.inf
    type_path = -1
    zx_1, zy_1, zx_2, zy_2, angle_Source, thtan, inner_angleSource = pointsTangentFromSource(x0[0], x0[1], center[0], center[1], R)
    # print("First tangent point", zx_1, zy_1)
    # print("Second tangent point", zx_2, zy_2)
    # print("inner_angleSource", inner_angleSource)
    if( xi**2 + yi**2 <= R**2 ): #xHat is outside the circle
        # If this is the case then the target is INSIDE the circle, we have two options: reached byu creeping ray or reached with 2 segments of lines (Snell's law)
        tau_optOriginal, opt_thetaOriginal = insideTwoSegmentLine(xi, yi, angle_Source, thtan, x0, center, R, eta1, eta2, eps) # reached with 2 segments of straight lines
        tau_optCreepingRay, opt_thetaShedRay, angle_z = insideCreepingRay(xi, yi, zx_1, zy_1, zx_2, zy_2, thtan, angle_Source, x0, center, R, eta1, eta2, eps) # reached by going around the circle
        tau = min(tau_optOriginal, tau_optCreepingRay)
        if(tau_optOriginal < tau_optCreepingRay):
            tau = tau_optOriginal
            px_opt, py_opt = paramCircle(opt_thetaOriginal, x0, center, R)
            path_taken = [[x0[0], px_opt, xi], [x0[1], py_opt, xi]]
            type_path = 1
        else:
            tau = tau_optCreepingRay
            px_opt, py_opt = paramCircle(opt_thetaShedRay, x0, center, R)
            path_taken = [[x0[0], px_opt, xi], [x0[1], py_opt, xi], opt_thetaShedRay, angle_z]
            type_path = 2
    else: 
        tau3 = sqrt( (xi-x0[0])**2 + (yi - x0[1])**2 )
        if(regA1Bool(xi, yi, zx_1, zy_1, inner_angleSource, tau3, x0, center, R,)): # if xhat is directly accesible from x0 just via reg1
            tau = tau3
            type_path = 3
            path_taken = [[x0[0], xi], [x0[1], yi]]
        else: # if xhat is in x0 but the ray has to go trough reg3 to get to xhat
            tau4, theta1, px_1, py_1, theta2, px_2, py_2 = outsideThroughCircle(xi, yi, angle_Source, thtan, x0, center, R, eta1, eta2)
            tau5, zx, zy, px, py, opt_theta = outsideShedRay(xi, yi, zx_1, zy_1, zx_2, zy_2, thtan, angle_Source, x0, center, R, eta1)
            if (tau4 < tau5):
                tau = tau4
                type_path = 4
                path_taken = [ [x0[0], px_1, px_2, xi], [x0[1], py_1, py_2, yi] ]
            else:
                tau = tau5
                type_path = 5
                px, py = paramCircle(opt_theta, x0, center, R)
                path_taken = [ [x0[0], zx, px, xi], [x0[1], zy, py, yi], opt_theta ]
    if path:
        return tau, path_taken, type_path
    else:
        return tau, type_path
