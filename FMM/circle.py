
import colorcet as cc
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
from math import pi


eta0 = 1
eta1 = 4

r = 3
R = 2

xsrc, ysrc = -2.5, -2.5
xhat, yhat = -0.75, 0

rsrc = np.sqrt(xsrc**2 + ysrc**2)
thsrc = np.arctan2(ysrc, xsrc)
thtan = np.arccos(R/rsrc)

xtanp, ytanp = R*np.cos(thsrc + thtan), R*np.sin(thsrc + thtan)
xtanm, ytanm = R*np.cos(thsrc - thtan), R*np.sin(thsrc - thtan)

uc, vc = -xsrc/rsrc, -ysrc/rsrc
rtan = np.sqrt((xtanp - xsrc)**2 + (ytanp - ysrc)**2)
utanp, vtanp = (xtanp - xsrc)/rtan, (ytanp - ysrc)/rtan

# Angle between ray from (xsrc, ysrc) to center of circle and tangents
thtansrc = np.arccos(uc*utanp + vc*vtanp)
print(thtansrc)
print(pi - pi/2 - thtan)

Theta = np.linspace(0, 2*np.pi)
Xcirc, Ycirc = R*np.cos(Theta), R*np.sin(Theta)

def get_opt_on_circ(x, y):
    def f(theta):
        xth, yth = R*np.cos(theta), R*np.sin(theta)
        return eta0*np.sqrt((xth - xsrc)**2 + (yth - ysrc)**2) \
            + eta1*np.sqrt((x - xth)**2 + (y - yth)**2)
    bounds = (thsrc - thtan, thsrc + thtan)
    # import ipdb; ipdb.set_trace()
    eps = np.finfo(np.float64).resolution
    res = scipy.optimize.minimize_scalar(f, None, bounds, method='bounded', tol=eps)
    thopt = res.x
    return R*np.cos(thopt), R*np.sin(thopt), thopt, f(thopt)

N = 129
X, Y = np.meshgrid(np.linspace(-r, r, N),
                   np.linspace(-r, r, N),
                   indexing='xy')
Tau = np.empty_like(X)
Tau[...] = np.nan
for i, j in it.product(range(N), range(N)):
    x, y, tau = X[i, j], Y[i, j], np.inf

    taudirect = np.sqrt((x - xsrc)**2 + (y - ysrc)**2)
    u, v = (x - xsrc)/taudirect, (y - ysrc)/taudirect

    if x**2 + y**2 <= R**2:
        xopt, yopt, thopt, tauopt = get_opt_on_circ(x, y)
        tau = min(tau, tauopt)
    elif abs(np.arccos(u*uc + v*vc)) >= thtansrc or taudirect <= rtan:
        tau = min(tau, taudirect)

    Tau[i, j] = tau

xopt, yopt, thopt, tauopt = get_opt_on_circ(xhat, yhat)

plt.figure()
plt.contourf(X, Y, Tau, levels=11, cmap=cc.cm.blues)
plt.colorbar()
plt.plot(Xcirc, Ycirc, c='k', linewidth=1, zorder=2)
plt.plot([xsrc, xsrc + 20*(xtanp - xsrc)], [ysrc, ysrc + 20*(ytanp - ysrc)],
         c='k', linewidth=1, linestyle='--', zorder=2)
plt.plot([xsrc, xsrc + 20*(xtanm - xsrc)], [ysrc, ysrc + 20*(ytanm - ysrc)],
         c='k', linewidth=1, linestyle='--', zorder=2)
plt.scatter(xtanp, ytanp, s=20, c='k', zorder=2)
plt.scatter(xtanm, ytanm, s=20, c='k', zorder=2)
plt.scatter(xsrc, ysrc, s=20, c='k', zorder=2)
plt.scatter(xhat, yhat, s=20, c='k', zorder=2)
plt.scatter(xopt, yopt, s=20, c='k', zorder=2)
plt.plot([xsrc, xopt, xhat], [ysrc, yopt, yhat], c='k', linewidth=1, zorder=2)
plt.xlim(-r, r)
plt.ylim(-r, r)
plt.gca().set_aspect('equal')
plt.show()
