

import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

#plt.ion()

eta0 = 1
eta1 = 10

r = 3
R = 2

xsrc, ysrc = -2.5, -2.5
xhat, yhat = 0, -1

Theta = np.linspace(0, 2*np.pi)
Xcirc, Ycirc = R*np.cos(Theta), R*np.sin(Theta)

def get_opt_on_circ(x, y):
    def f(theta):
        xth, yth = R*np.cos(theta), R*np.sin(theta)
        return eta0*np.sqrt((xth - xsrc)**2 + (yth - ysrc)**2) \
            + eta1*np.sqrt((x - xth)**2 + (y - yth)**2)
    theta_opt = scipy.optimize.minimize_scalar(f, None, None).x
    return R*np.cos(theta_opt), R*np.sin(theta_opt), f(theta_opt)

xopt, yopt, tau_opt = get_opt_on_circ(xhat, yhat)

N = 65
X, Y = np.meshgrid(np.linspace(-r, r, N),
                   np.linspace(-r, r, N),
                   indexing='xy')
Tau = np.empty_like(X)
for i, j in it.product(range(N), range(N)):
    x, y = X[i, j], Y[i, j]
    tau_direct = eta0*np.sqrt((x - xsrc)**2 + (y - ysrc)**2)
    x_opt, y_opt, tau_opt = get_opt_on_circ(x, y)

    # 1. if tau_direct ray doesn't intersect circle, use tau_direct
    # 2. if tangent vectors on boundary point into circle, use tau_opt
    # 3. are there are any cases we've missed at this point?

    Tau[i, j] = tau_direct

plt.figure()
plt.contourf(X, Y, Tau, levels=11)
plt.colorbar()
plt.plot(Xcirc, Ycirc, c='k', linewidth=1, zorder=2)
plt.scatter(xsrc, ysrc, s=20, c='k', zorder=2)
plt.scatter(xhat, yhat, s=20, c='k', zorder=2)
plt.scatter(xopt, yopt, s=20, c='k', zorder=2)
plt.plot([xsrc, xopt, xhat], [ysrc, yopt, yhat], c='k', linewidth=1, zorder=2)
plt.xlim(-r, r)
plt.ylim(-r, r)
plt.gca().set_aspect('equal')
plt.show()
