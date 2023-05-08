
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import sqrt, pi, cos, sin
import intermediateTests as itt

import matplotlib.colors as clr
from mpl_toolkits.mplot3d import axes3d


colormap2 = "cet_linear_worb_100_25_c53_r"
colormap2_r = "cet_linear_worb_100_25_c53"


r = 10

N = 51
N_interp = 50


theta = np.linspace(0, 2*np.pi, N, endpoint=False)

x = r*np.cos(theta)
y = r*np.sin(theta)

tx = -r*np.sin(theta)
ty = r*np.cos(theta)

dx = x[1:] - x[:-1]
dy = y[1:] - y[:-1]

h = np.sqrt(dx**2 + dy**2)
h = h[0] # Because they are equispaced, this is just norm(xk - xkM1)


p = lambda t: np.array([1, t, t**2, t**3])
dp = lambda t: np.array([0, 1, 2*t, 3*t**2])

V = lambda h: np.array([p(0), p(h), dp(0), dp(h)])

c00, c10, c01, c11 = np.linalg.inv(V(h)).T



T = np.linspace(0, h, 50) # Between each node compute 50 points

x_inter = np.empty(len(x)*N_interp)
y_inter = np.empty(len(y)*N_interp)
x_cord = np.empty(len(x)*N_interp)
y_cord = np.empty(len(y)*N_interp)
89
for i in range(-1, N-1):
    # Set params
    xFrom = np.array([x[i], y[i]])
    Bfrom = np.array([tx[i], ty[i]])
    Bfrom = Bfrom*h/norm(Bfrom)
    xTo = np.array([x[i+1], y[i+1]])
    Bto = np.array([tx[i+1], ty[i+1]])
    Bto = Bto*h/norm(Bto)
    X = lambda t: xFrom[0]*c00@p(t) + xTo[0]*c10@p(t) + tx[i]*c01@p(t) + tx[i+1]*c11@p(t)
    Y = lambda t: xFrom[1]*c00@p(t) + xTo[1]*c10@p(t) + ty[i]*c01@p(t) + ty[i+1]*c11@p(t)
    # Interpolate in that interval
    k = i
    for j in range(0, N_interp):
        X_interp = itt.hermite_boundary(T[j], xFrom, Bfrom, xTo, Bto)
        x_inter[k] = X_interp[0]
        y_inter[k] = X_interp[1]
        x_cord[k] = X(T[j])
        y_cord[k] = Y(T[j])
        k += 1



# Plot

plt.figure(figsize=(800/96, 800/96), dpi=96)
fig = plt.gcf()
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#00e4ff",fill=False)
ax.add_patch(circle_b)
plt.plot(x_inter, y_inter, c = "black", linewidth = 0.5, label = "Mine")
plt.plot(x_cord, y_cord, c = "#00ff6c", linewidth = 0.5, label = "Chordial")
plt.scatter(x, y, marker = "*", c = "#2c00ff")
plt.ylim([-10,10])
plt.xlim([-10,10])


        
