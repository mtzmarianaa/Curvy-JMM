import intermediateTests as itt
import matplotlib.pyplot as plt
from numpy.linalg import norm


run stepWithPython "./updates/update1132.json"

self = triFan
xSource = np.array([-15, -10])


itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)

plt.plot([self.xHat[0], xSource[0]], [self.xHat[1], xSource[1]])
plt.quiver([self.x0[0], self.x1[0]], [self.x0[1], self.x1[1]], [self.grad0[0], self.grad1[0]], [self.grad0[1], self.grad1[1]])
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)

plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title("With T0:" + str(self.T0) + "  and T1: " + str(self.T1) + " get  THAT: " + str(self.opti_fVal) + "  true: " + str(norm(xSource - self.xHat)) )






run stepWithPython "./updates/update1113.json"

self = triFan
xSource = np.array([-15, -10])


itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)

plt.plot([self.xHat[0], xSource[0]], [self.xHat[1], xSource[1]])
plt.quiver([self.x0[0], self.x1[0]], [self.x0[1], self.x1[1]], [self.grad0[0], self.grad1[0]], [self.grad0[1], self.grad1[1]])
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)

plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title("With T0:" + str(self.T0) + "  and T1: " + str(self.T1) + " get  THAT: " + str(self.opti_fVal) + "  true: " + str(norm(xSource - self.xHat)) )




run stepWithPython "./updates/update1000.json"

self = triFan
xSource = np.array([-15, -10])


itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)

plt.plot([self.xHat[0], xSource[0]], [self.xHat[1], xSource[1]])
plt.plot([self.x0[0], xSource[0]], [self.x0[1], xSource[1]], linestyle='dashed')
plt.plot([self.x1[0], xSource[0]], [self.x1[1], xSource[1]], linestyle='dashed')
plt.quiver([self.x0[0], self.x1[0], self.xHat[0]], [self.x0[1], self.x1[1], self.xHat[1]],
           [self.grad0[0], self.grad1[0], self.lastGrad[0]], [self.grad0[1], self.grad1[1], self.lastGrad[1]])
plt.scatter(self.path[0:2, 0], self.path[0:2, 1], marker = "*")


plt.title("With T0:" + str(self.T0) + "  and T1: " + str(self.T1) + " get  THAT: " + str(self.opti_fVal) + "  true: " + str(norm(xSource - self.xHat)) )






params2 = np.copy(self.optiParams)
params2[0] = 0.3
params2[1] = 1.0
fTest = fObj_generalized(params2, self.x0, self.T0, self.grad0,
                         self.x1, self.T1, self.grad1, self.xHat,
                         self.listIndices, self.listxk, self.listB0k,
                         self.listBk, self.listBkBk1,
                         self.optiIndCrTop, self.optiParamsCrTop,
                         self.optiIndStTop, self.optiParamsStTop)


itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = params2, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)

plt.plot([self.xHat[0], xSource[0]], [self.xHat[1], xSource[1]])
plt.plot([self.x0[0], xSource[0]], [self.x0[1], xSource[1]], linestyle='dashed')
plt.plot([self.x1[0], xSource[0]], [self.x1[1], xSource[1]], linestyle='dashed')

plt.quiver([self.x0[0], self.x1[0]], [self.x0[1], self.x1[1]], [self.grad0[0], self.grad1[0]], [self.grad0[1], self.grad1[1]])


plt.title("Trying new parameters, fFound:" + str(self.opti_fVal) + "  new value found:" + str(fTest) )










self = triFan
xSource = np.array([-15, -10])


itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)

plt.quiver([self.x0[0], self.x1[0], self.xHat[0]], [self.x0[1], self.x1[1], self.xHat[1]],
           [self.grad0[0], self.grad1[0], self.lastGrad[0]], [self.grad0[1], self.grad1[1], self.lastGrad[1]])
plt.scatter(self.path[0:2, 0], self.path[0:2, 1], marker = "*")
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)

plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title("With T0:" + str(self.T0) + "  and T1: " + str(self.T1) + " get  THAT: " + str(self.opti_fVal) + "  true: " + str(norm(xSource - self.xHat)) )
