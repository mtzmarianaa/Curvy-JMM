import intermediateTests as itt
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from analyticSol_circle import trueSolution


# Find which json file contains the desired update

def printWhichUpdate(xHat, path = "./updates/", maxF = 16560):
    '''
    Given json files print the name of those files whose xHat is the desired xHat
    '''
    for i in range(maxF):
        fName = path + "update" + str(i) + ".json"
        f = open(fName)
        params_dict = json.load(f)
        if( np.all( params_dict["xHat"] == xHat ) ):
            print("json file with desired xHat: ", i)
    

plt.ion()


# Read necessary information
H = "H4"
path_general = '/Users/marianamartinez/Documents/Curvy-JMM/JMM/'
path_information = path_general + H + "/" + H
eik_vals = np.fromfile(path_information + "_ComputedValuesFast.bin")
eik_coords = np.genfromtxt(path_information + "_MeshPoints.txt", delimiter=",")
eik_grads = np.fromfile(path_information + "_ComputedGradientsFast.bin")
eik_grads = eik_grads.reshape(len(eik_coords), 2)
triangles_points = np.genfromtxt(path_information + "_Faces.txt", delimiter=",")


self = triFan
xSource = np.array([-15, -10])
eta1 = 1.0
eta2 = 1.452
center = np.array([0,0])
R = 10.0

# Compute the true solution so that we can plot it as well

trueEik, typeSol, trueGrad = trueSolution( self.xHat[0], self.xHat[1], xSource, center, R, eta1, eta2  )
relErr = abs(self.opti_fVal - trueEik)/trueEik
relAngle = np.arccos( np.dot(self.lastGrad, trueGrad)/(norm(self.lastGrad)*norm(trueGrad) ) )

trueEik0, typeSol0, trueGrad0 = trueSolution( self.x0[0], self.x0[1], xSource, center, R, eta1, eta2  )
trueEik1, typeSol1, trueGrad1 = trueSolution( self.x1[0], self.x1[1], xSource, center, R, eta1, eta2  )

itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)

plt.quiver([self.x0[0], self.x1[0], self.xHat[0]], [self.x0[1], self.x1[1], self.xHat[1]],
           [self.grad0[0], self.grad1[0], self.lastGrad[0]], [self.grad0[1], self.grad1[1], self.lastGrad[1]])
plt.scatter(self.path[0:2, 0], self.path[0:2, 1], marker = "*")

plt.quiver( [self.xHat[0], self.x0[0], self.x1[0]],
            [self.xHat[1], self.x0[1], self.x1[1]],
            [trueGrad[0], trueGrad0[0], trueGrad1[0]],
            [trueGrad[1], trueGrad0[1], trueGrad1[1]], color = "#ff6200")
plt.plot( [self.xHat[0]-trueGrad[0], self.xHat[0]], [self.xHat[1]-trueGrad[1], self.xHat[1]], c = "#c84d00", linewidth = 0.3) 
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title(" THAT: " + "{0:0.5f}".format(self.opti_fVal) + "  true: " + "{0:0.5f}".format(trueEik) + " relErrEik: " + "{0:0.5f}".format(relErr) + " angleErr: " + "{0:0.5f}".format(relAngle) )






# Change parameters

newParams = np.copy(self.optiParams)
# newParams[0] = 0.7
# newParams[1] = 0.5
# newParams[2] = 0.5
newParams[3] = 1.0

f_test = fObj_generalized(newParams, self.x0, self.T0, self.grad0,
                            self.x1, self.T1, self.grad1, self.xHat,
                            self.listIndices, self.listxk, self.listB0k,
                            self.listBk, self.listBkBk1,
                            self.optiIndCrTop, self.optiParamsCrTop,
                            self.optiIndStTop, self.optiParamsStTop)

itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = newParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)


newPath, newGrads = getPathGradEikonal(newParams, self.listIndices,
                                       self.listxk, self.listB0k,
                                       self.listBk, self.listBkBk1,
                                       self.optiIndCrTop, self.optiParamsCrTop,
                                       self.optiIndStTop, self.optiParamsStTop)

gradLast = newGrads[-1]
k = len(newGrads) - 1
while( norm(gradLast) == 0 and k > -1):
    gradLast = newGrads[k]
    k = k -1

relAngle = np.arccos( np.dot(gradLast, trueGrad)/(norm(gradLast)*norm(trueGrad) ) )
plt.quiver([self.x0[0], self.x1[0], self.xHat[0]], [self.x0[1], self.x1[1], self.xHat[1]],
           [self.grad0[0], self.grad1[0], gradLast[0]], [self.grad0[1], self.grad1[1], gradLast[1]])
plt.scatter(newPath[0:2, 0], newPath[0:2, 1], marker = "*")

plt.quiver( [self.xHat[0], self.x0[0], self.x1[0]],
            [self.xHat[1], self.x0[1], self.x1[1]],
            [trueGrad[0], trueGrad0[0], trueGrad1[0]],
            [trueGrad[1], trueGrad0[1], trueGrad1[1]], color = "#ff6200")
plt.plot( [self.xHat[0]-trueGrad[0], self.xHat[0]], [self.xHat[1]-trueGrad[1], self.xHat[1]], c = "#c84d00", linewidth = 0.3) 
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title(" THAT: " + "{0:0.5f}".format(f_test) + "  true: " + "{0:0.5f}".format(trueEik) + " relErrEik: " + "{0:0.5f}".format(relErr) + " angleErr: " + "{0:0.5f}".format(relAngle) )



