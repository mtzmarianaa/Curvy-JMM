import intermediateTests as itt
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from analyticSol_circle import trueSolution
import json


plt.ion()

# Find which json file contains the desired update

def printWhichUpdate(xHat, path = "./updates/", maxF = 16560):
    '''
    Given json files print the name of those files whose xHat is the desired xHat
    '''
    jsonFiles = []
    for i in range(maxF):
        fName = path + "update" + str(i) + ".json"
        f = open(fName)
        params_dict = json.load(f)
        if( np.all( params_dict["xHat"] == xHat ) ):
            print("json file with desired xHat: ", i)
            jsonFiles.append(i)
    return jsonFiles
    



# Read necessary information
H = "H0"
path_general = '/Users/marianamartinez/Documents/Curvy-JMM/JMM/'
path_information = path_general + H + "/" + H
eik_vals = np.fromfile(path_information + "_ComputedValuesFast.bin")
eik_coords = np.genfromtxt(path_information + "_MeshPoints.txt", delimiter=",")
eik_grads = np.fromfile(path_information + "_ComputedGradientsFast.bin")
eik_grads = eik_grads.reshape(len(eik_coords), 2)
triangles_points = np.genfromtxt(path_information + "_Faces.txt", delimiter=",")


indHat = 517
xHat = eik_coords[indHat]
xSource = np.array([-15, -10])
eta1 = 1.0
eta2 = 1.45
center = np.array([0,0])
R = 10.0

# Compute the true solution so that we can plot it as well

trueEik, typeSol, trueGrad = trueSolution( xHat[0], xHat[1], xSource, center, R, eta1, eta2  )


self = triFan
trueEik, typeSol, trueGrad = trueSolution( self.xHat[0], self.xHat[1], xSource, center, R, eta1, eta2  )
relErr = abs(self.opti_fVal - trueEik)/trueEik
relAngle = np.arccos( np.dot(self.lastGrad, trueGrad)/(norm(self.lastGrad)*norm(trueGrad) ) )

trueEik, typeSol, trueGrad = trueSolution( self.xHat[0], self.xHat[1], xSource, center, R, eta1, eta2  )
trueEik0, typeSol0, trueGrad0 = trueSolution( self.x0[0], self.x0[1], xSource, center, R, eta1, eta2  )
trueEik1, typeSol1, trueGrad1 = trueSolution( self.x1[0], self.x1[1], xSource, center, R, eta1, eta2  )

itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)
plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
plt.quiver([self.x0[0], self.x1[0], self.xHat[0]], [self.x0[1], self.x1[1], self.xHat[1]],
           [self.grad0[0], self.grad1[0], self.lastGrad[0]], [self.grad0[1], self.grad1[1], self.lastGrad[1]],
           color = "#00fff7", alpha = 0.8, label = "Computed")

plt.plot( [self.xHat[0]- self.lastGrad[0], self.xHat[0]], [self.xHat[1]-  self.lastGrad[1],
                                                           self.xHat[1]], c = "#00fff7", linewidth = 0.3)

plt.scatter(self.path[0:2, 0], self.path[0:2, 1], marker = "*")

plt.quiver( [self.xHat[0], self.x0[0], self.x1[0]],
            [self.xHat[1], self.x0[1], self.x1[1]],
            [trueGrad[0], trueGrad0[0], trueGrad1[0]],
            [trueGrad[1], trueGrad0[1], trueGrad1[1]],
            color = "#ff6200", alpha = 0.8, label = "True")
plt.plot( [self.xHat[0]-trueGrad[0], self.xHat[0]], [self.xHat[1]-trueGrad[1], self.xHat[1]], c = "#ff6200", linewidth = 0.3) 
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title(" THAT: " + "{0:0.5f}".format(self.opti_fVal) + "  true: " + "{0:0.5f}".format(trueEik) +
          " relErrEik: " + "{0:0.5f}".format(relErr) + " angleErr: " + "{0:0.5f}".format(relAngle) +
          " JSON 1" )

plt.legend()




# Plot what was stored
ind0 = 101
ind1 = 102
x0 = eik_coords[ind0]
grad0 = eik_grads[ind0]
x1 = eik_coords[ind1]
grad1 = eik_grads[ind1]
# Compute the true solution and gradients of the actual parents
trueEik0, typeSol0, trueGrad0 = trueSolution( x0[0], x0[1], xSource, center, R, eta1, eta2  )
trueEik1, typeSol1, trueGrad1 = trueSolution( x1[0], x1[1], xSource, center, R, eta1, eta2  )

relAngle = np.arccos( np.dot(eik_grads[indHat], trueGrad)/(norm(eik_grads[indHat])*norm(trueGrad) ) )
relErr = abs(eik_vals[indHat] - trueEik)/trueEik

relErr0 = abs(eik_vals[ind0] - trueEik0)/trueEik0
relErr1 = abs(eik_vals[ind1] - trueEik1)/trueEik1

plt.figure()
plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
plt.scatter( xHat[0], xHat[1], marker = "*")
plt.scatter( eik_coords[ind0, 0], eik_coords[ind0, 1], marker = "o", label = "x0")
plt.scatter( eik_coords[ind1, 0], eik_coords[ind1, 1], marker = "s", label = "x1")
plt.quiver( xHat[0], xHat[1], eik_grads[indHat, 0], eik_grads[indHat, 1],
           color = "#00fff7", alpha = 0.8, label = "Computed")
plt.plot( [xHat[0]- eik_grads[indHat, 0], xHat[0]], [xHat[1]-  eik_grads[indHat, 1], xHat[1]],
          c = "#00fff7", linewidth = 0.3, label = "Computed")
plt.quiver( [xHat[0], x0[0], x1[0]],
            [xHat[1], x0[1], x1[1]],
            [trueGrad[0], trueGrad0[0], trueGrad1[0]],
            [trueGrad[1], trueGrad0[1], trueGrad1[1]],
            color = "#ff6200", alpha = 0.8, label = "True")
plt.plot( [xHat[0]-trueGrad[0], xHat[0]], [xHat[1]-trueGrad[1], xHat[1]],
          c = "#ff6200", linewidth = 0.3, label = "True") 
ax = plt.gca()
ax.set_aspect("equal")
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.xlim([-18, 18])
plt.ylim([-18, 24])
plt.legend()
plt.title(" THAT: " + "{0:0.5f}".format(eik_vals[indHat]) + "  true: " + "{0:0.5f}".format(trueEik) + " relErrEik: " + "{0:0.5f}".format(relErr) + " angleErr: " + "{0:0.5f}".format(relAngle) + " ACCEPTED" )









# IF THE MARCHER USED A SIMPLE TWO POINT UPDATE USING C TRY BUILDING THAT TRIANGLE FAN
# AND SOLVE IT USING PYTHON - SEE IF WE GET SIMILAR RESUTLS

indInside = eta2
listxk = [ list(x0), list(x1), list(xHat) ]
listB0k = 2*[ [0.0, 0.0] ]
listBk = 2*[ [0.0, 0.0] ]
listBkBk1 = 2*[ [0.0, 0.0] ]

triInfo = '{"x0": '+ str(list(x0)) + ', "T0": '+ str(eik_vals[ind0]) + ', "grad0": '+ str(list(eik_grads[ind0])) + ', "x1": '+ str(list(x1)) + ', "T1": '+ str(eik_vals[ind1]) + ', "grad1": '+ str(list(eik_grads[ind1])) + ', "xHat": '+ str(list(xHat)) +', "listIndices": '+ str(3*[indInside]) +', "listxk": '+ str(listxk) + ', "listB0k": '+ str(listB0k)+', "listBk": '+ str(listBk) +', "listBkBk1": '+ str(listBkBk1) +', "plotBefore": 0, "plotAfter": 0, "plotOpti": 0 } '

params_dict = json.loads(triInfo)
nRegions = len(params_dict['listxk']) -2
triFan = oP.triangleFan(nRegions)
dict_out = triFan.outputJSON(triInfo)
self = triFan
indHat = 519
xHat = eik_coords[indHat]
xSource = np.array([-15, -10])
eta1 = 1.0
eta2 = 1.45
center = np.array([0,0])
R = 10.0

# Compute the true solution so that we can plot it as well

trueEik, typeSol, trueGrad = trueSolution( xHat[0], xHat[1], xSource, center, R, eta1, eta2  )
relErr = abs(self.opti_fVal - trueEik)/trueEik
relAngle = np.arccos( np.dot(self.lastGrad, trueGrad)/(norm(self.lastGrad)*norm(trueGrad) ) )

trueEik0, typeSol0, trueGrad0 = trueSolution( self.x0[0], self.x0[1], xSource, center, R, eta1, eta2  )
trueEik1, typeSol1, trueGrad1 = trueSolution( self.x1[0], self.x1[1], xSource, center, R, eta1, eta2  )

itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)
plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
plt.quiver([self.x0[0], self.x1[0], self.xHat[0]], [self.x0[1], self.x1[1], self.xHat[1]],
           [self.grad0[0], self.grad1[0], self.lastGrad[0]], [self.grad0[1], self.grad1[1], self.lastGrad[1]],
           color = "#00fff7", alpha = 0.8, label = "Computed")

plt.plot( [self.xHat[0]- self.lastGrad[0], self.xHat[0]], [self.xHat[1]-  self.lastGrad[1],
                                                           self.xHat[1]], c = "#00fff7", linewidth = 0.3)

plt.scatter(self.path[0:2, 0], self.path[0:2, 1], marker = "*")

plt.quiver( [self.xHat[0], self.x0[0], self.x1[0]],
            [self.xHat[1], self.x0[1], self.x1[1]],
            [trueGrad[0], trueGrad0[0], trueGrad1[0]],
            [trueGrad[1], trueGrad0[1], trueGrad1[1]],
            color = "#ff6200", alpha = 0.8, label = "True")
plt.plot( [self.xHat[0]-trueGrad[0], self.xHat[0]], [self.xHat[1]-trueGrad[1], self.xHat[1]], c = "#ff6200", linewidth = 0.3) 
ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title(" THAT: " + "{0:0.5f}".format(self.opti_fVal) + "  true: " + "{0:0.5f}".format(trueEik) +
           " absErrEik: " + "{0:0.5f}".format(abs(self.opti_fVal - trueEik)) + " relErrEik: " +
          "{0:0.5f}".format(relErr) + " angleErr: " + "{0:0.5f}".format(relAngle) +
          " opti with python" )

plt.legend()






# Change parameters

newParams = np.copy(self.optiParams)
newParams[1] = 0
newParamsCrTop = [0.1, 0.9]
newIndCrTop = [1]

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

plt.triplot( eik_coords[:, 0], eik_coords[:, 1], triangles_points, '-', c = "#d4bdff", lw = 0.3 )
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

relErr = abs(f_test - trueEik)/trueEik
relAngle = np.arccos( np.dot(gradLast, trueGrad)/(norm(gradLast)*norm(trueGrad) ) )

plt.title(" THAT: " + "{0:0.5f}".format(f_test) + "  true: " + "{0:0.5f}".format(trueEik) + " relErrEik: " + "{0:0.5f}".format(relErr) + " angleErr: " + "{0:0.5f}".format(relAngle) + " params change" )





