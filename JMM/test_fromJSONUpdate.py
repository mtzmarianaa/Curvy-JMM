import intermediateTests as itt
import matplotlib.pyplot as plt


run stepWithPython "./updates/update10501.json"

self = triFan

itt.plotFann(self.x0, self.listB0k, self.listxk, self.listBk,
             params = self.optiParams, indCrTop = self.optiIndCrTop,
             paramsCrTop = self.optiParamsCrTop, indStTop = self.optiIndStTop,
             paramsStTop = self.optiParamsStTop, listBkBk1 = self.listBkBk1)

ax = plt.gca()
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)

plt.xlim([-18, 18])
plt.ylim([-18, 24])

plt.title("With T0:" + str(self.T0) + "  and T1: " + str(self.T1) + " get  THAT: " + str(self.opti_fVal) )

