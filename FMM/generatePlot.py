import matplotlib.pyplot as plt
import numpy as np
colormap = plt.cm.get_cmap('plasma')
x = np.fromfile('Outputs/data.bin')
x = x.reshape(200, 200)
colors = colormap(x)
sm = plt.cm.ScalarMappable(cmap=colormap)
plt.figure()
plt.imshow(x, cmap=colormap)
plt.colorbar(sm)
plt.show()