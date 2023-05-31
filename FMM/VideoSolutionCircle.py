import matplotlib.pyplot as plt
import numpy as np
from math import pi
import matplotlib.tri as tri
import random
import os
import cv2

plt.ion()
saveFigures = True
nx, ny = 36*10, 42*10
my_dpi = 96
fps = 25
Hs = "H21"
Amp = 1    # wave's amplitude
T = 1   # wave's period
omega = 2*pi/T   # wave's angular frequency
indicesRef = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1,
              1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 
              1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2,
              2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 3.4, 2.45, 2.5,
              2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3]
path_to_figures = '/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestBaseSnow/' + Hs + "/" + Hs 
indicesRef = [1.452*t for t in indicesRef] #get the actual different indices of refraction
nIndicesRef = len(indicesRef)
nameFilesEnd = ["_ind"+str(i) + ".bin" for i in range(nIndicesRef)]
colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)

# we read the coordinates and the triangles
coords = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + Hs + "/" + Hs + "_MeshPoints.txt", delimiter=",")
nPoints = len(coords)
triangles = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + Hs + "/" + Hs + "_Faces.txt", delimiter=",")
triang = tri.Triangulation(coords[:, 0], coords[:, 1], triangles)

# for the quiver plot we set the random seed - we don't want to plot all the gradients because they are too many, we would see anything
random.seed(a = 166297, version=2)
indicesGradients = random.sample(range(nPoints), round(nPoints/50))

# grid to plot
xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
true_solGrid = np.zeros(xi.shape)
type_solution = np.zeros(xi.shape)

figuresLevelSet = []
figuresWaveReal = []
figuresWaveImaginary = []
figuresGradients = []

# Generate the plots (both the quiver and the linear interpolation plots)

for i in range(nIndicesRef):
    # We first read the eikonal and the gradients
    # eik_vals = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + Hs + "/" +Hs + "_ComputedValues" + nameFilesEnd[i])
    # eik_grads = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/" + Hs + "/" +Hs +  "_ComputedGradients" + nameFilesEnd[i])
    # eik_grads = eik_grads.reshape(nPoints, 2) # reshape
    
    # # Linear Interpolation
    # interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
    # zi_lin = interp_lin(xi, -yi+6)
    # zi_linP = interp_lin(xi, yi)
    # # Compute the real and imaginary part of the wave with amplitude Amp
    # wave_real = np.cos( np.multiply(omega, zi_linP) )
    # wave_real = np.multiply(Amp, wave_real)
    # wave_imaginary = np.sin( np.multiply(omega, zi_linP) )
    # wave_imaginary = np.multiply(Amp, wave_imaginary)
    
    # # Figure with the level sets
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-18,18)
    # ax.set_ylim(-18, 24)
    # im2_5 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2, levels = 30)
    # plt.title("Linear interpolation, circle, eta1 = 1, eta2 = " +  str(indicesRef[i]) )
    # plt.show(block = False)
    # plt.colorbar(im2_5)
    nameFig = path_to_figures + '_LinearInt_n'+ str(i) +'.png'
    # plt.savefig(nameFig, dpi=my_dpi * 10)
    figuresLevelSet += [nameFig]
    
    # Real part of the wave
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-18,18)
    # ax.set_ylim(-18, 24)
    # im2_5 = plt.imshow(wave_real, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    # plt.title("Real part wave, circle, eta1 = 1, eta2 = " +  str(indicesRef[i]) )
    # plt.show(block = False)
    # plt.colorbar(im2_5)
    nameFig = path_to_figures + '_WaveReal_n'+ str(i) +'.png'
    # plt.savefig(nameFig, dpi=my_dpi * 10)
    figuresWaveReal += [nameFig]
    
    # Imaginary part of the wave
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-18,18)
    # ax.set_ylim(-18, 24)
    # im2_5 = plt.imshow(wave_imaginary, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    # plt.title("Real part wave, circle, eta1 = 1, eta2 = " +  str(indicesRef[i]) )
    # plt.show(block = False)
    # plt.colorbar(im2_5)
    nameFig = path_to_figures + '_WaveImaginary_n'+ str(i) +'.png'
    # plt.savefig(nameFig, dpi=my_dpi * 10)
    figuresWaveImaginary += [nameFig]
    
    # Figure with the gradients
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # ax = plt.gca()
    # ax.set_xlim(-18,18)
    # ax.set_ylim(-18, 24)
    # im2_13 = plt.imshow( zi_linP, cmap = colormap2, extent=[-18,18,-18,24], origin='lower'  )
    # plt.quiver(coords[indicesGradients, 0], coords[indicesGradients, 1], eik_grads[indicesGradients, 0], eik_grads[indicesGradients, 1])
    # plt.title("Linear interpolation and computed eikonal gradient, circle, eta1 = 1, eta2 = " +  str(indicesRef[i]) )
    # plt.show(block = False)
    # plt.colorbar(im2_13)
    nameFig = path_to_figures + '_Gradients_n'+ str(i) +'.png'
    # plt.savefig(nameFig, dpi=my_dpi * 10)
    figuresGradients += [nameFig]

######    
# Then we read those images and put them together in a video
im = cv2.imread(figuresGradients[0])
width, height, layers = im.shape
width = round(width/4)
height = round(height/4)

# For the level sets
fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
video=cv2.VideoWriter(path_to_figures + 'LevelSets_videoS.mp4', fourcc, fps,(width,height))
for j in range(0,nIndicesRef):
    im = cv2.imread(figuresLevelSet[j])
    im = cv2.resize(im, (width, height))
    video.write(im)
cv2.destroyAllWindows()
video.release()

# For the real part of the wave
fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
video=cv2.VideoWriter(path_to_figures + 'WaveReal_videoS.mp4', fourcc, fps,(width,height)) 
for j in range(0,nIndicesRef):
    im = cv2.imread(figuresWaveReal[j])
    im = cv2.resize(im, (width, height))
    video.write(im)
cv2.destroyAllWindows()
video.release()

# For the imaginary part of the wave
fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
video=cv2.VideoWriter(path_to_figures + 'WaveImaginary_videoS.mp4', fourcc, fps,(width,height))
for j in range(0,nIndicesRef):
    im = cv2.imread(figuresWaveImaginary[j])
    im = cv2.resize(im, (width, height))
    video.write(im)
cv2.destroyAllWindows()
video.release()

# For the level sets
fourcc = cv2.VideoWriter_fourcc(*'mp4v') 
video=cv2.VideoWriter(path_to_figures + 'Gradients_videoS.mp4', fourcc, fps,(width,height))
for j in range(0,nIndicesRef):
    im = cv2.imread(figuresGradients[j])
    im = cv2.resize(im, (width, height))
    video.write(im)
cv2.destroyAllWindows()
video.release()
    
        