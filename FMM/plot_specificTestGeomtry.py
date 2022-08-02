# SCRIPT TO VISUALIZE ERRORS (can I say this?)
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri
import matplotlib.animation as animation
import imageio
import os

colormap2 = plt.cm.get_cmap('magma')
sm2 = plt.cm.ScalarMappable(cmap=colormap2)

nx = 36*200
ny = 42*200
my_dpi=96


def rotate(angle):
    ax.view_init(azim=angle)
    
    
def buildPath(indexStart, paths):
    path = [paths[indexStart]] #how it starts
    previousInPath = paths[indexStart]
    while(previousInPath != paths[previousInPath]):
        previousInPath = paths[previousInPath]
        path.extend(  [previousInPath]  )
    return path

######################################################
######################################################
######################################################
#### 

# Global variables
eik_coords = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/TestIndex/MeshPoints.txt", delimiter=",")
triangles = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/TestIndex/Faces.txt", delimiter=",")

xi, yi = np.meshgrid(np.linspace(-18, 18, nx), np.linspace(-18, 24, ny))
# We need a triangulation object thing
triang = tri.Triangulation(eik_coords[:, 0], eik_coords[:, 1], triangles)

figsContours = []
figsLinInt = []

factors = [1/1.452, 1/1.348, 1/1.24, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4, 4.25, 4.5, 4.75, 5]

for i in range(1, 21):
    # eik_vals = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/TestIndex/ComputedValues_"+ str(i) + ".bin")
    # eik_paths = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/TestIndex/Paths_20"+ ".bin", dtype=np.uint32)
    # eik_grads = np.fromfile("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/TestIndex/ComputedGradients_"+ str(i) + ".bin");
    # eik_grads = eik_grads.reshape(len(eik_coords), 2)
    # # path1 = buildPath(1, eik_paths)
    # # path2 = buildPath(2, eik_paths)
    # # path3 = buildPath(3, eik_paths)
    # # path4 = buildPath(4, eik_paths)
    # nu1 = 1.348*factors[i-1]
    # nu2 = 1.452*factors[i-1]
    
    # # We plot the computed eikonal values in 3D
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # ax = plt.axes(projection='3d')
    # ax.scatter(eik_coords[:, 0], eik_coords[:, 1], eik_vals, c= eik_vals, cmap=colormap2)
    # plt.title("Computed eikonal values, test geometry nu1="+ str(nu1) + ', nu2='+str(nu2))
    # plt.show(block = False)
    # rot_animation = animation.FuncAnimation(fig, rotate, frames=np.arange(0, 362, 2), interval=100)
    # rot_animation.save('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/TestIndex/ComputedValues_'+str(i)+'.gif', dpi=80, writer='Pillow')
    
    # # We interpolate the solution 
    # # To be able to use LinearTriInterpolator
    # interp_lin = tri.LinearTriInterpolator(triang, eik_vals)
    # zi_lin = interp_lin(xi, -yi+6)
    
    # #Now we can plot according to the computed eikonal values
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # im_bar1 = plt.contourf(xi, 6-yi, zi_lin, cmap = colormap2, levels = 25)
    # plt.scatter(eik_coords[:, 0], eik_coords[:, 1], marker = '.' , c = eik_vals, cmap = colormap2)
    # plt.title("Linear interpolation, test geometry nu1="+ str(nu1) + ', nu2='+str(nu2))
    # plt.show(block = False)
    # plt.colorbar(im_bar1)
    figName_Contour = '/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/TestIndex/LinearInt_Contour_'+str(i)+'.png'
    figsContours.extend([figName_Contour])
    # plt.savefig(figName_Contour, dpi=my_dpi * 10)
    
    # # We plot it with imshow
    # fig = plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
    # plt.axis('equal')
    # im_bar2 = plt.imshow( zi_lin, cmap = colormap2, extent=[-18,18,-18,24]  )
    # # plt.plot( eik_coords[path1, 0], eik_coords[path1, 1], marker='.',linestyle='dashed', c = '#040072' )
    # # plt.plot( eik_coords[path2, 0], eik_coords[path2, 1], marker='.',linestyle='dashed', c = '#05009c' )
    # # plt.plot( eik_coords[path3, 0], eik_coords[path3, 1], marker='.',linestyle='dashed', c = '#0600c6' )
    # # plt.plot( eik_coords[path4, 0], eik_coords[path4, 1], marker='.',linestyle='dashed', c = '#0700e5' )
    # plt.title("Linear interpolation, test geometry nu1="+ str(nu1) + ', nu2='+str(nu2))
    # plt.show(block = False)
    # plt.colorbar(im_bar2)
    figName_LinInt = '/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/TestIndex/LinearInt_'+str(i)+'.png'
    figsLinInt.extend([figName_LinInt])
    # plt.savefig(figName_LinInt, dpi=my_dpi * 10)


# Build GIF of the contours
print('creating gif\n')
images = []
for filename in figsContours:
    images.append(imageio.imread(filename))
imageio.mimsave('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/TestIndex/ContoursChange.gif', images)

    
# Build GIF of the imshows
print('creating gif\n')
images = []
for filename in figsLinInt:
    images.append(imageio.imread(filename))
imageio.mimsave('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/Mesh_generation/TestIndex/LinIntChange.gif', images)
        



#plt.show()