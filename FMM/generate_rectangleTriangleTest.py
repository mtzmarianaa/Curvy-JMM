# Generation of  test geometry (just a square with an inverted triangle)
import matplotlib.pyplot as plt
import numpy as np
import meshpy.triangle as triangle
import numpy.linalg as la
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from math import floor
import pandas as pd

my_dpi=96
h = 10
h_string = str(h)
# h_stringPrev = "H0_1"

# previousPoints = np.genfromtxt("/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/" + h_stringPrev + "/" + h_stringPrev + "_MeshPoints.txt", delimiter=",")
# previousPoints = previousPoints
# edges_square = [  (p[0], p[1]) for p in previousPoints  ]
# print("Amount of previous points:  ", len(edges_square))

def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]


def needs_refinement(vertices, area):
        max_area = h
        return bool(area > max_area)

edges_square = [ (0,-2) ,(-10, -10), (10, -10), (10, 10), (-10, 10)]

facets = round_trip_connect(1, 4  ) 

edges_square += [(0, 0)]

facets += [(3, 5), (5, 4)]


# We add the inverted triangle part

#Set information for the mesh

info = triangle.MeshInfo()
info.set_points(edges_square)
info.set_facets(facets)

mesh = triangle.build(info,volume_constraints= True, refinement_func=needs_refinement, attributes=True, generate_faces=True)

mesh_points = np.array(mesh.points) # These are the points we want to export
mesh_tris = np.array(mesh.elements) # These are the faces we want to export
mesh_neigTriangles = np.array(mesh.neighbors)

N_points = len(mesh_points)
mesh_neigh = []
MaxN = 0
# Create the list of lists
for p in range(N_points):
    list_p = []
    for t in range(len(mesh_tris)):
        list_p += [point for point in mesh_tris[t, :] if p in mesh_tris[t, :] and point != p and point not in list_p]
    mesh_neigh.append( list_p )
    
# Now we want an array which has a list of the indices of the faces that are incident on each vertex

mesh_IncidentFaces = []
for p in range(N_points):
    list_faces = []
    for t in range(len(mesh_tris)):
        if (p in mesh_tris[t, :]):
            list_faces += [t]
    mesh_IncidentFaces.append(list_faces)
    
# We add the label for the indicator function to see in which region each face is in

faces_label = []
colors = []
for fi in range(len(mesh_tris)):
    # we need to get the coordinates of the 3 points that define that upper triangle
    p1 = mesh_points[mesh_tris[fi, 0] , :]
    p2 = mesh_points[mesh_tris[fi, 1] , :]
    p3 = mesh_points[mesh_tris[fi, 2] , :]
    if (p1[1]>= abs( p1[0] ) ) and ( p2[1]>= abs( p2[0]) ) and (p3[1]>= abs(p3[0])):
        faces_label += [2]
        colors += [0.3]
    else:
        faces_label += [1]
        colors += [0.5]
        

plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
plt.title('Delaunay triangulation of test square')
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/H0_1/H0_1_TriangulationWhite.png', dpi=my_dpi * 10)
plt.show(block = False)

viridis = plt.get_cmap('magma', 256)
new_colors = viridis(np.linspace(0, 1, 256))
blue = new_colors[15]
purple = new_colors[45]
mid = floor(len(new_colors)/2)
new_colors = [blue]*mid
new_colors += [purple]*mid
newcmp = ListedColormap(new_colors)


plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
plt.gca().set_aspect('equal')
plt.tripcolor(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, colors, cmap = newcmp)
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#00fffb')
plt.title('Delaunay triangulation of test square, H0_1, h=' + h_string )
plt.savefig('/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquareTriangle/H0_1/H0_1_Triangulation.png', dpi=my_dpi * 10)
plt.show(block = False)


plt.show()




# # # Now we save this triangulation to a bin file so that we can read it later from C

np.savetxt('TestTriangleSquare/H0_1/H0_1_BoundaryPoints.txt', np.array(edges_square), delimiter =', ', fmt = '%.8f' )

facets_arr = np.array(facets)
np.savetxt('TestTriangleSquare/H0_1/H0_1_Facets.txt', facets_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestTriangleSquare/H0_1/H0_1_MeshPoints.txt', mesh_points, delimiter =', ', fmt = '%.8f' )
print("Amount of new points:  ", len(mesh_points))

np.savetxt('TestTriangleSquare/H0_1/H0_1_Faces.txt', mesh_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestTriangleSquare/H0_1/H0_1_NeighTriangles.txt', mesh_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

# Save the list of lists into a txt file
separator = "," 

with open("TestTriangleSquare/H0_1/H0_1_Neigh.txt", "w") as out_file:
    for l in mesh_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestTriangleSquare/H0_1/H0_1_IncidentFaces.txt", "w") as out_file:
    for l in mesh_IncidentFaces:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestTriangleSquare/H0_1/H0_1_FacesLabel.txt", "w") as out_file:
    for l in faces_label:
        out_string = separator.join(str(l)) + "\n"
        out_file.write(out_string)