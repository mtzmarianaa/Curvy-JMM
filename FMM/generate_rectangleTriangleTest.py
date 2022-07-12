# Generation of  test geometry (just a square with an inverted triangle)
import matplotlib.pyplot as plt
import numpy as np
import meshpy.triangle as triangle
import numpy.linalg as la
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from math import floor

def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]

def needs_refinement(vertices, area):
        bary = np.sum(np.array(vertices), axis=0) /2
        max_area = 1 + (la.norm(bary, np.inf) - 1) * 0.1
        return bool(area > max_area)

edges_square = [(-5, -5), (5, -5), (5, 5), (-5, 5)]

facets = round_trip_connect(0, len(edges_square)-1  ) 

edges_square += [(0, 0)]

facets += [(2, 4), (4, 3)]

edges_square = edges_square + [(0, -2)]

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
        

plt.figure(1)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
plt.title('Delaunay triangulation of test square')
plt.show(block = False)

viridis = plt.get_cmap('magma', 256)
new_colors = viridis(np.linspace(0, 1, 256))
blue = new_colors[15]
purple = new_colors[45]
mid = floor(len(new_colors)/2)
new_colors = [blue]*mid
new_colors += [purple]*mid
newcmp = ListedColormap(new_colors)


plt.figure(2)
fig = plt.gcf()
plt.gca().set_aspect('equal')
plt.tripcolor(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, colors, cmap = newcmp)
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#00fffb')
plt.title('Delaunay triangulation of test square')
plt.show(block = False)


plt.show()




# # # Now we save this triangulation to a bin file so that we can read it later from C

np.savetxt('TestTriangleSquare/BoundaryPoints.txt', np.array(edges_square), delimiter =', ', fmt = '%.8f' )

facets_arr = np.array(facets)
np.savetxt('TestTriangleSquare/Facets.txt', facets_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestTriangleSquare/MeshPoints.txt', mesh_points, delimiter =', ', fmt = '%.8f' )

np.savetxt('TestTriangleSquare/Faces.txt', mesh_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestTriangleSquare/NeighTriangles.txt', mesh_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

# Save the list of lists into a txt file
separator = "," 

with open("TestTriangleSquare/Neigh.txt", "w") as out_file:
    for l in mesh_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestTriangleSquare/IncidentFaces.txt", "w") as out_file:
    for l in mesh_IncidentFaces:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestTriangleSquare/FacesLabel.txt", "w") as out_file:
    for l in faces_label:
        out_string = separator.join(str(l)) + "\n"
        out_file.write(out_string)