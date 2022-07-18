# Generation of super baby test geometry (just a square)
import matplotlib.pyplot as plt
import numpy as np
import meshpy.triangle as triangle
import numpy.linalg as la
from math import sqrt

def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]

# def needs_refinement(vertices, area):
#         bary = np.sum(np.array(vertices), axis=0) /3
#         max_area = 0.1 + (la.norm(bary, np.inf)-1) * 0.1
#         return bool(area > max_area)

# def needs_refinement(vertices, area ):
#     vert_origin, vert_destination, vert_apex = vertices
#     bary_x = (vert_origin.x + vert_destination.x + vert_apex.x) / 3
#     bary_y = (vert_origin.y + vert_destination.y + vert_apex.y) / 3
#     dist_center = sqrt( (bary_x)**2 + (bary_y)**2 )
#     max_area = abs( 0.05 * (dist_center-0.5) ) + 0.01
#     return area > max_area
    
def needs_refinement(vertices, area):
        max_area = 0.05
        return bool(area > max_area)

edges_square = [(-5, -5), (5, -5), (5, 5), (-5, 5)]

facets = [(0, 1), (1, 2), (2, 3), (3, 0)]

edges_square = edges_square + [(0, 0)]

#Set information for the mesh

info = triangle.MeshInfo()
info.set_points(edges_square)
info.set_facets(facets)

mesh = triangle.build(info,volume_constraints= True, refinement_func=needs_refinement)

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
for fi in range(len(mesh_tris)):
    faces_label += [1]

plt.figure(1)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
plt.title('Delaunay triangulation of test square')
plt.show()



# # Now we save this triangulation to a bin file so that we can read it later from C

np.savetxt('TestSquare/BoundaryPoints.txt', np.array(edges_square), delimiter =', ', fmt = '%.8f' )

facets_arr = np.array(facets)
np.savetxt('TestSquare/Facets.txt', facets_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestSquare/MeshPoints.txt', mesh_points, delimiter =', ', fmt = '%.8f' )

np.savetxt('TestSquare/Faces.txt', mesh_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestSquare/NeighTriangles.txt', mesh_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

# Save the list of lists into a txt file
separator = "," 

with open("TestSquare/Neigh.txt", "w") as out_file:
    for l in mesh_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestSquare/IncidentFaces.txt", "w") as out_file:
    for l in mesh_IncidentFaces:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestSquare/FacesLabel.txt", "w") as out_file:
    for l in faces_label:
        out_string = separator.join(str(l)) + "\n"
        out_file.write(out_string)