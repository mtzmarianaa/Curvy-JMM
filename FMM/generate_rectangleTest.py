# Generation of super baby test geometry (just a square)
import matplotlib.pyplot as plt
import numpy as np
import meshpy.triangle as triangle
import numpy.linalg as la
from math import sqrt

h = 0.01
H_name = "H5"
my_dpi = 96
saveFigure = False

def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]

def needs_refinement(vertices, area):
        max_h = h
        return bool(area > max_h)

edges_square = [(0,0), (-5, -5), (5, -5), (5, 5), (-5, 5)]

facets = [(1, 2), (2, 3), (3, 4), (4, 1)]


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
plt.title('Delaunay triangulation of test square ' + H_name)
if(saveFigure):
    plt.savefig("/Users/marianamartinez/Documents/NYU-Courant/FMM-bib/Figures/TestSquare/" + H_name + "/" + H_name + "_Triangulation.png", dpi = my_dpi)
plt.show()



# # Now we save this triangulation to a bin file so that we can read it later from C

np.savetxt('TestSquare/'+ H_name + '/' + H_name +'_BoundaryPoints.txt', np.array(edges_square), delimiter =', ', fmt = '%.8f' )

facets_arr = np.array(facets)
np.savetxt('TestSquare/' + H_name + '/' + H_name +'_Facets.txt', facets_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestSquare/' + H_name + '/' + H_name +'_MeshPoints.txt', mesh_points, delimiter =', ', fmt = '%.8f' )

np.savetxt('TestSquare/' + H_name + '/' + H_name +'_Faces.txt', mesh_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestSquare/' + H_name + '/' + H_name +'_NeighTriangles.txt', mesh_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

# Save the list of lists into a txt file
separator = "," 

with open('TestSquare/' + H_name + '/' + H_name +'_Neigh.txt', "w") as out_file:
    for l in mesh_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open('TestSquare/' + H_name + '/' + H_name +'_IncidentFaces.txt', "w") as out_file:
    for l in mesh_IncidentFaces:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open('TestSquare/' + H_name + '/' + H_name +'_FacesLabel.txt', "w") as out_file:
    for l in faces_label:
        out_string = separator.join(str(l)) + "\n"
        out_file.write(out_string)
