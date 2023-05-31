
# Generation of test geometry and producing its triangular mesh
# I also generated a triangular mesh in a square because I'm going to 
# first test my fmm method here (very naivly)
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt, cos, sin, floor
import meshpy.triangle as triangle
import numpy.linalg as la
from matplotlib.colors import ListedColormap

############# Generate the two circles

def arch(t):
    '''
    Parametrization of the base cirlce (base of the snowman) for the test geometry.
    t can be a number or an array or numbers. If t is a number the output are the
    coordinates [x, y]. If t is an array the output is a n, 2 array, the first 
    column corresponds to the coordinates of x, the second to the coordinates of y.
    '''
    points = []
    for i in range(len(t)):
        points += [ (10*cos(t[i]), 10*sin(t[i]) -2  ) ]
    return points

Nbase = 2**4

t = np.linspace(0, 1, num = Nbase )
t = np.multiply(t, pi)

points_arch = arch(t)




############# Generate the triangular mesh

def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]

def needs_refinement(vertices, area):
        bary = np.sum(np.array(vertices), axis=0) / 3
        max_area = 1 + (la.norm(bary, np.inf) - 1) * 0.1
        return bool(area > max_area)


# The points on the outline of the snowman
points = [(-10, 10), (10, 10), (10, -2), (10, -10), (-10, -10), (-10, -2)]
# Add the facets
facets = round_trip_connect(0, len(points) -1 ) 

# Then we need to add the arch
points += points_arch[1:-1]
facets += connect(6, len(points)-1)

facets += [ (2,6), (5, len(points)-1) ]



# # # Set the information for the triangle mesh
info = triangle.MeshInfo()
info.set_points(points)
info.set_facets(facets)

# # # # Build the mesh

mesh = triangle.build(info,volume_constraints= True, refinement_func=needs_refinement)

mesh_points = np.array(mesh.points) # These are the points we want to export
mesh_tris = np.array(mesh.elements) # These are thee faces we want to export
mesh_neigTriangles = np.array(mesh.neighbors)
# Now we want to create an np array where the rows are the #of the point, the columns each one of its neighbours
# This is the most naive way of creating such thing, might be useful to optimize it later (?)
# We look in the row of the mesh.elements file
N_points = len(mesh_points)
mesh_neigh = []
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
    print("\npoints")
    p1 = mesh_points[mesh_tris[fi, 0] , :]
    print(p1)
    point1 = (p1[0], p1[1])
    n1 =  sqrt( p1[0]**2 + (p1[1]+2)**2  )
    print(n1)
    p2 = mesh_points[mesh_tris[fi, 1] , :]
    print(p2)
    n2 = sqrt( p2[0]**2 + (p2[1]+2)**2  )
    print(n2)
    point2 = (p2[0], p2[1])
    p3 = mesh_points[mesh_tris[fi, 2] , :]
    print(p3)
    n3 = sqrt( p3[0]**2 + (p3[1]+2)**2  )
    print(n3)
    point3 = (p3[0], p3[1])
    if ( n1 >= 10 and n2 >= 10 and n3 >= 10 and p1[1]>= 0 and p2[1] >= 0 and p3[1] >= 0 ):
        faces_label += [2]
        print("2")
        colors += [0.3]
    elif( (n1+n2+n3)/3 >= 10 and p1[1]>= 0 and p2[1]>= 0 and p3[1]>= 0 ):
        faces_label += [2]
        colors += [0.3]
    elif( p1[1] == 10 and p2[1] == 10 and p2[0] >= -2 and p1[0] >= -2 ):
        faces_label += [2]
        colors += [0.3]
    elif( p1[1] == 10 and p3[1] == 10 and p3[0] >= -2 and p1[0] >= -2 ):
        faces_label += [2]
        colors += [0.3]
    elif( p2[1] == 10 and p3[1] == 10 and p3[0] >= -2 and p2[0] >= -2 ):
        faces_label += [2]
        colors += [0.3]
    else:
        faces_label += [1]
        colors += [0.5]
        print("1")
  

# # #Plot
plt.figure(1)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
plt.title('Delaunay triangulation of test arch')
plt.show(block=False)

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
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.tripcolor(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, faces_label, cmap = newcmp)
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
plt.title('Delaunay triangulation of test arch')
plt.show(block=False)


plt.show()



# # # Now we save this triangulation to a bin file so that we can read it later from C

np.savetxt('TestArch/BoundaryPoints.txt', np.array(points), delimiter =', ', fmt = '%.8f' )

facets_arr = np.array(facets)
np.savetxt('TestArch/Facets.txt', facets_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestArch/MeshPoints.txt', mesh_points, delimiter =', ', fmt = '%.8f' )

np.savetxt('TestArch/Faces.txt', mesh_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestArch/NeighTriangles.txt', mesh_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

# Save the list of lists into a txt file
separator = "," 

with open("TestArch/Neigh.txt", "w") as out_file:
    for l in mesh_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestArch/IncidentFaces.txt", "w") as out_file:
    for l in mesh_IncidentFaces:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestArch/FacesLabel.txt", "w") as out_file:
    for l in faces_label:
        out_string = separator.join(str(l)) + "\n"
        out_file.write(out_string)