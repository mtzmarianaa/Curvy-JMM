# Generation of test geometry and producing its triangular mesh
# I also generated a triangular mesh in a square because I'm going to 
# first test my fmm method here (very naivly)
import matplotlib.pyplot as plt
import numpy as np
from math import pi, sqrt
import meshpy.triangle as triangle
from matplotlib.patches import Arc
import numpy.linalg as la

############# Generate the two circles

def circ_base(t):
    '''
    Parametrization of the base cirlce (base of the snowman) for the test geometry.
    t can be a number or an array or numbers. If t is a number the output are the
    coordinates [x, y]. If t is an array the output is a n, 2 array, the first 
    column corresponds to the coordinates of x, the second to the coordinates of y.
    '''
    points = np.zeros((len(t), 2))
    points[:, 0] = 10*np.cos(t)
    points[:, 1] = 10*np.sin(t)
    return points

def circ_top(s):
    '''
    Parametrization of the top circle (head of the snowman) for the test geometry.
    '''
    points = np.zeros((len(s), 2))
    points[:, 0] = 5*sqrt(2)*np.cos(s)
    points[:, 1] = 5*sqrt(2) + 5*sqrt(2)*np.sin(s)
    return points

Nbase = 2**6
Ntop = 2**5

t = np.linspace(1/4, 9/4, num = Nbase , endpoint = False)
s = np.linspace(0, 1, num = Ntop )
t = np.multiply(t, pi)
s = np.multiply(s, pi)
t_arc1 = np.where( t >= pi/4 )
t_arc2 = np.where(t <= 3*pi/4)
t_arc = np.intersect1d(t_arc1, t_arc2)

points_base = circ_base(t)
points_top = circ_top(s)

plt.figure(1)
plt.plot( points_base[:, 0], points_base[:, 1], '-.', c='#6800ff', alpha = 0.3 )
plt.plot(points_top[:, 0], points_top[:, 1], '-.', c='#00f3ff', alpha = 0.3)
plt.show(block=False)



############# Generate the triangular mesh

def round_trip_connect(start, end):
    return [(i, i + 1) for i in range(start, end)] + [(end, start)]

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]

def needs_refinement(vertices, area):
        bary = np.sum(np.array(vertices), axis=0) / 3
        max_area = 1 + (la.norm(bary, np.inf) - 1) * 0.1
        return bool(area > max_area)

points_arch = points_base[t_arc, :  ]
points_nonArch = np.delete(points_base, t_arc, axis=0)

# The points on the outline of the snowman
points = [ (points_nonArch[i, 0], points_nonArch[i, 1]) for i in range( len(points_nonArch) ) ]
points.extend([ (points_top[i, 0], points_top[i, 1]) for i in range( Ntop ) ])
# Add the facets
facets = round_trip_connect(0, len(points) -1 ) 

# Then we need to add the arch
points.extend([ (points_arch[i, 0], points_arch[i, 1]) for i in range(1, len(points_arch)-1 ) ])




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
MaxN = 0
# Create the list of lists
for p in range(N_points):
    list_p = []
    for t in range(len(mesh_tris)):
        list_p += [point for point in mesh_tris[t, :] if p in mesh_tris[t, :] and point != p and point not in list_p]
    mesh_neigh.append( list_p )
    if MaxN < len(list_p): # to pad the lists so that we can save them as a np.array and have them as a 2d array in C
        MaxN = len(list_p)



# # #Plot
plt.figure(2)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
ax.add_patch(Arc((0,  5*sqrt(2)), (5*sqrt(2))*2,  (5*sqrt(2))*2, theta1=0.0, theta2=180, edgecolor="#000536", lw=1.5))
plt.title('Delaunay triangulation of test geometry')
plt.show(block=False)

r = 5*sqrt(2)
# # # # Zoom in 1
h1 = 1

plt.figure(3)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
ax.add_patch(Arc((0,  5*sqrt(2)), (5*sqrt(2))*2,  (5*sqrt(2))*2, theta1=0.0, theta2=180, edgecolor="#000536", lw=1.5))
plt.title('Delaunay triangulation of test geometry')
plt.xlim((r - h1, r+h1))
plt.ylim((r - h1, r+h1))
plt.show(block=False)

# # # # Zoom in 2
h2 = 0.25

plt.figure(4)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
ax.add_patch(Arc((0,  5*sqrt(2)), (5*sqrt(2))*2,  (5*sqrt(2))*2, theta1=0.0, theta2=180, edgecolor="#000536", lw=1.5))
plt.title('Delaunay triangulation of test geometry')
plt.xlim((r - h2, r+h2))
plt.ylim((r - h2, r+h2))
plt.show(block=False)

# # # # Zoom in 3
h3 = 0.01

plt.figure(5)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
ax.add_patch(Arc((0,  5*sqrt(2)), (5*sqrt(2))*2,  (5*sqrt(2))*2, theta1=0.0, theta2=180, edgecolor="#000536", lw=1.5))
plt.title('Delaunay triangulation of test geometry')
plt.xlim((r - h3, r+h3))
plt.ylim((r - h3, r+h3))
plt.show(block=False)


# # # # Zoom in 4
h4 = 1e-3

plt.figure(6)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
ax.add_patch(Arc((0,  5*sqrt(2)), (5*sqrt(2))*2,  (5*sqrt(2))*2, theta1=0.0, theta2=180, edgecolor="#000536", lw=1.5))
plt.title('Delaunay triangulation of test geometry')
plt.xlim((r - h4, r+h4))
plt.ylim((r - h4, r+h4))
plt.show(block=False)

# # # # Zoom in 5
h5 = 1e-5

plt.figure(7)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_points[:, 0], mesh_points[:, 1], mesh_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
ax.add_patch(Arc((0,  5*sqrt(2)), (5*sqrt(2))*2,  (5*sqrt(2))*2, theta1=0.0, theta2=180, edgecolor="#000536", lw=1.5))
plt.title('Delaunay triangulation of test geometry')
plt.xlim((r - h5, r+h5))
plt.ylim((r - h5, r+h5))
plt.show(block=False)


# # Now we save this triangulation to a bin file so that we can read it later from C

np.savetxt('MeshInfo/BoundaryPoints.txt', np.array(points), delimiter =', ' )

facets_arr = np.array(facets)
np.savetxt('MeshInfo/Facets.txt', facets_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('MeshInfo/MeshPoints.txt', mesh_points, delimiter =', ' )

np.savetxt('MeshInfo/Faces.txt', mesh_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('MeshInfo/NeighTriangles.txt', mesh_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

# Save the list of lists into a txt file
separator = "," 

with open("MeshInfo/Neigh.txt", "w") as out_file:
    for l in mesh_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)

##############################################################################
##############################################################################
##############################################################################
##############################################################################
###### SQUARE ON THE OUTSIDE

edges_square = [(-15, -15), (15, -15), (15, 18), (-15, 18)]

# Add them to the points
points_square = points + edges_square
# Add the facets
facets_square = facets + round_trip_connect(len(points), len(points_square)-1) 

# Add them to the triangulation (this is a new triangulation)
info_square = triangle.MeshInfo()
info_square.set_points(points_square)
info_square.set_facets(facets_square)

mesh_square = triangle.build(info_square,volume_constraints= True, refinement_func=needs_refinement)

mesh_square_points = np.array(mesh_square.points) # These are the points we want to export
mesh_square_tris = np.array(mesh_square.elements) # These are thee faces we want to export
mesh_square_neigTriangles = np.array(mesh_square.neighbors) # These are the neighbors we want to export
# Now we want to create an np array where the rows are the #of the point, the columns each one of its neighbours
# This is the most naive way of creating such thing, might be useful to optimize it later (?)
# We look in the row of the mesh.elements file
N_points_square = len(mesh_square_points)
mesh_square_neigh = []
MaxN_square = 0
# Create the list of lists
for p in range(N_points_square):
    list_p = []
    for t in range(len(mesh_square_tris)):
        list_p += [point for point in mesh_square_tris[t, :] if p in mesh_square_tris[t, :] and point != p and point not in list_p]
    mesh_square_neigh.append( list_p )
    if MaxN < len(list_p): # to pad the lists so that we can save them as a np.array and have them as a 2d array in C
        MaxN = len(list_p)




# # #Plot
plt.figure(8)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
ax.add_patch(Arc((0,  5*sqrt(2)), (5*sqrt(2))*2,  (5*sqrt(2))*2, theta1=0.0, theta2=180, edgecolor="#000536", lw=1.5))
plt.title('Delaunay triangulation of test geometry with rectangle')
plt.show(block=False)


## Save them as well

np.savetxt('MeshInfo/BoundaryPoints_Sq.txt', np.array(points_square), delimiter =', ' )

facets_square_arr = np.array(facets)
np.savetxt('MeshInfo/Facets_Sq.txt', facets_square_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('MeshInfo/MeshPoints_Sq.txt', mesh_square_points, delimiter =', ' )

np.savetxt('MeshInfo/Faces_Sq.txt', mesh_square_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('MeshInfo/NeighTriangles_Sq.txt', mesh_square_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

separator = "," 

with open("MeshInfo/Neigh_Sq.txt", "w") as out_file:
    for l in mesh_square_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)

plt.show()