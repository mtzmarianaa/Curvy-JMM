# Generate test geometry but just its base (i.e. just the circle) specifying h# Generation of test geometry and producing its triangular mesh
# I also generated a triangular mesh in a square because I'm going to 
# first test my fmm method here (very naivly)
import matplotlib.pyplot as plt
import numpy as np
from numpy.linalg import norm
from math import pi, sqrt, floor, ceil
import meshpy.triangle as triangle
from matplotlib.patches import Arc
import numpy.linalg as la
from matplotlib.colors import ListedColormap
import colorcet as cc



colormap2 = "cet_linear_worb_100_25_c53_r"
colormap2_r = "cet_linear_worb_100_25_c53"




plt.ion()

h = 0.0075
h_string = str(h)
eta1 = 1.0
eta2 = 1.45
currentH = "H12"
path_figures = "/Users/marianamartinez/Documents/Documents - Mariana’s MacBook Pro/NYU-Courant/FMM-bib/Figures/TestBaseSnow/"
path_info = '/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/'

############# Generate the two circles

my_dpi=96

def circ_base(t):
    '''
    Parametrization of the base cirlce (base of the snowman) for the test geometry.
    t can be a number or an array or numbers. If t is a number the output are the
    coordinates [x, y]. If t is an array the output is a n, 2 array, the first 
    column corresponds to the coordinates of x, the second to the coordinates of y.
    '''
    points = np.zeros((len(t), 2))
    tan_points = np.zeros((len(t), 2))
    points[:, 0] = 10*np.cos(t)
    points[:, 1] = 10*np.sin(t)
    tan_points[:, 0] = -np.sin(t)
    tan_points[:, 1] = np.cos(t)
    return points, tan_points



############# Useful functions

def round_trip_connect(start, end, boundary_tan = False, tans = None, hc = None):
    facets = [(i, i + 1) for i in range(start, end)] + [(end, start)]
    if( boundary_tan ):
        # This means that we also want to get the boundary tan struct
        B = []
        for i in range(start, end):
            tans_i = tans[i-1]/norm(tans[i-1])*hc[i-1]
            tans_i1 = tans[i]/norm(tans[i])*hc[i]
            B.append( (i-1, tans_i[0], tans_i[1], tans_i1[0], tans_i1[1] ) )
        tans_n = tans[end-1]/norm(tans[end-1])*hc[end -1]
        tans1 = tans[start-1]/norm(tans[start-1])*hc[start -1]
        B.append( (end -1, tans_n[0], tans_n[1], tans1[0], tans1[1] ) )
        return facets, B
    else:
        return facets

def connect(start, end):
    return [(i, i + 1) for i in range(start, end)]

def needs_refinement(vertices, area):
    max_area = h
    return bool(area > max_area)
    
def nPoints(h):
    '''
    With this function we know how many points we need in the top of the snowman and in the bottom, depending on which h we want
    WE ASSUME THAT WE ARE DISTRIBUTING UNIFORMLY THESE POINTS
    '''
    nBottom = ceil(15*pi/h) 
    return nBottom

Nbase = nPoints(h)

t = np.linspace(3/4, 11/4, num = Nbase , endpoint = False)
t = np.multiply(t, pi)

points_base, tan_base = circ_base(t)
hc = [norm( points_base[0, :] - points_base[1, :] )]*len(points_base) # All are the same because the points are equispaced on the circle
    

# The points on the outline of the snowman and the tangent
points =  [(-15, -10)] # start with x0 so that in each mesh we know its index
points.extend([ (points_base[i, 0], points_base[i, 1]) for i in range( len(points_base) ) ] )
# Add the facets
facets, mesh_square_boundaryCurve = round_trip_connect(1, len(points) -1, boundary_tan = True, tans = tan_base, hc = hc ) # we don't want to connnect x0



##############################################################################
##############################################################################
##############################################################################
##############################################################################
###### SQUARE ON THE OUTSIDE

edges_square = [(-18, -18), (18, -18), (18, 24), (-18, 24)]

# Add them to the points 
points_square = points + edges_square
# Add the facets
facets_square = facets + round_trip_connect(len(points), len(points_square)-1) 

# Add them to the triangulation (this is a new triangulation)
info_square = triangle.MeshInfo()
info_square.set_points(points_square)
info_square.set_facets(facets_square)

mesh_square = triangle.build(info_square,volume_constraints= True, refinement_func=needs_refinement, attributes=True, generate_faces=True)
mesh_square_points = np.array(mesh_square.points) # These are the points we want to export
mesh_square_tris = np.array(mesh_square.elements) # These are thee faces we want to export
mesh_square_tris = np.sort(mesh_square_tris, axis = 1) # But we want them in order
mesh_square_neigTriangles = np.array(mesh_square.neighbors) # These are the neighbors we want to export
mesh_square_neigTriangles = np.sort(mesh_square_neigTriangles, axis = 1) # We also want these sorted
mesh_square_facets = np.array(mesh_square.facets) # These are the edges we want to export
mesh_square_boundaryCurve = np.array(mesh_square_boundaryCurve) # Because our C code adds the NULL automatically

print("Triangulation done, number of points: ", len(mesh_square_points))

# Now we want to create an np array where the rows are the #of the point, the columns each one of its neighbours
# This is the most naive way of creating such thing, might be useful to optimize it later (?)
# We look in the row of the mesh.elements file
N_points_square = len(mesh_square_points)
mesh_square_neigh = [[] for i in range(N_points_square)]
MaxN_square = 0


# Create the list of lists

for t in range(len(mesh_square_tris)):
    # We iterate through each triangle
    p0 = mesh_square_tris[t,0]
    p1 = mesh_square_tris[t,1]
    p2 = mesh_square_tris[t,2]
    if p0 not in mesh_square_neigh[p1]:
        mesh_square_neigh[p1].append(p0)
    if p0 not in mesh_square_neigh[p2]:
        mesh_square_neigh[p2].append(p0)
    if p1 not in mesh_square_neigh[p0]:
        mesh_square_neigh[p0].append(p1)
    if p1 not in mesh_square_neigh[p2]:
        mesh_square_neigh[p2].append(p1)
    if p2 not in mesh_square_neigh[p0]:
        mesh_square_neigh[p0].append(p2)
    if p2 not in mesh_square_neigh[p1]:
        mesh_square_neigh[p1].append(p2)

print("List of neighbors done")


# Create the list for the incident faces and edges

mesh_IncidentFaces_sq = [[] for i in range(N_points_square)]
facets_created = []

for t in range(len(mesh_square_tris)):
    # Iterate through each triangle
    p0 = mesh_square_tris[t, 0]
    p1 = mesh_square_tris[t, 1]
    p2 = mesh_square_tris[t, 2]
    # For the incident faces
    mesh_IncidentFaces_sq[p0].append(t)
    mesh_IncidentFaces_sq[p1].append(t)
    mesh_IncidentFaces_sq[p2].append(t)
    # For the edges
    facets_created.append((p0, p1))
    facets_created.append((p0, p2))
    facets_created.append((p1,p2))

# Combine the facets_created with mesh_square_facets
facets_created = np.array(facets_created)
mesh_square_facets = np.concatenate((mesh_square_facets, facets_created), axis=0)
indices = np.unique(mesh_square_facets, return_index=True, axis=0)[1] # Index of the unique rows
indices = np.sort(indices)
mesh_square_facets = mesh_square_facets[indices]

print("List of incident faces done")
    
# Create the labels per face


faces_label = []
colors = []
r = 5*sqrt(2)
for fi in range(len(mesh_square_tris)):
    # we need to get the coordinates of the 3 points that define that upper triangle
    p1 =  mesh_square_points[mesh_square_tris[fi, 0] , :]
    p2 =  mesh_square_points[mesh_square_tris[fi, 1] , :]
    p3 =  mesh_square_points[mesh_square_tris[fi, 2] , :]
    if ( sqrt(p1[0]**2 + p1[1]**2)<= 10 ) and ( sqrt(p2[0]**2 + p2[1]**2)<= 10 ) and (sqrt(p3[0]**2 + p3[1]**2)<= 10):
        faces_label += [eta2]
        colors += [0.3]
    elif( (sqrt(p1[0]**2 + p1[1]**2) + sqrt(p2[0]**2 + p2[1]**2) + sqrt(p3[0]**2 + p3[1]**2))<= 30  ):
        faces_label += [eta2]
        colors += [0.3]
    else:
        faces_label += [eta1]
        colors += [0.1]

print("Lists of faces labels done")




# #Plot
plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
ax = fig.gca()
plt.gca().set_aspect('equal')
plt.triplot(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, '-.', lw=0.5, c='#6800ff')
circle_b = plt.Circle((0, 0), 10, color="#000536",fill=False)
ax.add_patch(circle_b)
plt.title('Delaunay triangulation of test geometry, H = '+h_string)
plt.savefig(path_figures + currentH + '/'+ currentH + '_TriangulationWhite.png', dpi=my_dpi * 10)
plt.xlim([-18, 18])
plt.ylim([-18, 24])
plt.show(block=False)


viridis = plt.get_cmap(colormap2, 256)
new_colors = viridis(np.linspace(0, 1, 256))
blue = new_colors[0]
purple = new_colors[200]
mid = floor(len(new_colors)/2)
new_colors = [blue]*mid
new_colors += [purple]*mid
newcmp = ListedColormap(new_colors)


plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
plt.gca().set_aspect('equal')
plt.tripcolor(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, faces_label, cmap = colormap2)
plt.triplot(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, '-.', lw=0.5, c='#00fffb')
plt.title('Delaunay triangulation of test geometry, H = '+h_string)
plt.xlim([-18, 18])
plt.ylim([-18, 24])
plt.show(block = False)
plt.savefig(path_figures + currentH + '/'+ currentH + '_Triangulation.png', dpi=my_dpi * 10)


# With the tangents
plt.figure(figsize=(800/my_dpi, 800/my_dpi), dpi=my_dpi)
fig = plt.gcf()
plt.gca().set_aspect('equal')
plt.tripcolor(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, faces_label, cmap = colormap2)
plt.triplot(mesh_square_points[:, 0], mesh_square_points[:, 1], mesh_square_tris, '-.', lw=0.5, c='#00fffb')
plt.quiver(points_base[:,0], points_base[:,1], mesh_square_boundaryCurve[:, 1],
           mesh_square_boundaryCurve[:, 2], width = 0.005, scale = 1.0, scale_units = 'xy')
plt.title('Delaunay triangulation of test geometry, H = '+h_string)
plt.xlim([-18, 18])
plt.ylim([-18, 24])
plt.show(block = False)
plt.savefig(path_figures + currentH + '/'+ currentH + '_TriangulationTangents.png', dpi=my_dpi * 10)


## Save them as well

np.savetxt(path_info + currentH + '/'+ currentH +'_Edges.txt', mesh_square_facets.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_MeshPoints.txt', mesh_square_points, delimiter =', ', fmt = '%.8f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_Faces.txt', mesh_square_tris.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_NeighTriangles.txt', mesh_square_neigTriangles.astype(int), delimiter =', ', fmt ='%.0f')

np.savetxt(path_info + currentH + '/'+ currentH +'_BoundaryCurve.txt', mesh_square_boundaryCurve, delimiter =', ', fmt = '%.8f' )

np.savetxt(path_info + currentH + '/'+ currentH +'_Indices.txt', np.array(faces_label), delimiter =', ', fmt = '%.8f' )

separator = "," 

with open(path_info + currentH + '/'+ currentH + "_Neigh.txt", "w") as out_file:
    for l in mesh_square_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open(path_info + currentH + '/'+ currentH + "_IncidentFaces.txt", "w") as out_file:
    for l in mesh_IncidentFaces_sq:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
plt.show()
