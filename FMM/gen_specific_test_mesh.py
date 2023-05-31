# Generate test for the artificial triangle update
import numpy as np

meshPoints = [(1, -1), (0,0), (2,0), (2,1), (1,2), (0, 3)]
faces = [(0, 1, 2), (5, 1, 2), (2, 5, 4), (2, 3, 4)]
faces = np.array(faces)
facets = [(0, 1), (1, 5), (5, 4), (4, 3), (3, 2), (2, 0)]
N_points_square = len(meshPoints)
mesh_square_neigh = []
MaxN_square = 0
# Create the list of lists
for p in range(N_points_square):
    list_p = []
    for t in range(len(faces)):
        list_p += [point for point in faces[t, :] if p in faces[t, :] and point != p and point not in list_p]
    mesh_square_neigh.append( list_p )
    
mesh_IncidentFaces_sq = []
for p in range(N_points_square):
    list_faces = []
    for t in range(len(faces )):
        if (p in faces[t, :]):
            list_faces += [t]
    mesh_IncidentFaces_sq.append(list_faces)
    
faces_label = [1]*N_points_square


np.savetxt('TestBaseSnow/TestArtificial/BoundaryPoints.txt', np.array(meshPoints), delimiter =', ', fmt = '%.8f' )

facets_arr = np.array(facets)
np.savetxt('TestBaseSnow/TestArtificial/Facets.txt', facets_arr.astype(int), delimiter =', ', fmt ='%.0f' )

np.savetxt('TestBaseSnow/TestArtificial/MeshPoints.txt', meshPoints, delimiter =', ', fmt = '%.8f' )

np.savetxt('TestBaseSnow/TestArtificial/Faces.txt', faces.astype(int), delimiter =', ', fmt ='%.0f' )

separator = "," 

with open("TestBaseSnow/TestArtificial/Neigh.txt", "w") as out_file:
    for l in mesh_square_neigh:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestBaseSnow/TestArtificial/IncidentFaces.txt", "w") as out_file:
    for l in mesh_IncidentFaces_sq:
        out_string = separator.join(str(x) for x in l) + "\n"
        out_file.write(out_string)
        
with open("TestBaseSnow/TestArtificial/FacesLabel.txt", "w") as out_file:
    for l in faces_label:
        out_string = separator.join(str(l)) + "\n"
        out_file.write(out_string)