#pragma once

#include "coord.h"
#include "faces.h"
#include "neighbors.h"
#include "files_methods.h"

typedef struct {
    // ask if everything in here is usefull/necessary or if its too much
    // ask if this should be an opaque type (according to me it shouldn't)
    // inspiration: what we might need + what might be useful to plot the mesh using triplot in Python
  coordS *points; // these are ALL the coordinates + number of points in the mesh
  neighborsRS *neighbors; // for each point i, its neighbors (i.e. there is a face that includes both points)
  neighborsRS *incidentFaces; // for each point i, its incident faces
  facesS *faces; // these are the faces, the triangles in the mesh
  int nPoints;
  int *indexRegions; // this is the "indicator function" of the type of region that each FACE is in
} triMesh_2Ds;

void triMesh_2Dalloc(triMesh_2Ds **triM_2D);

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D);

void triMesh2_init(triMesh_2Ds *triM_2D, coordS *points, neighborsRS *neighbors, neighborsRS *incidentFaces, facesS *faces, int nPoints, int *indexRegions);
 
void triMesh2_init_from_meshpy(triMesh_2Ds *triM_2D, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathFaces, char const *pathIndexRegions);

void printGeneralInfoMesh(triMesh_2Ds *triM_2D);

void printEverythingInMesh(triMesh_2Ds *triM_2D);

int regionBetweenTwoPoints(triMesh_2Ds *triM_2D, int index_from, int index_to);

int faceBetween3Points(triMesh_2Ds *triM_2D, int index1, int index2, int index3);

void twoTrianglesFromEdge(triMesh_2Ds *triM_2D, int index0, int index1, int possibleTriangles[2], int possibleThirdVertices[2]);

void pointWhereRegionChanges(triMesh_2Ds *triM_2D, int x0_ind, int x1_ind, int xHat, int directionToStart, int *xChange_ind, double *angle_xHat, double *angle_xChange);