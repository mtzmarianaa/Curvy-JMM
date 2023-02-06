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
  double (*boundary_tan)[2]; // tangent to the boundary, vector of zeros if point not on boundary
  int (*boundary_chain)[2]; // if not -1 immediate neighbours to a point on the boundary
} triMesh_2Ds;

typedef struct{
  // information regarding the triangle fan in an update
  //(since xHat and x1 are not necessarily neighbors)
  int nChanges; // number of changes in the index of refraction
  int updatable; // if the update with x0, x1, xHat should be considered or not
  // notice that if a point in the fan is valid then this update should not be considered
  int xChange_ind; // the index where the edge x0xChange changes region
  double angle_xHat; // the angle of the triangle xHat x0 x1
  double angle_xChange; // the angle of the triangle x2 x0 x1
  double indexRef_01; // index of refraction in the triangle x2 x0 x1
  double indexRef_02; // index of refraction in the triangle xHat x0 x2
} infoTriangleFan;

void triMesh_2Dalloc(triMesh_2Ds **triM_2D);

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D);

void infoTriangleFan_alloc(infoTriangleFan **infoOut);

void infoTriangleFan_dalloc(infoTriangleFan **infoOut);

void triMesh2_init(triMesh_2Ds *triM_2D, coordS *points, neighborsRS *neighbors, neighborsRS *incidentFaces, facesS *faces, int nPoints, int *indexRegions, double (*boundary_tan)[2], int (*boundary_chain)[2]);
 
void triMesh2_init_from_meshpy(triMesh_2Ds *triM_2D, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathFaces, char const *pathIndexRegions, char const *pathBoundaryTan, char const *pathBoundaryChain);

void printGeneralInfoMesh(triMesh_2Ds *triM_2D);

void printEverythingInMesh(triMesh_2Ds *triM_2D);

double regionBetweenTwoPoints(triMesh_2Ds *triM_2D, int index_from, int index_to);

int faceBetween3Points(triMesh_2Ds *triM_2D, int index1, int index2, int index3);

void twoTrianglesFromEdge(triMesh_2Ds *triM_2D, int index0, int index1, int possibleTriangles[2], int possibleThirdVertices[2]);

void pointWhereRegionChanges(triMesh_2Ds *triM_2D, int *current_states, int x0_ind, int x1_ind, int xHat, int directionToStart, infoTriangleFan *infoOut);
