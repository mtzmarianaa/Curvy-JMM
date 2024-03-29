#pragma once

#include "neighbors.h"
#include "files_methods.h"

#include <stdlib.h>

typedef struct {
  // Boundary curve struct in 2D
  size_t i_Edge; // to which edge these tangents belong to
  double B[2][2]; // B1 and B2, tangents at the end point of the edge i
} boundaryCurve;

typedef struct mesh2S{
  // Mesh struct in 2D that also captures information
  // of the boundary (which edges form part of a
  // boundary and the boundary tangent at those points)
  double (*points)[2]; // points on the mesh ALL OF THEM, not only those on the boundaries
  size_t nPoints; // number of points on the mesh
  size_t (*faces)[3]; // indices of the points who form faces SORTED
  size_t nFaces; // number of faces in the mesh
  size_t (*edges)[2]; // pairs of indices of points who form an edge on the mesh (not necesarily on the boundary) SORTED
  size_t nEdges; // number of edges in the mesh
  size_t (*edgesInFace)[3]; // the 3 edge indices that define each face
  neighborsRS *neighbors; // for each point i, its neighbors (i.e. there is a face that includes both points)
  neighborsRS *incidentFaces; // for each point i, its incident faces
  double *eta; // indices of refraction on each face
  boundaryCurve *h_i; // boundaryCurve struct for each edge
} mesh2S;

typedef struct triangleFan {
  size_t nRegions;
  double x0[2];
  double x1[2];
  double xHat[2];
  int  *listFaces; // list of the triangle indices in this triangle fan FILLED WITH -1 IF REDUCED TRIANGLE
  double *listIndices; // in a triangle fan the etas are going to be called indices
  int *listEdges; // edges in this triangle fan FILLED WITH -1 IF REDUCED TRIANGLE length 2*nRegions + 1
  double (*listxk)[2]; // list of points x0, x1, ..., xHat in the triangle fan
  double (*listB0k)[2]; // length nRegions + 1
  double (*listBk)[2]; // length nRegions + 1
  double (*listBkBk1)[2]; // length 2*nRegions
} triangleFanS;


void mesh2_alloc(mesh2S **mesh2);

void mesh2_dealloc(mesh2S **mesh2);

void boundaryCurve_alloc(boundaryCurve **h_i);

void boundaryCurve_dealloc(mesh2S **h_i);

void triangleFan_alloc(triangleFanS **triFan);

void triangleFan_dalloc(triangleFanS **triFan);

void boundaryCurve_init(boundaryCurve *h_i, size_t i_Edge, double B[2][2]);

void mesh2_init(mesh2S *mesh2, double (*points)[2], size_t nPoints,
		size_t (*faces)[3], size_t nFaces, size_t (*edges)[2],
		size_t (*edgesInFace)[3], size_t nEdges, neighborsRS *neighbors,
		neighborsRS *incidentFaces, double *eta, boundaryCurve *h_i) ;

void boundaryCurve_init_from_meshpy(boundaryCurve *h_i, size_t nEdges, char const *pathBoundary);

void mesh2_init_from_meshpy(mesh2S *mesh2, char const *pathPoints, char const *pathFaces,
			    char const *pathEdges, char const *pathEdgesInFace,
			    char const *pathNeighbors, char const *pathIncidentFaces,
			    char const *pathIndices, char const *pathBoundary) ;

void printGeneralInfoMesh(mesh2S *mesh2);

void printEverythingInMesh(mesh2S *mesh2);

void twoTrianglesFromEdge(mesh2S *mesh2, size_t index0, size_t index1,
			  int possibleTriangles[2], size_t possibleThirdVertices[2] );

double minEtaFromTwoPoints(mesh2S *mesh2, size_t index0, size_t index1);

int faceBetween3Points(mesh2S *mesh2, size_t index0, size_t index1, size_t index2);


void triangleFan_initFromIndices(triangleFanS *triFan, mesh2S *mesh2, size_t nRegions,
				 size_t index0, size_t index1,
				 size_t indexHat, size_t *listIndicesNodes);


void triangleFan_init(triangleFanS *triFan, size_t nRegions, double x0[2],
		      double x1[2], double xHat[2],
		      int *listFaces, double *listIndices,
		      int *listEdges, double (*listxk)[2],
		      double (*listB0k)[2], double (*listBk)[2],
		      double (*listBkBk1)[2]);

size_t allSameTriangles(triangleFanS *triFan);

double thetaSnells(double grad[2], double normal[2], double eta1, double eta2);

void gradFromSnells(mesh2S *mesh2, triangleFanS *triFan, size_t index0, size_t index1,
		    double grad0[2], double grad1[2],
		    double grad0Snell[2], double grad1Snell[2])  ;

void oneGradFromSnells(mesh2S *mesh2, double point[2], double xHat[2],
		       double gradOutside[2],
		       double tanChange[2], double etaOutside, double etaInside,
		       double gradSnell[2]) ;

void getTangentChangeReg(mesh2S *mesh2, size_t indexPoint, size_t firstTriangle,
			 double tanChange[2], double etaOutside);

size_t pointOnBoundary(mesh2S *mesh2, size_t indexPoint);

void printEverythingTriFan(triangleFanS *triFan);



