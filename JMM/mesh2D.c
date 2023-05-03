
/* TRIANGLE 2D MESH STRUCTURE
This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given
*/

#include "mesh2D.h"
#include "linAlg.h"
#include "SoSFunction.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


void mesh2_alloc(mesh2S **mesh2) {
  *mesh2 = malloc(sizeof(mesh2S));
  assert(*mesh2 != NULL);
}

void mesh2_dealloc(mesh2S **mesh2) {
  free(*mesh2);
  *mesh2 = NULL;
}

void boundaryCurve_alloc(boundaryCurve **h_i) {
  *h_i = malloc(sizeof(boundaryCurve));
  assert(*h_i != NULL);
}

void boundaryCurve_dealloc(mesh2S **h_i) {
  free(*h_i);
  *h_i = NULL;
}

void boundaryCurve_init(boundaryCurve *h_i, size_t i_Edge, double B[2][2]) {
  h_i->i_Edge = i_Edge;
  h_i->B = B;
}

void mesh2_init(mesh2S *mesh2, double (*points)[2], size_t nPoints, size_t (*faces)[3], size_t nFaces, size_t (*edges)[2], size_t nEdges, neighborsRS *neighbors, neighborsRS *incidentFaces, double *eta, boundaryCurve *h_i) {
  mesh2->points = points;
  mesh2->nPoints = nPoints;
  mesh2->faces = faces;
  mesh2->nFaces = nFaces;
  mesh2->edges = edges;
  mesh2->nEdges = nEdges;
  mesh2->neighbors = neighbors;
  mesh2->indicentFaces = incidentFaces;
  mesh2->eta = eta;
  mesh2->h_i = h_i;
}

void mesh2_init_from_meshpy(mesh2D *mesh2, char const *pathPoints, char const *pathFaces, char const *pathEdges, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathIndices, char const *pathBoundary) {
}


