
/* TRIANGLE 2D MESH STRUCTURE
This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given
*/

#include "mesh2D.h"
#include "linAlg.h"
#include "SoSFunction.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

void boundaryCurve_init_from_meshpy(boundaryCurve *h_i, size_t nEdges, char const *pathBoundary) {
  // structure of the file from meshpy: first 
}

void mesh2_init_from_meshpy(mesh2D *mesh2, char const *pathPoints, char const *pathFaces, char const *pathEdges, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathIndices, char const *pathBoundary) {
  // loads information from files generated from python meshpy
  size_t nPoints, nFaces, nEdges, (*faces)[3], (*edges)[2];
  double *eta, (*points)[2];
  neighborsRS *neighbors, *incidentFaces;
  boundaryCurve *h_i;

  // set sizes of things inside the struct
  nPoints = numLinesInFile(pathPoints); // set number of points
  nFaces = numLinesInFile(pathFaces); // set number of faces
  nEdges = numLinesInFile(pathEdges); // set number of edges
  mesh2->nPoints = nPoints;
  mesh2->nFaces = nFaces;
  mesh2->nEdges = nEdges;

  // mallocs
  points = malloc(nPoints*2*sizeof(double));
  faces = malloc(nFaces*3*sizeof(size_t));
  edges = malloc(nEdges*2*sizeof(size_t));
  neighbors = malloc(nPoints*sizeof(neighborsRS));
  incidentFaces = malloc(nPoints*sizeof(neighborsRS));
  h_i = malloc(nEdges*sizeof(boundaryCurve));

  // Now that everything is alloc'ed and declared we can open the files

  // for the points
  read_n2Files_double(points[0], pathPoints);
  printf("Read the points\n");
  mesh2->points = points;

  // for the faces
  read_n3File_int(faces[0], pathFaces);
  printf("Read the faces\n");
  mesh2->faces = faces;

  // for the edges
  read_n2File_int(edges[0], pathEdges);
  printf("Read the edges\n");
  mesh2->edges = edges;

  // read the neighbors
  neighbors_init(neighbors, pathNeighbors, nPoints);
  printf("Read neighbors\n");
  mesh2->neighbors = neighbors;
  neighbors_init(incidentFaces, pathIncidentFaces, nPoints);
  printf("Read indicent faces to points\n");
  mesh2->incidentFaces = incidentFaces;

  // init the boundary curve struct
  

}


void printGeneralInfoMesh(mesh2S *mesh2) {
  printf("\n  GENERAL INFORMATION ABOUT THIS MESH \n\n");
  printf("    Number of points: %u\n", mesh2->nPoints);
  printf("    Number of faces: %u\n", mesh2->nFaces);
  printf("    Number of edges: %u\n", mesh2->nEdges);
}


