
/* TRIANGLE 2D MESH STRUCTURE
This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given
*/

#include "mesh2D.h"
#include "linAlg.h"


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
  h_i->B[0][0] = B[0][0];
  h_i->B[0][1] = B[0][1];
  h_i->B[1][0] = B[1][0];
  h_i->B[1][1] = B[1][1];
}

void mesh2_init(mesh2S *mesh2, double (*points)[2], size_t nPoints, size_t (*faces)[3], size_t nFaces, size_t (*edges)[2], size_t nEdges, neighborsRS *neighbors, neighborsRS *incidentFaces, double *eta, boundaryCurve *h_i) {
  mesh2->points = points;
  mesh2->nPoints = nPoints;
  mesh2->faces = faces;
  mesh2->nFaces = nFaces;
  mesh2->edges = edges;
  mesh2->nEdges = nEdges;
  mesh2->neighbors = neighbors;
  mesh2->incidentFaces = incidentFaces;
  mesh2->eta = eta;
  mesh2->h_i = h_i;
}

void boundaryCurve_init_from_meshpy(boundaryCurve *h_i, size_t nEdges, char const *pathBoundary) {
  // structure of the file from meshpy: first the number of the edge, then
  // the tangent with respect to the first node (2 places), then the tangent with
  // respect to the second node (2 places). Considering sorted nodes
  FILE *fp;
  char *line = NULL;
  int i = 0;
  int numEdge;
  size_t len = 0; // length of each line of the file
  ssize_t read; // reading each line in the file
  double row[5]; // temporary row
  row[0] = 0;
  row[1] = 0;
  row[2] = 0;
  row[3] = 0;
  row[4] = 0;

  printf("\ntrying to open boundary tangent file\n");
  // Check if the file exists under that path
  fp = fopen(pathBoundary, "r");
  // if the file doesnt exist of for some reason it can't be opened:
  if( fp == NULL ) {
    printf("No such file, sorry");
    exit(EXIT_FAILURE);
  }
  else {
    printf("\nFile successfully found\n");
  }

  // scan such file
  fp = fopen(pathBoundary, "r");
  while( (read = getline(&line, &len, fp) != -1 ) ) {
    separateARowDb(line, 5, row);
    numEdge = (int)row[0];
    if( numEdge == i ){
      // meaning tha the i-th edge is a curvy boundary
      h_i[i].i_Edge = i;
      h_i[i].B[0][0] = row[1];
      h_i[i].B[0][1] = row[2];
      h_i[i].B[1][0] = row[3];
      h_i[i].B[1][1] = row[4];
      i ++;
    }
    else {
      // meaning that the i-th edge is not a curvy boundary
      h_i[i].i_Edge = i;
      h_i[i].B[0][0] = 0;
      h_i[i].B[0][1] = 0;
      h_i[i].B[1][0] = 0;
      h_i[i].B[1][1] = 0;
      i ++;
    }
  }

  // if we've gone through the file and there are still more edges, add zeros
  while( i < nEdges ){
    // this means that the i-th edge is not in the file and its not a curvy boundary
    h_i[i].i_Edge = i;
    h_i[i].B[0][0] = 0;
    h_i[i].B[0][1] = 0;
    h_i[i].B[1][0] = 0;
    h_i[i].B[1][1] = 0;
    i ++;
  }

  if( i != (int)nEdges ){
    printf("\nError with number of edges, i = %d, nEdges = %zu\n", i, nEdges);
    exit(EXIT_FAILURE);
  }
      
}

void mesh2_init_from_meshpy(mesh2S *mesh2, char const *pathPoints, char const *pathFaces, char const *pathEdges, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathIndices, char const *pathBoundary) {
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
  eta = malloc(nFaces*sizeof(double)); // indices of refraction

  // Now that everything is alloc'ed and declared we can open the files

  // for the points
  read_n2File_double(points[0], pathPoints);
  printf("Read the points\n");
  mesh2->points = points;

  // for the faces
  read_n3File_size_t(faces[0], pathFaces);
  printf("Read the faces\n");
  mesh2->faces = faces;

  // for the edges
  read_n2File_size_t(edges[0], pathEdges);
  printf("Read the edges\n");
  mesh2->edges = edges;

  // read the neighbors structs
  // point neighbors
  neighbors_init(neighbors, pathNeighbors, nPoints);
  printf("Read neighbors\n");
  mesh2->neighbors = neighbors;
  // indices neighbors
  neighbors_init(incidentFaces, pathIncidentFaces, nPoints);
  printf("Read indicent faces to points\n");
  mesh2->incidentFaces = incidentFaces;
  

  // init the boundary curve struct
  boundaryCurve_init_from_meshpy(h_i, nEdges, pathBoundary);
  mesh2->h_i = h_i;

  // for the indices of refraction
  readDbColumn(pathIndices, eta);
  mesh2->eta = eta;

}


void printGeneralInfoMesh(mesh2S *mesh2) {
  printf("\n  GENERAL INFORMATION ABOUT THIS MESH \n\n");
  printf("    Number of points: %zu\n", mesh2->nPoints);
  printf("    Number of faces: %zu\n", mesh2->nFaces);
  printf("    Number of edges: %zu\n", mesh2->nEdges);
}

void printEverythingInMesh(mesh2S *mesh2) {
  int i;
  printf("\n\n---------------------------------------\n");
  printf("---------------------------------------\n");
  printf("\n EVERYTHING CONTAINED IN THIS MESH \n\n");
  printf("\n\n---------------------------------------\n");
  printf("POINTS\n");
  printf("Number of points in the mesh:  %zu.\n", mesh2->nPoints);
  printf("Such points are the following: \n");
  for( i = 0; i < mesh2->nPoints; i++) {
    printf("    %lf   |    %lf   \n", mesh2->points[i][0], mesh2->points[i][1] );
  }
  printf("\n\n---------------------------------------\n");
  printf("FACES\n");
  for( i = 0; i < mesh2->nFaces; i++) {
    printf("Face %d,  points:    %zu   | %zu   | %zu   \n", i, mesh2->faces[i][0], mesh2->faces[i][1], mesh2->faces[i][2]);
  }
  printf("\n\n---------------------------------------\n");
  printf("EDGES\n");
  for( i = 0; i < mesh2->nEdges; i++) {
    printf("Edge %d,  points:    %zu   | %zu   \n", i, mesh2->edges[i][0], mesh2->edges[i][1] );
  }
  printf("\n\n---------------------------------------\n");
  printf("NEIGHBORS\n");
  printf("The neighbors for each indexed point in this mesh are the following: \n");
  printAllNeighbors(mesh2->neighbors, mesh2->nPoints);
  printf("\n\n---------------------------------------\n");
  printf("INCIDENT FACES\n");
  printf("The incident faces for each indexed point in this mesh: \n");
  printAllNeighbors(mesh2->incidentFaces, mesh2->nPoints);
  printf("\n\n---------------------------------------\n");
  printf("FACES INDICES\n");
  for(i = 0; i < mesh2->nFaces; i++){
      printf("Face %d, index %f\n", i, mesh2->eta[i]);
  }
  printf("\n\n---------------------------------------\n");
  printf("TANGENTS TO BOUNDARY\n");
  for(i = 0; i < mesh2->nEdges; i++){
    printf("Edge number %d,    tangent to boundary 0:  %lf | %lf\n", i, mesh2->h_i[i].B[0][0], mesh2->h_i[i].B[0][1] );
    printf("                   tangent to boundary 1:  %lf | %lf\n", mesh2->h_i[i].B[1][0], mesh2->h_i[i].B[1][1] );
  }
}



