
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

void triangleFan_alloc(triangleFanS **triFan) {
  *triFan = malloc(sizeof(triangleFanS));
  assert(*triFan != NULL);
}

void triangleFan_dalloc(triangleFanS **triFan) {
  free(*triFan);
  *triFan = NULL;
}

void boundaryCurve_init(boundaryCurve *h_i, size_t i_Edge, double B[2][2]) {
  h_i->i_Edge = i_Edge;
  h_i->B[0][0] = B[0][0];
  h_i->B[0][1] = B[0][1];
  h_i->B[1][0] = B[1][0];
  h_i->B[1][1] = B[1][1];
}

void mesh2_init(mesh2S *mesh2, double (*points)[2], size_t nPoints,
		size_t (*faces)[3], size_t nFaces, size_t (*edges)[2],
		size_t (*edgesInFace)[3], size_t nEdges, neighborsRS *neighbors,
		neighborsRS *incidentFaces, double *eta, boundaryCurve *h_i) {
  mesh2->points = points;
  mesh2->nPoints = nPoints;
  mesh2->faces = faces;
  mesh2->nFaces = nFaces;
  mesh2->edges = edges;
  mesh2->nEdges = nEdges;
  mesh2->edgesInFace = edgesInFace;
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
    while( i < numEdge ) {
      // meaning that the i-th edge is not a curvy boundary
      h_i[i].i_Edge = i;
      h_i[i].B[0][0] = 0;
      h_i[i].B[0][1] = 0;
      h_i[i].B[1][0] = 0;
      h_i[i].B[1][1] = 0;
      i ++;
    }
    if( numEdge == i ){
      // meaning tha the i-th edge is a curvy boundary
      h_i[i].i_Edge = i;
      h_i[i].B[0][0] = row[1];
      h_i[i].B[0][1] = row[2];
      h_i[i].B[1][0] = row[3];
      h_i[i].B[1][1] = row[4];
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

void mesh2_init_from_meshpy(mesh2S *mesh2, char const *pathPoints, char const *pathFaces,
			    char const *pathEdges, char const *pathEdgesInFace,
			    char const *pathNeighbors, char const *pathIncidentFaces,
			    char const *pathIndices, char const *pathBoundary) {
  // loads information from files generated from python meshpy
  size_t nPoints, nFaces, nEdges, (*faces)[3], (*edges)[2], (*edgesInFace)[3];
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
  edgesInFace = malloc(nFaces*3*sizeof(size_t));
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

  // for the edges in face
  read_n3File_size_t(edgesInFace[0], pathEdgesInFace);
  printf("Read the edges that conform each face\n");
  mesh2->edgesInFace = edgesInFace;

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
  printf("EDGES IN FACES\n");
  for( i = 0; i < mesh2->nFaces; i++) {
    printf("Face %d,  edges in this face:    %zu   | %zu   | %zu   \n", i, mesh2->edgesInFace[i][0], mesh2->edgesInFace[i][1], mesh2->edgesInFace[i][2] );
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


// we need a few methods to initialize triangle fans in the correct way
// notice that given a point we can define two triangle fans

void twoTrianglesFromEdge(mesh2S *mesh2, size_t index0, size_t index1,
			  size_t possibleTriangles[2], size_t possibleThirdVertices[2] ) {
  // given a triangle mesh and two indices indices defining two nodes in the mesh,
  // output the two adjacent triangles to that edge and the two other possible third vertices
  size_t j, currentTriangle;
  j = 0;
  for( int i = 0; i<mesh2->incidentFaces[index0].len; i++ ){
    // circle around all the incident faces to the node associated with index0
    currentTriangle = mesh2->incidentFaces[index0].neis_i[i]; // current incident face to x0
    if( (mesh2->faces[currentTriangle][0] == index0) & (mesh2->faces[currentTriangle][1] == index1) ) {
      possibleTriangles[j] = currentTriangle;
      possibleThirdVertices[j] = mesh2->faces[currentTriangle][2];
      j ++;
    }
    else if( (mesh2->faces[currentTriangle][1] == index0) & (mesh2->faces[currentTriangle][0] == index1) ) {
      possibleTriangles[j] = currentTriangle;
      possibleThirdVertices[j] = mesh2->faces[currentTriangle][2];
      j ++ ;
    }
    else if( (mesh2->faces[currentTriangle][0] == index0) & (mesh2->faces[currentTriangle][2] == index1) ) {
      possibleTriangles[j] = currentTriangle;
      possibleThirdVertices[j] = mesh2->faces[currentTriangle][1];
      j ++;
    }
    else if( (mesh2->faces[currentTriangle][2] == index0) & (mesh2->faces[currentTriangle][0] == index1) ) {
      possibleTriangles[j] = currentTriangle;
      possibleThirdVertices[j] = mesh2->faces[currentTriangle][1];
      j ++;
    }
    else if( (mesh2->faces[currentTriangle][1] == index0) & (mesh2->faces[currentTriangle][2] == index1) ) {
      possibleTriangles[j] = currentTriangle;
      possibleThirdVertices[j] = mesh2->faces[currentTriangle][0];
      j ++;
    }
    else if( (mesh2->faces[currentTriangle][2] == index0) & (mesh2->faces[currentTriangle][1] == index1) ) {
      possibleTriangles[j] = currentTriangle;
      possibleThirdVertices[j] = mesh2->faces[currentTriangle][0];
      j ++;
    }
  }
  if( j == 0 ){
    printf("\n\nCant find incident faces\n\n");
    assert(j != 0);
  }
}


double minEtaFromTwoPoints(mesh2S *mesh2, size_t index0, size_t index1) {
  size_t possibleTriangles[2], possibleThirdVertices[2];
  // using the function defined before
  twoTrianglesFromEdge(mesh2, index0, index1, possibleTriangles, possibleThirdVertices);
  // compare the etas
  double eta0, eta1;
  eta0 = mesh2->eta[possibleTriangles[0]];
  eta1 = mesh2->eta[possibleTriangles[1]];
  if( eta0 < eta1 ){
    return eta0;
  }
  else {
    return eta1;
  }
}

size_t faceBetween3Points(mesh2S *mesh2, size_t index0, size_t index1, size_t index2) {
  // given 3 points that share a face it outputs the index of such face found
  // if there is no such face it outputs -1
  size_t faceIndex;
  int faceFound;
  faceFound = -1; // because we haven't found the desired face yet
  size_t currentFace;
  for( int i = 0; i<mesh2->incidentFaces[index0].len; i++) {
    currentFace = mesh2->incidentFaces[index0].neis_i[i];
    if( (mesh2->faces[currentFace][0] == index1) | (mesh2->faces[currentFace][1] == index1) | (mesh2->faces[currentFace][2] == index1) ) {
      // index0 and index1 share a face/triangle
      if( (mesh2->faces[currentFace][0] == index2) | (mesh2->faces[currentFace][1] == index2) | (mesh2->faces[currentFace][2] == index2) ){
	// index0, index1, index2 share this face/triangle
	faceIndex = currentFace;
	faceFound = 0; // because we've found a face
	break;
      }
    }
  }
  if( faceFound == -1 ){
    printf("\n\nIndex0: %zu, Index1: %zu, Index2: %zu   DONT SHARE A FACE\n\n", index0, index1, index2);
    assert(faceFound != 1);
  }
  return faceIndex;
}

void triangleFan_initFromIndices(triangleFanS *triFan, mesh2S *mesh2, size_t nRegions,
				 size_t index0, size_t index1,
				 size_t indexHat, size_t *listIndicesNodes) {
  // from information from eikgrid we set up a triangle fan
  size_t *listFaces, faceBetweenPoints, *listEdges;
  double x0[2], x1[2], x2[2], xHat[2], *listIndices;
  double (*listxk)[2], (*listB0k)[2], (*listBk)[2], (*listBkBk1)[2];
  triFan->nRegions = nRegions;
  x0[0] = mesh2->points[index0][0];
  x0[1] = mesh2->points[index0][1];
  x1[0] = mesh2->points[index1][0];
  x1[1] = mesh2->points[index1][1];
  x2[0] = mesh2->points[listIndicesNodes[2]][0];
  x2[1] = mesh2->points[listIndicesNodes[2]][1];
  xHat[0] = mesh2->points[indexHat][0];
  xHat[1] = mesh2->points[indexHat][1];
  listFaces = malloc(nRegions*sizeof(size_t));
  listIndices = malloc((2*nRegions + 1)*sizeof(double));
  listEdges = malloc((2*nRegions + 1)*sizeof(size_t));
  listxk = malloc(2*(nRegions + 2)*sizeof(double));
  listB0k = malloc(2*(nRegions + 1)*sizeof(double));
  listBk = malloc(2*(nRegions + 1)*sizeof(double));
  listBkBk1 = malloc(2*(2*nRegions)*sizeof(double));
  // start filling the information in
  double xk[2], xk1[2];
  size_t indexk, indexk1, edge0, edge1, edge2;
  size_t edge0k, edge0k1, edgekk1, possibleTriangles[2], possibleThirdVertices[2];
  int j, i;
  for( i = 0; i <(nRegions + 2); i++){
    listxk[i][0] = mesh2->points[listIndicesNodes[i]][0];
    listxk[i][1] = mesh2->points[listIndicesNodes[i]][1];
  }

  
  
  for( i = 0; i<nRegions; i++){
    j = 2*i; // useful for the BkBk1
    // iterate
    indexk = listIndicesNodes[i+1];
    indexk1 = listIndicesNodes[i+2];
    // get the coordinates of xk and xk1
    xk[0] = mesh2->points[indexk][0];
    xk[1] = mesh2->points[indexk][1];
    xk1[0] = mesh2->points[indexk1][0];
    xk1[1] = mesh2->points[indexk1][1];
    // get the index of the triangle x0 xk xk1
    faceBetweenPoints = faceBetween3Points(mesh2, index0, indexk, indexk1);
    listFaces[i] = faceBetweenPoints; // add this information
    listIndices[i] = mesh2->eta[faceBetweenPoints]; // add this index of refraction
    // We need to add the index of refraction of the triangle on the outside
    twoTrianglesFromEdge(mesh2, indexk, indexk1, possibleTriangles, possibleThirdVertices);
    // we don't want faceBetweenPoints, we want the other one
    if( possibleTriangles[0] != faceBetweenPoints){
      listIndices[nRegions + i+1] = mesh2->eta[possibleTriangles[0]];
    }
    else{
      listIndices[nRegions + i+1] = mesh2->eta[possibleTriangles[1]];
    }
    // get the 3 that make up this triangle
    edge0 = mesh2->edgesInFace[faceBetweenPoints][0];
    edge1 = mesh2->edgesInFace[faceBetweenPoints][1];
    edge2 = mesh2->edgesInFace[faceBetweenPoints][2];
    printf("edge0 %zu\n", edge0);
    printf("edge1 %zu \n", edge1);
    printf("edge2 %zu \n", edge2);
    // we need to know which edge is which AND IN WHICH DIRECTION THEY ARE SAVED
    // EDGE0
    if( (mesh2->edges[edge0][0] == index0) & (mesh2->edges[edge0][1] == indexk) ) {
      // we are in the x0 xk edge
      listEdges[i] = edge0;
      listB0k[i][0] = mesh2->h_i[edge0].B[0][0];
      listB0k[i][1] = mesh2->h_i[edge0].B[0][1];
      listBk[i][0] = mesh2->h_i[edge0].B[1][0];
      listBk[i][1] = mesh2->h_i[edge0].B[1][1];
    }
    else if( (mesh2->edges[edge0][1] == index0) & (mesh2->edges[edge0][0] == indexk) ) {
      // we are in the xk x0 edge
      listEdges[i] = edge0;
      listB0k[i][0] = mesh2->h_i[edge0].B[1][0];
      listB0k[i][1] = mesh2->h_i[edge0].B[1][1];
      listBk[i][0] = mesh2->h_i[edge0].B[0][0];
      listBk[i][1] = mesh2->h_i[edge0].B[0][1];
    }
    else if( (mesh2->edges[edge0][0] == index0) & (mesh2->edges[edge0][1] == indexk1) ) {
      // we are in the x0 xk1 edge
      listEdges[i+1] = edge0;
      listB0k[i+1][0] = mesh2->h_i[edge0].B[0][0];
      listB0k[i+1][1] = mesh2->h_i[edge0].B[0][1];
      listBk[i+1][0] = mesh2->h_i[edge0].B[1][0];
      listBk[i+1][1] = mesh2->h_i[edge0].B[1][1];
    }
    else if( (mesh2->edges[edge0][1] == index0) & (mesh2->edges[edge0][0] == indexk1) ) {
      // we are in the xk1 x0 edge
      listEdges[i+1] = edge0;
      listB0k[i+1][0] = mesh2->h_i[edge0].B[1][0];
      listB0k[i+1][1] = mesh2->h_i[edge0].B[1][1];
      listBk[i+1][0] = mesh2->h_i[edge0].B[0][0];
      listBk[i+1][1] = mesh2->h_i[edge0].B[0][1];
    }
    else if( (mesh2->edges[edge0][1] == indexk) & (mesh2->edges[edge0][0] == indexk1) ) {
      // we are in the xk xk1 edge
      listEdges[nRegions + i + 1] = edge0;
      listBkBk1[j][0] = mesh2->h_i[edge0].B[1][0];
      listBkBk1[j][1] = mesh2->h_i[edge0].B[1][1];
      listBkBk1[j+1][0] = mesh2->h_i[edge0].B[0][0];
      listBkBk1[j+1][1] = mesh2->h_i[edge0].B[0][1];
    }
    else if( (mesh2->edges[edge0][0] == indexk) & (mesh2->edges[edge0][1] == indexk1) ) {
      // we are in the xk1 xk edge
      listEdges[nRegions + i + 1] = edge0;
      listBkBk1[j+1][0] = mesh2->h_i[edge0].B[1][0];
      listBkBk1[j+1][1] = mesh2->h_i[edge0].B[1][1];
      listBkBk1[j][0] = mesh2->h_i[edge0].B[0][0];
      listBkBk1[j][1] = mesh2->h_i[edge0].B[0][1];
    }

    // EDGE1
    if( (mesh2->edges[edge1][0] == index0) & (mesh2->edges[edge1][1] == indexk) ) {
      // we are in the x0 xk edge
      listEdges[i] = edge1;
      listB0k[i][0] = mesh2->h_i[edge1].B[0][0];
      listB0k[i][1] = mesh2->h_i[edge1].B[0][1];
      listBk[i][0] = mesh2->h_i[edge1].B[1][0];
      listBk[i][1] = mesh2->h_i[edge1].B[1][1];
    }
    else if( (mesh2->edges[edge1][1] == index0) & (mesh2->edges[edge1][0] == indexk) ) {
      // we are in the xk x0 edge
      listEdges[i] = edge1;
      listB0k[i][0] = mesh2->h_i[edge1].B[1][0];
      listB0k[i][1] = mesh2->h_i[edge1].B[1][1];
      listBk[i][0] = mesh2->h_i[edge1].B[0][0];
      listBk[i][1] = mesh2->h_i[edge1].B[0][1];
    }
    else if( (mesh2->edges[edge1][0] == index0) & (mesh2->edges[edge1][1] == indexk1) ) {
      // we are in the x0 xk1 edge
      listEdges[i+1] = edge1;
      listB0k[i+1][0] = mesh2->h_i[edge1].B[0][0];
      listB0k[i+1][1] = mesh2->h_i[edge1].B[0][1];
      listBk[i+1][0] = mesh2->h_i[edge1].B[1][0];
      listBk[i+1][1] = mesh2->h_i[edge1].B[1][1];
    }
    else if( (mesh2->edges[edge1][1] == index0) & (mesh2->edges[edge1][0] == indexk1) ) {
      // we are in the xk1 x0 edge
      listEdges[i+1] = edge1;
      listB0k[i+1][0] = mesh2->h_i[edge1].B[1][0];
      listB0k[i+1][1] = mesh2->h_i[edge1].B[1][1];
      listBk[i+1][0] = mesh2->h_i[edge1].B[0][0];
      listBk[i+1][1] = mesh2->h_i[edge1].B[0][1];
    }
    else if( (mesh2->edges[edge1][1] == indexk) & (mesh2->edges[edge1][0] == indexk1) ) {
      // we are in the xk xk1 edge
      listEdges[nRegions + i + 1] = edge1;
      listBkBk1[j][0] = mesh2->h_i[edge1].B[1][0];
      listBkBk1[j][1] = mesh2->h_i[edge1].B[1][1];
      listBkBk1[j+1][0] = mesh2->h_i[edge1].B[0][0];
      listBkBk1[j+1][1] = mesh2->h_i[edge1].B[0][1];
    }
    else if( (mesh2->edges[edge1][0] == indexk) & (mesh2->edges[edge1][1] == indexk1) ) {
      // we are in the xk1 xk edge
      listEdges[nRegions + i + 1] = edge1;
      listBkBk1[j+1][0] = mesh2->h_i[edge1].B[1][0];
      listBkBk1[j+1][1] = mesh2->h_i[edge1].B[1][1];
      listBkBk1[j][0] = mesh2->h_i[edge1].B[0][0];
      listBkBk1[j][1] = mesh2->h_i[edge1].B[0][1];
    }


    // EDGE2
    if( (mesh2->edges[edge2][0] == index0) & (mesh2->edges[edge2][1] == indexk) ) {
      // we are in the x0 xk edge
      listEdges[i] = edge2;
      listB0k[i][0] = mesh2->h_i[edge2].B[0][0];
      listB0k[i][1] = mesh2->h_i[edge2].B[0][1];
      listBk[i][0] = mesh2->h_i[edge2].B[1][0];
      listBk[i][1] = mesh2->h_i[edge2].B[1][1];
    }
    else if( (mesh2->edges[edge2][1] == index0) & (mesh2->edges[edge2][0] == indexk) ) {
      // we are in the xk x0 edge
      listEdges[i] = edge2;
      listB0k[i][0] = mesh2->h_i[edge2].B[1][0];
      listB0k[i][1] = mesh2->h_i[edge2].B[1][1];
      listBk[i][0] = mesh2->h_i[edge2].B[0][0];
      listBk[i][1] = mesh2->h_i[edge2].B[0][1];
    }
    else if( (mesh2->edges[edge2][0] == index0) & (mesh2->edges[edge2][1] == indexk1) ) {
      // we are in the x0 xk1 edge
      listEdges[i+1] = edge2;
      listB0k[i+1][0] = mesh2->h_i[edge2].B[0][0];
      listB0k[i+1][1] = mesh2->h_i[edge2].B[0][1];
      listBk[i+1][0] = mesh2->h_i[edge2].B[1][0];
      listBk[i+1][1] = mesh2->h_i[edge2].B[1][1];
    }
    else if( (mesh2->edges[edge2][1] == index0) & (mesh2->edges[edge2][0] == indexk1) ) {
      // we are in the xk1 x0 edge
      listEdges[i+1] = edge2;
      listB0k[i+1][0] = mesh2->h_i[edge2].B[1][0];
      listB0k[i+1][1] = mesh2->h_i[edge2].B[1][1];
      listBk[i+1][0] = mesh2->h_i[edge2].B[0][0];
      listBk[i+1][1] = mesh2->h_i[edge2].B[0][1];
    }
    else if( (mesh2->edges[edge2][1] == indexk) & (mesh2->edges[edge2][0] == indexk1) ) {
      // we are in the xk xk1 edge
      listEdges[nRegions + i + 1] = edge2;
      listBkBk1[j][0] = mesh2->h_i[edge2].B[1][0];
      listBkBk1[j][1] = mesh2->h_i[edge2].B[1][1];
      listBkBk1[j+1][0] = mesh2->h_i[edge2].B[0][0];
      listBkBk1[j+1][1] = mesh2->h_i[edge2].B[0][1];
    }
    else if( (mesh2->edges[edge2][0] == indexk) & (mesh2->edges[edge2][1] == indexk1) ) {
      // we are in the xk1 xk edge
      listEdges[nRegions + i + 1] = edge2;
      listBkBk1[j+1][0] = mesh2->h_i[edge2].B[1][0];
      listBkBk1[j+1][1] = mesh2->h_i[edge2].B[1][1];
      listBkBk1[j][0] = mesh2->h_i[edge2].B[0][0];
      listBkBk1[j][1] = mesh2->h_i[edge2].B[0][1];
    }
    


    
  }


  // add information regarding the edge x0 xHat
  size_t indexHatM1;
  indexHatM1 = listIndices[nRegions];
  twoTrianglesFromEdge(mesh2, index0, indexHat, possibleTriangles, possibleThirdVertices);
  if( possibleTriangles[0] != listFaces[nRegions - 1] ){
    // this is the triangle we want
    listIndices[nRegions] = mesh2->eta[possibleTriangles[0]];
  }
  else{
    listIndices[nRegions] = mesh2->eta[possibleTriangles[1]];
  }
  // for the edges
  edge0 = mesh2->edgesInFace[listFaces[nRegions-1]][0];
  edge1 = mesh2->edgesInFace[listFaces[nRegions-1]][1];
  edge2 = mesh2->edgesInFace[listFaces[nRegions-1]][2];
  printf("\n\nedge0 %zu\n", edge0);
  printf("edge1 %zu \n", edge1);
  printf("edge2 %zu \n", edge2);
  if( (mesh2->edges[edge0][0] == index0) & (mesh2->edges[edge0][1] == indexHat)) {
    listEdges[nRegions] = edge0;
    listB0k[nRegions][0] = mesh2->h_i[edge0].B[0][0];
    listB0k[nRegions][1] = mesh2->h_i[edge0].B[0][1];
    listBk[nRegions][0] = mesh2->h_i[edge0].B[1][0];
    listBk[nRegions][1] = mesh2->h_i[edge0].B[1][1];
  }
  else if( (mesh2->edges[edge0][0] == indexHat) & (mesh2->edges[edge0][1] == index0)) {
    listEdges[nRegions] = edge0;
    listBk[nRegions][0] = mesh2->h_i[edge0].B[0][0];
    listBk[nRegions][1] = mesh2->h_i[edge0].B[0][1];
    listB0k[nRegions][0] = mesh2->h_i[edge0].B[1][0];
    listB0k[nRegions][1] = mesh2->h_i[edge0].B[1][1];
  }
  if( (mesh2->edges[edge1][0] == index0) & (mesh2->edges[edge1][1] == indexHat)) {
    listEdges[nRegions] = edge1;
    listB0k[nRegions][0] = mesh2->h_i[edge1].B[0][0];
    listB0k[nRegions][1] = mesh2->h_i[edge1].B[0][1];
    listBk[nRegions][0] = mesh2->h_i[edge1].B[1][0];
    listBk[nRegions][1] = mesh2->h_i[edge1].B[1][1];
  }
  else if( (mesh2->edges[edge1][0] == indexHat) & (mesh2->edges[edge1][1] == index0)) {
    listEdges[nRegions] = edge1;
    listBk[nRegions][0] = mesh2->h_i[edge1].B[0][0];
    listBk[nRegions][1] = mesh2->h_i[edge1].B[0][1];
    listB0k[nRegions][0] = mesh2->h_i[edge1].B[1][0];
    listB0k[nRegions][1] = mesh2->h_i[edge1].B[1][1];
  }
  if( (mesh2->edges[edge2][0] == index0) & (mesh2->edges[edge2][1] == indexHat)) {
    listEdges[nRegions] = edge2;
    listB0k[nRegions][0] = mesh2->h_i[edge2].B[0][0];
    listB0k[nRegions][1] = mesh2->h_i[edge2].B[0][1];
    listBk[nRegions][0] = mesh2->h_i[edge2].B[1][0];
    listBk[nRegions][1] = mesh2->h_i[edge2].B[1][1];
  }
  else if( (mesh2->edges[edge2][0] == indexHat) & (mesh2->edges[edge2][1] == index0)) {
    listEdges[nRegions] = edge2;
    listBk[nRegions][0] = mesh2->h_i[edge2].B[0][0];
    listBk[nRegions][1] = mesh2->h_i[edge2].B[0][1];
    listB0k[nRegions][0] = mesh2->h_i[edge2].B[1][0];
    listB0k[nRegions][1] = mesh2->h_i[edge2].B[1][1];
  }
  

  // finally assign everything
  triFan->x0[0] = x0[0];
  triFan->x0[1] = x0[1];
  triFan->x1[0] = x1[0];
  triFan->x1[1] = x1[1];
  triFan->xHat[0] = xHat[0];
  triFan->xHat[1] = xHat[1];
  triFan->listFaces = listFaces;
  triFan->listIndices = listIndices;
  triFan->listEdges = listEdges;
  triFan->listxk = listxk;
  triFan->listB0k = listB0k;
  triFan->listBk = listBk;
  triFan->listBkBk1 = listBkBk1;
  
}


void printEverythingTriFan(triangleFanS *triFan) {
  // we want to print all the information regarding a triangle fan
  int i;
  printf("\n\nInformation regarding this triangle fan:\n");
  printf("Number of regions: %zu\n", triFan->nRegions);
  printf("x0: %lf  |  %lf\n", triFan->x0[0], triFan->x0[1]);
  printf("x1: %lf  |  %lf\n", triFan->x1[0], triFan->x1[1]);
  printf("xHat: %lf  |  %lf\n", triFan->xHat[0], triFan->xHat[1]);
  printf("\n\nList faces: \n");
  for( i = 0; i<triFan->nRegions; i++){
    printf("%zu, ", triFan->listFaces[i]);
  }
  
  printf("\n\nList indices: \n");
  for( i = 0; i<(2*triFan->nRegions + 1); i++){
    printf("%fl, ", triFan->listIndices[i]);
  }
  
  printf("\n\nList Edges: \n");
  for( i = 0; i<(2*triFan->nRegions + 1); i++){
    printf("%zu, ", triFan->listEdges[i]);
  }

  printf("\n\nList xk: \n");
  for( i = 0; i<(triFan->nRegions+2); i++){
    printf("(%fl,  %fl) , ", triFan->listxk[i][0], triFan->listxk[i][1] );
  }

  printf("\n\nList B0k: \n");
  for( i = 0; i<(triFan->nRegions+1); i++){
    printf("(%fl,  %fl) , ", triFan->listB0k[i][0], triFan->listB0k[i][1] );
  }

  printf("\n\nList Bk: \n");
  for( i = 0; i<(triFan->nRegions+1); i++){
    printf("(%fl,  %fl) , ", triFan->listBk[i][0], triFan->listBk[i][1] );
  }

  printf("\n\nList Bkk1: \n");
  for( i = 0; i<(2*triFan->nRegions); i++){
    printf("(%fl,  %fl) , ", triFan->listBkBk1[i][0], triFan->listBkBk1[i][1] );
  }
  
}



