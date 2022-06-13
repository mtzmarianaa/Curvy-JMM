/* TRIANGLE 2D MESH STRUCTURE

This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given

*/

#include "triMesh_2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct triMesh_2D {
  int nPoints; // number of points in the triangle mesh
  int nFaces; // number of faces in the triangle mesh
  double *mesh_points; // nPoints x 2 array with the coordinates of the points in the mesh
  int *faces; // M x 3 array with the faces in the mesh (a face is characterized by 3 points, their indeces)
  int *neighbors; // M x 3 array with the neighbor information per triangle (each triangle has 3 triangle neighbours)
  double *boundary_points; // N x 2 array with the coordinates of the points in the mesh on the border
  int *facets; // Nf x 3 array with the facets in the mesh, the faces that correspond to the boundaries
} ;

void triMesh_2Dalloc(triMesh_2Ds **triM_2D) {
    *triM_2D = malloc(sizeof(triMesh_2Ds));
    assert(*triM_2D != NULL);
}

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D) {
    free(*triM_2D);
    *triM_2D = NULL;
}

void triMesh2_init_from_meshpy(triMesh_2Ds **triM_2D, int nPoints, int nFaces, char const *pathPoints,
char const *pathFaces, char const *pathNeighbors, char const *pathBPoints, char const *pathFacets){
    // there are a lot of files needed to be opened
    
}
