/* TRIANGLE 2D MESH STRUCTURE

This is the 2d triangle mesh structure. It assumes that an output from 
meshpy is given

*/

#include "triMesh_2D.h"
#include "coord.h"
#include "facets.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct triMesh_2D {
    // ask if everything in here is usefull/necessary or if its too much
    // inspiration: what we might need + what might be useful to plot the mesh using triplot in Python
  coordS *points; // these are ALL the coordinates + number of points in the mesh
  coordS *boundaryPoints;  // these are just the coordinates of the boundary points + number of boundary points
  facetsS *facets; // these are the "instructions on how to connect indexed dots"
  int nFaces; // number of faces in the triangle mesh
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
