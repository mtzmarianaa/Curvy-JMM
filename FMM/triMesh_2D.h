#pragma once

#include "coord.h"
#include "facets.h"
#include "faces.h"
#include "neighbors.h"
#include "files_methods.h"

typedef struct triMesh_2D triMesh_2Ds;

void triMesh_2Dalloc(triMesh_2Ds **triM_2D);

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D);

void triMesh2_init(triMesh_2Ds *triM_2D, coordS *points, neighborsRS *neighbors, coordS *boundaryPoints, facetsS *facets, facesS *faces, int nPoints);

void triMesh2_init_from_meshpy(triMesh_2Ds *triM_2D, char const *pathPoints, char const *pathNeighbors, char const *pathBoundaryPoints, char const *pathFacets, char const *pathFaces);

void printGeneralInfoMesh(triMesh_2Ds *triM_2D);

void printEverythingInMesh(triMesh_2Ds *triM_2D);