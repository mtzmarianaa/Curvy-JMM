#pragma once

typedef struct triMesh_2D triMesh_2Ds;

void triMesh_2Dalloc(triMesh_2Ds **triM_2D);

void triMesh_2Ddalloc(triMesh_2Ds **triM_2D);

void triMesh2_init_from_meshpy(triMesh_2Ds **triM_2D, int nPoints, int nFaces, char const *pathPoints,
char const *pathFaces, char const *pathNeighbors, char const *pathBPoints, char const *pathFacets);