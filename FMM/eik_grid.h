#pragma once

#include "triMesh_2D.h"

typedef struct eik_grid eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, int *start, int nStart, triMesh_2Ds *triM_2D) ;

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathBoundaryPoints, char const *pathFacets, char const *pathFaces);

void printGeneralInfo(eik_gridS *eik_g);

double onePointUpdate_eikValue(eik_gridS *eik_g, int indexFrom, int indexTo);

double twoPointUpdate_eikValue(eik_gridS *eik_g, int x0_ind, int x1_ind, int xHat_ind);

void addNeighbors_fromAccepted(eik_gridS *eik_g, int index_accepted);
