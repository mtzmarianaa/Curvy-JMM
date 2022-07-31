#pragma once

#include "triMesh_2D.h"
#include "linAlg.h"

typedef struct eik_grid eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, int *start, int nStart, triMesh_2Ds *triM_2D) ;

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathBoundaryPoints, char const *pathFacets, char const *pathFaces, char const *pathIndexRegions);

void printGeneralInfo(eik_gridS *eik_g);

void onePointUpdate_eikValue(eik_gridS *eik_g, int indexFrom, int indexTo, double *That1, int *regionIndex);

void twoPointUpdate_eikValue(eik_gridS *eik_g, int x0_ind, int x1_ind, int xHat_ind, double xlam[2], double *That2, int *regionIndex);

void addNeighbors_fromAccepted(eik_gridS *eik_g, int index_accepted);

void update_afterAccepted(eik_gridS *eik_g, int index_accepted);

void popAddNeighbors(eik_gridS *eik_g);

int currentMinIndex(eik_gridS *eik_g);

int nStillInQueue(eik_gridS *eik_g);

void saveComputedValues(eik_gridS *eik_g, const char *pathFile);

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile);

void savePathsTaken(eik_gridS *eik_g, const char *pathFile);
