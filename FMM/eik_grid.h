#pragma once

#include "triMesh_2D.h"
#include "linAlg.h"

typedef struct eik_grid eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, int *start, int nStart, triMesh_2Ds *triM_2D) ;

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathFaces, char const *pathIndexRegions);

void printGeneralInfo(eik_gridS *eik_g);

void printAllInfoMesh(eik_gridS *eik_g);

void simple_Update(double x0[2], double x1[2], double xHat[2], double T0, double T1, double indexRef, double *That2, double *lambda);

void initializePointsNear(eik_gridS *eik_g, double rBall);



// void popAddNeighbors(eik_gridS *eik_g);

int currentMinIndex(eik_gridS *eik_g);

int nStillInQueue(eik_gridS *eik_g);

void saveComputedValues(eik_gridS *eik_g, const char *pathFile);

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile);

void saveComputedParents(eik_gridS *eik_g, const char *pathFile);

void saveComputedLambdas(eik_gridS *eik_g, const char *pathFile);
