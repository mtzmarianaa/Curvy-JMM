#pragma once

#include "triMesh_2D.h"
#include "linAlg.h"
#include "priority_queue.h"

typedef struct eik_grid {
  size_t *start; // the index of the point that is the source (could be multiple, that's why its a pointer)
  size_t nStart; // number of points in start
  mesh2S *mesh2; // NEW MESH struct
  double *eik_vals; // the current Eikonal values for each indexed point in the mesh
  double (*grads)[2]; // gradient of the eikonal
  p_queue *p_queueG; // priority queue struct
  size_t *current_states; // 0 far, 1 trial, 2 valid
} eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, size_t *start, size_t nStart, mesh2S *mesh2);

void eik_grid_initFromFile(eik_gridS *eik_g, size_t *start, size_t nStart,
			   char const *pathPoints, char const *pathFaces,
			    char const *pathEdges, char const *pathEdgesInFace,
			    char const *pathNeighbors, char const *pathIncidentFaces,
			    char const *pathIndices, char const *pathBoundary) ;

void printGeneralInfo(eik_gridS *eik_g);

void printAllInfoMesh(eik_gridS *eik_g);

void simple_Update(double x0[2], double x1[2], double xHat[2], double T0, double T1, double indexRef, double *That2, double *lambda);

void twoStepUpdate(double x0[2], double x1[2], double x2[2], double xHat[2], double T0, double T1, double indexRef_01, double indexRef_02, double *That_step2, double *lambda, double *mu);

void approximateEikonalGradient(double x0[2], double x1[2], double xHat[2], double parameterization, double indexRefraction, double grad[2]);

void initializePointsNear(eik_gridS *eik_g, double rBall);

void updateCurrentValues(eik_gridS *eik_g, int indexToBeUpdated, int parent1, int parent2, double param, double TFound, double indexRefraction);

void updateCurrentValues3(eik_gridS *eik_g, int indexToBeUpdated, int parent0, int parent1, int parent2, double lambda, double mu, double TFound, double indexRefraction);

void addNeighbors_fromAccepted(eik_gridS *eik_g, int indexAccepted);

void popAddNeighbors(eik_gridS *eik_g);

int currentMinIndex(eik_gridS *eik_g);

int nStillInQueue(eik_gridS *eik_g);

void saveComputedValues(eik_gridS *eik_g, const char *pathFile);

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile);

void saveComputedParents(eik_gridS *eik_g, const char *pathFile);

void saveComputedLambdas(eik_gridS *eik_g, const char *pathFile);

// For the methods associated with the cubic hermite interpolation:

void simple_UpdateCubic(eik_gridS *eik_g, double x0[2], double x1[2], double xHat[2], double T0,
			double T1, double grad0[2], double grad1[2],
			int indexAccepted, int xHat_ind, double *indexRef, double *That2, double *lambda);

void twoStepUpdateCubic(eik_gridS *eik_g, double x0[2], double x1[2], double x2[2], double xHat[2], double T0, double T1,
			double grad0[2], double grad1[2], int indexAccepted, int xHat_ind,
			double *indexRef_01, double *indexRef_02, double *That_step2, double *lambda, double *mu);

void addNeighbors_fromAcceptedCubic(eik_gridS *eik_g, int indexAccepted);

void popAddNeighborsCubic(eik_gridS *eik_g);

