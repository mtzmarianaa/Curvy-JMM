#pragma once

#include "mesh2D.h"
#include "linAlg.h"
#include "priority_queue.h"
#include "updates_2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

typedef struct triangleFan {
  size_t nRegions;
  double *params; // mu1, lam2, mu2, ..., lamn, mun, lamn1, mun1 =1 should be size 2*nRegions + 1
  double x0[2];
  double T0;
  double grad0[2];
  double x1[2];
  double T1;
  double grad1[2];
  double xHat[2];
  double *listIndices; // in a triangle fan the etas are going to be called indices
  double (*listxk)[2]; // list of points x0, x1, ..., xHat in the triangle fan
  double (*listB0k)[2];
  double (*listBk)[2];
  double (*listBkBk1)[2];
  double *paramsCrTop; // initialize and then fill in after optimizing
  double *paramsStTop; // initialize and then fill in after optimizing
  double THat; // T(xHat) found after optimizing
} triangleFanS;

typedef struct eik_grid {
  int *start; // the index of the point that is the source (could be multiple, that's why its a pointer)
  int nStart; // number of points in start
  mesh2S mesh2; // NEW MESH struct
  double *eik_vals; // the current Eikonal values for each indexed point in the mesh
  double (*eik_grad)[2]; // this is a pointer to a list of the gradients of the eikonal
  p_queue *p_queueG; // priority queue struct
  int *current_states; // 0 far, 1 trial, 2 valid
  int (*parents_path)[3]; // this are the two parent nodes (their indices) from which each node has been updated
  double *lambdas; // lambdas from which (using the two parents) the node was updated
  double *mus; // mus from the two step update (if needed)
  int *type_update; // type of update used for that node
} eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, int *start, int nStart, triMesh_2Ds *triM_2D);

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints,
			   char const *pathNeighbors, char const *pathIncidentFaces,
			   char const *pathFaces, char const *pathIndexRegions,
			   char const *pathBoundaryTan, char const *pathBoundaryChain);

void printGeneralInfo(eik_gridS *eik_g);

void printAllInfoMesh(eik_gridS *eik_g);

void initializePointsNear(eik_gridS *eik_g, double rBall);

void approximateEikonalGradient(double xA[2], double xB[2], double xHat[2],
				double parameterization, double indexRefraction, double grad[2]);

void updateCurrentValues(eik_gridS *eik_g, int indexToBeUpdated, int parent0, int parent1,
			 double param, double TFound, double indexRefraction, int typeUpdate);

void updateCurrentValues3(eik_gridS *eik_g, int indexToBeUpdated, int parent0, int parent1,
			  int parent2, double lambda, double mu, double TFound,
			  double indexRefraction, int typeUpdate);

void addNeighbors_fromAccepted(eik_gridS *eik_g, int indexAccepted);


void popAddNeighbors(eik_gridS *eik_g);

int currentMinIndex(eik_gridS *eik_g);

int nStillInQueue(eik_gridS *eik_g);

void saveComputedValues(eik_gridS *eik_g, const char *pathFile);

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile);

void saveComputedParents(eik_gridS *eik_g, const char *pathFile);

void saveComputedLambdas(eik_gridS *eik_g, const char *pathFile);

void saveComputedMus(eik_gridS *eik_g, const char *pathFile);

void saveComputedTypesUpdate(eik_gridS *eik_g, const char *pathFile);


