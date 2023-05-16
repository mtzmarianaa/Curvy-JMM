#pragma once

#include "mesh2D.h"
#include "linAlg.h"
#include "priority_queue.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <json-c/json.h> // used for reading the json type string from python


// things to call out python optimizer
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>

typedef struct fanUpdate {
  // geometric and eikonal information for a triangle fan
  triangleFanS *triFan; // all the geometric information regarding this triangle fan
  double *params; // mu1, lam2, mu2, ..., lamn, mun, lamn1, mun1 =1 should be size 2*nRegions + 1
  double T0;
  double grad0[2];
  double T1;
  double grad1[2];
  size_t nIndCrTop; // number of indices for CrTop
  size_t *indCrTop;
  double *paramsCrTop; // initialize and then fill in after optimizing, length 2*nIndCrTop
  size_t nIndStTop; // number of indices for StTop
  size_t *indStTop;
  double *paramsStTop; // initialize and then fill in after optimizing, length 2*nIndStTop
  double THat; // T(xHat) found after optimizing
  double (*grads)[2]; // gradients computed using all params, paramsCrTop, paramsStTop, length 2*nRegions + 1 + 2*nIndCrTop + 2*nIndStTop
  double (*path)[2]; // path computed using all params, paramsCrTop, paramsStTop, length 2*nRegions + 1 + 2*nIndCrTop + 2*nIndStTop
  double gradHat[2]; // gradient which is going to be used for xHat
} fanUpdateS;

typedef struct eik_grid {
  size_t *start; // the index of the point that is the source (could be multiple, that's why its a pointer)
  size_t nStart; // number of points in start
  mesh2S *mesh2; // NEW MESH struct
  double *eik_vals; // the current Eikonal values for each indexed point in the mesh
  p_queue *p_queueG; // priority queue struct
  size_t *current_states; // 0 far, 1 trial, 2 valid
  fanUpdateS *fanUpdate; // triangle fan of the optimal update associated to each point in the mesh (includes opti params)
} eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void fanUpdate_alloc(fanUpdateS **fanUpdate);

void fanUpdate_dalloc(fanUpdateS **fanUpdate);

void eik_grid_init( eik_gridS *eik_g, size_t *start, size_t nStart, mesh2S *mesh2);

void fanUpdate_init(fanUpdateS *fanUpdate, triangleFanS *triFan, double *params,
		    double T0, double grad0[2], double T1, double grad1[2],
		    size_t nIndCrTop, size_t *indCrTop, double *paramsCrTop,
		    size_t nIndStTop, size_t *indStTop, double *paramsStTop,
		    double THat, double (*grads)[2], double (*path)[2],
		    double gradHat[2]);

void fanUpdate_initPreOpti(fanUpdateS *fanUpdate, triangleFanS *triFan, double T0,
			   double grad0[2], double T1, double grad1[2]);

void eik_grid_initFromFile(eik_gridS *eik_g, size_t *start, size_t nStart, char const *pathPoints, char const *pathFaces,
			    char const *pathEdges, char const *pathEdgesInFace,
			    char const *pathNeighbors, char const *pathIncidentFaces,
			    char const *pathIndices, char const *pathBoundary) ;

void printGeneralInfo(eik_gridS *eik_g);

void printInfoFanUpdate(fanUpdateS *fanUpdate);

void printAllInfoMesh(eik_gridS *eik_g);

//void initializePointsNear(eik_gridS *eik_g, double rBall);

void findEdgesOnValidFront(eik_gridS *eik_g, size_t index0, int indices1[2], int indices2[2], int firstTriangles[2]);

void initTriFan(eik_gridS *eik_g, triangleFanS *triFan,
		size_t index0, size_t index1, size_t index2,
		size_t indexHat, size_t firstTriangle, double angleMax) ;

void createJSONinput(fanUpdateS *fanUpdate, char **input_json);

void createJSONFile(fanUpdateS *fanUpdate, char const *path);

void deserializeJSONoutput(fanUpdateS *fanUpdate, json_object *output_obj);

void optimizeTriangleFan_wPython(fanUpdateS *fanUpdate);

void updateOneWay(eik_gridS *eik_g, size_t index0, size_t index1, size_t index2,
		  int indexStop_int, size_t firstTriangle);

void addNeighbors_fromAccepted(eik_gridS *eik_g, size_t minIndex);

//void popAddNeighbors(eik_gridS *eik_g);

int currentMinIndex(eik_gridS *eik_g);

int nStillInQueue(eik_gridS *eik_g);

void saveComputedValues(eik_gridS *eik_g, const char *pathFile) ;

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile) ;

