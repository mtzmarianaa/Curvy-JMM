/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"
#include "SoSFunction.h"
#include "opti_method.h"
#include "linAlg.h"
#include "updates_2D.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <json-c/json.h> // used for reading the json type string from python

void eik_grid_alloc(eik_gridS **eik_g ) {
  *eik_g = malloc(sizeof(eik_gridS));
  assert(*eik_g != NULL);
}

void eik_grid_dealloc(eik_gridS **eik_g ) {
  free(*eik_g);
  *eik_g = NULL;
}

void triangleFan_alloc(triangleFanS **triFan) {
  *triFan = malloc(sizeof(triangleFanS));
  assert(*triFan != NULL);
}

void triangleFan_dealloc(triangleFanS **triFan) {
  free(*triFan);
  *triFan = NULL;
}

void eik_grid_init( eik_gridS *eik_g, size_t *start, size_t nStart, mesh2S *mesh2) 
{
  // the rest of the parameters, eik_vals, p_queueG, current_states are going to be assigned inside
  eik_g->start = start;
  eik_g->nStart = nStart;
  eik_g->mesh2 = mesh2;

  // we first set all the current eik_vals to infinity, set all the current_states to 0 (far)
  double *eik_vals;
  double (*eik_grad)[2]; // this is a pointer to a list of the gradients of the eikonal
  size_t *current_states;
  size_t *type_update;
  triangleFanS *triFans;
  
  eik_vals = malloc(mesh2->nPoints*sizeof(double)); 
  current_states = malloc(mesh2->nPoints*sizeof(int));
  eik_grad = malloc(2*mesh2->nPoints*sizeof(double)); // each gradient has two coordinates (for each point)
  parents_path = malloc(3*mesh2->nPoints*sizeof(int)); // each node has two parents by index
  lambdas = malloc(mesh2->nPoints*sizeof(double)); // each node was updated with one lambda which is a double
  mus = malloc(mesh2->nPoints*sizeof(double));
  type_update = malloc(mesh2->nPoints*sizeof(int)); // each node has a type of update associated to it
  
  for(int i = 0; i<mesh2->nPoints; i++){
    eik_vals[i] = INFINITY; // set them all to infinity
    current_states[i] = 0; // set them all to far
    eik_grad[i][0] = 0; // initialize all the gradients to zero
    eik_grad[i][1] = 0;
    parents_path[i][0] = i; // initialize both parents as themselves
    parents_path[i][1] = i;
    parents_path[i][2] = i;
    lambdas[i] = 0; // initialize all lambdas to 0 (meaning that its parent is itself, useful for the starting points)
    mus[i] = 0;
    type_update[i] = -1; // initialize the type of update for each node to -1
  }
  eik_g->eik_vals = eik_vals;
  eik_g->current_states = current_states;
  eik_g->eik_grad = eik_grad;
  eik_g->parents_path = parents_path;
  eik_g->lambdas = lambdas;
  eik_g->mus = mus;
  eik_g->type_update = type_update;

  // we initialize the priority queue, all the elements in start are going to be inserted and their current states set to 1
  // notice that we need to add both the index of the starting points AND their eikonal value (which is 0) to the p_queue struct
  p_queue *p_queueG;
  priority_queue_alloc(&p_queueG); // allocate
  priority_queue_init(p_queueG); // initiate
  for(int i = 0; i<nStart; i++){
    insert(p_queueG, 0, start[i]); // insert all the starting points with eikonal value 0
  }
  eik_g->p_queueG = p_queueG;
  assert(&eik_g != NULL); // eik_g should not be null
}

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathFaces, char const *pathIndexRegions, char const *pathBoundaryTan, char const *pathBoundaryChain) {
  // the only difference between this and the previous method is that in here we do need to initialize the Mesh structure
  mesh2S *mesh2;
  triMesh_2Dalloc(&mesh2);
  triMesh2_init_from_meshpy(mesh2, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions, pathBoundaryTan, pathBoundaryChain);
  // and then we can use the previous method
  eik_grid_init( eik_g, start, nStart, mesh2); // voila
}

void printGeneralInfo(eik_gridS *eik_g) {
  printf("\n\n\n\n     GENERAL INFORMATION ON THIS EIKONAL STRUCT     \n\n");
  printf("Number of starting points: %d \n", eik_g->nStart);
  printf("The starting points are: \n");
  for(int i = 0; i<eik_g->nStart; i++) {
    printf("|   %d   |", eik_g->start[i]);
  }
  printGeneralInfoMesh(eik_g->mesh2);
  printf("Current state of priority queue: \n");
  printeik_queue(eik_g->p_queueG);
  printf("\nCurrent Eikonal values: \n");
  for(int i = 0; i<eik_g->mesh2->nPoints; i++){
    double x[2];
    x[0] = eik_g->mesh2->points->x[i];
    x[1] = eik_g->mesh2->points->y[i];
    printf("Index   %d    ||  Coordinates:   (%fl, %fl)    ||  Eikonal:   %fl     ||  Current state:   %d ||  Eik gradient:   (%fl, %fl)||  Parents:   (%d, %d, %d)||  LamOpt:   %fl ||  MuOpt:   %fl  || Update: %d\n", i, x[0], x[1]  , eik_g->eik_vals[i] , eik_g->current_states[i], eik_g->eik_grad[i][0], eik_g->eik_grad[i][1], eik_g->parents_path[i][0], eik_g->parents_path[i][1], eik_g->parents_path[i][2], eik_g->lambdas[i], eik_g->mus[i], eik_g->type_update[i]);
  }
}

void printAllInfoMesh(eik_gridS *eik_g){
  printEverythingInMesh(eik_g->mesh2);
}



void initializePointsNear(eik_gridS *eik_g, double rBall) {
  // THIS JUST SETS THE CURRENT STATE TO VALID AND ADDS THE TRUE EIKONAL
  // given a ball of radius rBall around the initial points we initialize all the points inside those balls with the
  // true value of the eikonal (i.e. the distance times the index of refraction). We are going to assume that
  // all those balls belong to the same regions of indices of refraction
  double xMinxStart[2], xStart[2], xCurrent[2], normCurrent, initialIndexRefraction;
  int indexStart;
  for(int j = 0; j<eik_g->nStart; j++){
    deleteRoot(eik_g->p_queueG);
    indexStart = eik_g->start[j];
    xStart[0] = eik_g->mesh2->points->x[ indexStart ];
    xStart[1] = eik_g->mesh2->points->y[ indexStart];
    initialIndexRefraction = s_function_threeSections(xStart, eik_g->mesh2->indexRegions[ eik_g->mesh2->incidentFaces[indexStart].neis_i[0] ]);
    for(int i = 0; i<eik_g->mesh2->nPoints; i++){
      xCurrent[0] = eik_g->mesh2->points->x[i];
      xCurrent[1] = eik_g->mesh2->points->y[i];
      vec2_subtraction( xCurrent, xStart, xMinxStart );
      normCurrent = l2norm(xMinxStart);
      if(normCurrent < rBall ){
        if( eik_g->current_states[i] == 1 ){
          // if it was previously considered as trial we need to delete this from the queue directly
          delete_findIndex(eik_g->p_queueG, i);
        }
        // if this happens, this point is "close enough" to the starting point so that we can initialize it
        eik_g->current_states[i] = 2; // we initialized it directly
        eik_g->eik_vals[i] = initialIndexRefraction*normCurrent; // add their true Eikonal value
        eik_g->parents_path[i][0] = indexStart;
        eik_g->parents_path[i][1] = indexStart; // both of its parents are the starting point
	eik_g->parents_path[i][2] = indexStart;
        eik_g->lambdas[i] = 0; // could also be 1, same thing
        eik_g->eik_grad[i][0] = initialIndexRefraction*xMinxStart[0]/normCurrent; // update its gradient
        eik_g->eik_grad[i][1] = initialIndexRefraction*xMinxStart[1]/normCurrent;
        addNeighbors_fromAccepted(eik_g, i); // we add its neighbors
      }
    }
  }

}


void approximateEikonalGradient(double xA[2], double xB[2], double xHat[2], double parameterization, double indexRefraction, double grad[2]) {
  // given two points in a base with its corresponding parameterization we approximate the gradient of the Eikonal
  // depending on the type of update performed on xHat xA is x0 and xB could be either x1 or x2
  double MinParametrization, xBasePart1[2], xBasePart2[2], xBase[2], normBase, direction[2];
  MinParametrization = 1 - parameterization;
  scalar_times_2vec( MinParametrization, xA, xBasePart1 );
  scalar_times_2vec( parameterization, xB, xBasePart2 );
  vec2_addition( xBasePart1, xBasePart2, xBase );
  vec2_subtraction(xHat, xBase, direction);
  normBase = l2norm(direction);
  grad[0] = indexRefraction*direction[0]/normBase;
  grad[1] = indexRefraction*direction[1]/normBase;
}









void popAddNeighbors(eik_gridS *eik_g){
  // int nNeighs;
  int minIndex = indexRoot(eik_g->p_queueG);
  // double minEikVal = valueRoot(eik_g->p_queueG);
  // printf("Initially the queue looks like this: \n");
  // printeik_queue(eik_g->p_queueG);
  // printf("\n");
  deleteRoot(eik_g->p_queueG); // delete the root from the priority queue
  // printf("We've eliminated the root from the queue successfully\n");
  eik_g->current_states[minIndex] = 2; // set the newly accepted index to valid
  // printf("We've updated the current state of %d to valid\n", minIndex);
  addNeighbors_fromAccepted(eik_g, minIndex); // add neighbors from the recently accepted index
  // printf("The coordinates of the current minimum value just accepted are: (%fl,%fl)\n", eik_g->mesh2->points->x[minIndex], eik_g->mesh2->points->y[minIndex]);
  // printf("Its neighbors are: \n");
  // nNeighs = eik_g->mesh2->neighbors[minIndex].len; // get the number of neighbors of the minimum index
  // for(int i=0; i<nNeighs; i++){
  //   printf("|     %d     |", eik_g->mesh2->neighbors[minIndex].neis_i[i] );
  // }
  // printf("\n");
  // for(int i=0; i<nNeighs; i++){
  //   printf("|     (%fl, %fl)     |", eik_g->mesh2->points->x[eik_g->mesh2->neighbors[minIndex].neis_i[i]], eik_g->mesh2->points->y[eik_g->mesh2->neighbors[minIndex].neis_i[i]] );
  // }
  // printf("\n");
}

int currentMinIndex(eik_gridS *eik_g) {
  return indexRoot(eik_g->p_queueG);
}

int nStillInQueue(eik_gridS *eik_g) {
  return getSize(eik_g->p_queueG);
}

void saveComputedValues(eik_gridS *eik_g, const char *pathFile) {
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->eik_vals, sizeof(double), eik_g->mesh2->nPoints, fp);
  fclose(fp);
}

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile) {
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->mesh2->nPoints; ++i){
      fwrite(eik_g->eik_grad[i], sizeof(double), 2, fp);
    }

    fclose(fp);
}

void saveComputedParents(eik_gridS *eik_g, const char *pathFile) {
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->mesh2->nPoints; ++i){
      fwrite(eik_g->parents_path[i], sizeof(int), 3, fp);
    }

    fclose(fp);
}

void saveComputedLambdas(eik_gridS *eik_g, const char *pathFile) {
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->lambdas, sizeof(double), eik_g->mesh2->nPoints, fp);
  fclose(fp);
}

void saveComputedMus(eik_gridS *eik_g, const char *pathFile) {
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->mus, sizeof(double), eik_g->mesh2->nPoints, fp);
  fclose(fp);
}

void saveComputedTypesUpdate(eik_gridS *eik_g, const char *pathFile) {
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->type_update, sizeof(int), eik_g->mesh2->nPoints, fp);
  fclose(fp);
}
