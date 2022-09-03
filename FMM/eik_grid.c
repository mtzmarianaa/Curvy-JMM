/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"
#include "SoSFunction.h"
#include "opti_method.h"
#include "linAlg.h"
#include "path.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

struct eik_grid {
  int *start; // the index of the point that is the source (could be multiple, that's why its a pointer)
  int nStart; // number of points in start
  triMesh_2Ds *triM_2D; // mesh structure that includes coordinates, neighbors, incident faces, etc.
  double *eik_vals; // the current Eikonal values for each indexed point in the mesh
  double (*eik_grad)[2]; // this is a pointer to a list of the gradients of the eikonal
  p_queue *p_queueG; // priority queue struct
  int *current_states; // 0 far, 1 trial, 2 valid
  int (*parents_path)[2]; // this are the two parent nodes (their indices) from which each node has been updated
  double *lambdas; // lambdas from which (using the two parents) the node was updated
} ;

void eik_grid_alloc(eik_gridS **eik_g ) {
  *eik_g = malloc(sizeof(eik_gridS));
  assert(*eik_g != NULL);
}

void eik_grid_dealloc(eik_gridS **eik_g ) {
  free(*eik_g);
  *eik_g = NULL;
}

void eik_grid_init( eik_gridS *eik_g, int *start, int nStart, triMesh_2Ds *triM_2D) 
{
  // the rest of the parameters, eik_vals, p_queueG, current_states are going to be assigned inside
  eik_g->start = start;
  eik_g->nStart = nStart;
  eik_g->triM_2D = triM_2D;

  // we first set all the current eik_vals to infinity, set all the current_states to 0 (far)
  double *eik_vals;
  double (*eik_grad)[2]; // this is a pointer to a list of the gradients of the eikonal
  int (*parents_path)[2]; // this is a pointer to a list of the parents that are used to update that node (by indices)
  double *lambdas; // pointer to the lambdas used to update each node using their parents
  int *current_states;
  eik_vals = malloc(triM_2D->nPoints*sizeof(double)); 
  current_states = malloc(triM_2D->nPoints*sizeof(int));
  eik_grad = malloc(2*triM_2D->nPoints*sizeof(double)); // each gradient has two coordinates (for each point)
  parents_path = malloc(2*triM_2D->nPoints*sizeof(int)); // each node has two parents by index
  lambdas = malloc(triM_2D->nPoints*sizeof(double)); // each node was updated with one lambda which is a double
  for(int i = 0; i<triM_2D->nPoints; i++){
    eik_vals[i] = INFINITY; // set them all to infinity
    current_states[i] = 0; // set them all to far
    eik_grad[i][0] = 0; // initialize all the gradients to zero
    eik_grad[i][1] = 0;
    parents_path[i][0] = i; // initialize both parents as themselves
    parents_path[i][1] = i;
    lambdas[i] = 0; // initialize all lambdas to 0 (meaning that its parent is itself, useful for the starting points)
  }
  eik_g->eik_vals = eik_vals;
  eik_g->current_states = current_states;
  eik_g->eik_grad = eik_grad;
  eik_g->parents_path = parents_path;
  eik_g->lambdas = lambdas;

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

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathFaces, char const *pathIndexRegions) {
  // the only difference between this and the previous method is that in here we do need to initialize the Mesh structure
  triMesh_2Ds *triM_2D;
  triMesh_2Dalloc(&triM_2D);
  triMesh2_init_from_meshpy(triM_2D, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions);
  // and then we can use the previous method
  eik_grid_init( eik_g, start, nStart, triM_2D); // voila
}

void printGeneralInfo(eik_gridS *eik_g) {
  printf("\n\n\n\n     GENERAL INFORMATION ON THIS EIKONAL STRUCT     \n\n");
  printf("Number of starting points: %d \n", eik_g->nStart);
  printf("The starting points are: \n");
  for(int i = 0; i<eik_g->nStart; i++) {
    printf("|   %d   |", eik_g->start[i]);
  }
  printGeneralInfoMesh(eik_g->triM_2D);
  printf("Current state of priority queue: \n");
  printeik_queue(eik_g->p_queueG);
  printf("\nCurrent Eikonal values: \n");
  for(int i = 0; i<eik_g->triM_2D->nPoints; i++){
    double x[2];
    x[0] = eik_g->triM_2D->points->x[i];
    x[1] = eik_g->triM_2D->points->y[i];
    printf("Index   %d    ||  Coordinates:   (%fl, %fl)    ||  Eikonal:   %fl     ||  Current state:   %d ||  Eik gradient:   (%fl, %fl)||  Parents:   (%d, %d)||  LamOpt:   %fl \n", i, x[0], x[1]  , eik_g->eik_vals[i] , eik_g->current_states[i], eik_g->eik_grad[i][0], eik_g->eik_grad[i][1], eik_g->parents_path[i][0], eik_g->parents_path[i][1], eik_g->lambdas[i]);
  }
}

void printAllInfoMesh(eik_gridS *eik_g){
  printEverythingInMesh(eik_g->triM_2D);
}

void simple_Update(double x0[2], double x1[2], double xHat[2], double T0, double T1, double indexRef, double *That2, double *lambda) {
  // this is a simple two point update, meaning that there are no change in regions considered
  // uses optimization to find xLam, which is on the line defined from x0 to x1 parametrized by lambda

  // first we find the optimal lambda
  double tol, lambda0, lambda1;
  int maxIter;
  tol = 0.00001;
  maxIter = 30;
  lambda0 = 0;
  lambda1 = 1;

  *lambda = secant_2D(lambda0, lambda1, T0, T1, x0, x1, xHat, tol, maxIter, indexRef); // optimal lambda found
  // compute the approximation to That2
  *That2 = eikApproxLin(T1, T0, *lambda, x0, x1, xHat, indexRef);

}

void twoStepUpdate(double x0[2], double x1[2], double x2[2], double xHat[2], double T0, double T1, double indexRef_01, double indexRef_02, double *That_step2, double *mu){
  // a two step update, there is a change of the index of refraction from indexRef_01 to indexRef_02 from the triangles
  // x2 x0 x1 to the triangle xHat x0 x2
  double tol, optimizers[2];
  int maxIter;
  tol = 0.00001;
  maxIter = 30;
  // we use the projected gradient descent method to get both lambda and mu optimal
  projectedGradientDescent(optimizers, T0, T1, x0, x1, x2, xHat, tol, maxIter, indexRef_01, indexRef_02);
  // use them to get the minimum path and the approximated eikonal
  *That_step2 = eikApproxLin_2Regions(T0, T1, optimizers[0], optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
  // and extract the optimal mu from the optimizers found (we're going to use this to approximate the gradient)
  *mu = optimizers[1];
}

void approximateEikonalGradient(double x0[2], double x1[2], double parameterization, double indexRefraction, double grad[2]){
  // given two points in a base with its corresponding parameterization we approximate the gradient of the Eikonal
  double MinParametrization, xBasePart1[2], xBasePart2[2], xBase[2], normBase;
  MinParametrization = 1 - parameterization;
  scalar_times_2vec( MinParametrization, x0, xBasePart1 );
  scalar_times_2vec( parameterization, x1, xBasePart2 );
  vec2_addition( xBasePart1, xBasePart2, xBase );
  normBase = l2norm(xBase);
  grad[0] = indexRefraction*xBase[0]/normBase;
  grad[1] = indexRefraction*xBase[1]/normBase;
}


void initializePointsNear(eik_gridS *eik_g, double rBall) {
  // THIS JUST SETS THE CURRENT STATE TO VALID AND ADDS THE TRUE EIKONAL
  // given a ball of radius rBall around the initial points we initialize all the points inside those balls with the
  // true value of the eikonal (i.e. the distance times the index of refraction). We are going to assume that
  // all those balls belong to the same regions of indices of refraction
  // no neighbors are added or considered here
  double xMinxStart[2], xStart[2], xCurrent[2], normCurrent, initialIndexRefraction;
  int indexStart;
  for(int j = 0; j<eik_g->nStart; j++){
    indexStart = eik_g->start[j];
    xStart[0] = eik_g->triM_2D->points->x[ indexStart ];
    xStart[1] = eik_g->triM_2D->points->y[ indexStart];
    initialIndexRefraction = s_function_threeSections(xStart, eik_g->triM_2D->indexRegions[ eik_g->triM_2D->incidentFaces[indexStart].neis_i[0] ]);
    for(int i = 0; i<eik_g->triM_2D->nPoints; i++){
      xCurrent[0] = eik_g->triM_2D->points->x[i];
      xCurrent[1] = eik_g->triM_2D->points->y[i];
      vec2_substraction( xCurrent, xStart, xMinxStart );
      normCurrent = l2norm(xMinxStart);
      if(normCurrent < rBall ){
        // if this happens, this point is "close enough" to the starting point so that we can initialize it
        eik_g->current_states[i] = 2; // we initialized it directly
        eik_g->eik_vals[i] = initialIndexRefraction*normCurrent; // add their true Eikonal value
        eik_g->parents_path[i][0] = indexStart;
        eik_g->parents_path[i][1] = indexStart; // both of its parents are the starting point
        eik_g->lambdas[i] = 0; // could also be 1, same thing
        eik_g->eik_grad[i][0] = initialIndexRefraction*xMinxStart[0]/normCurrent; // update its gradient
        eik_g->eik_grad[i][1] = initialIndexRefraction*xMinxStart[1]/normCurrent;
      }
    }
  }
  // now we need to add the neighbors of all the points that are 

}


// void popAddNeighbors(eik_gridS *eik_g){
//   int nNeighs;
//   int minIndex = indexRoot(eik_g->p_queueG);
//   double minEikVal = valueRoot(eik_g->p_queueG);
//   printf("Initially the queue looks like this: \n");
//   printeik_queue(eik_g->p_queueG);
//   printf("\n");
//   deleteRoot(eik_g->p_queueG); // delete the root from the priority queue
//   eik_g->current_states[minIndex] = 2; // set the newly accepted index to valid
//   eik_g->eik_vals[minIndex] = minEikVal; // add the computed eikonal value to the list of eikonal values
//   addNeighbors_fromAccepted(eik_g, minIndex); // add neighbors from the recently accepted index
//   printf("The coordinates of the current minimum value just accepted are: (%fl,%fl)\n", eik_g->triM_2D->points->x[minIndex], eik_g->triM_2D->points->y[minIndex]);
//   printf("Its neighbors are: \n");
//   nNeighs = eik_g->triM_2D->neighbors[minIndex].len; // get the number of neighbors of the minimum index
//   for(int i=0; i<nNeighs; i++){
//     printf("|     %d     |", eik_g->triM_2D->neighbors[minIndex].neis_i[i] );
//   }
//   printf("\n");
//   for(int i=0; i<nNeighs; i++){
//     printf("|     (%fl, %fl)     |", eik_g->triM_2D->points->x[eik_g->triM_2D->neighbors[minIndex].neis_i[i]], eik_g->triM_2D->points->y[eik_g->triM_2D->neighbors[minIndex].neis_i[i]] );
//   }
//   printf("\n");
// }

int currentMinIndex(eik_gridS *eik_g) {
  return indexRoot(eik_g->p_queueG);
}

int nStillInQueue(eik_gridS *eik_g) {
  return getSize(eik_g->p_queueG);
}

void saveComputedValues(eik_gridS *eik_g, const char *pathFile){
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->eik_vals, sizeof(double), eik_g->triM_2D->nPoints, fp);
  fclose(fp);
}

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile){
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->triM_2D->nPoints; ++i){
      fwrite(eik_g->eik_grad[i], sizeof(double), 2, fp);
    }

    fclose(fp);
}

void saveComputedParents(eik_gridS *eik_g, const char *pathFile){
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->triM_2D->nPoints; ++i){
      fwrite(eik_g->parents_path[i], sizeof(int), 2, fp);
    }

    fclose(fp);
}

void saveComputedLambdas(eik_gridS *eik_g, const char *pathFile){
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->lambdas, sizeof(double), eik_g->triM_2D->nPoints, fp);
  fclose(fp);
}
