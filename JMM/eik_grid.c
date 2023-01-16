/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"
#include "SoSFunction.h"
#include "opti_method.h"
#include "linAlg.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

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
  int (*parents_path)[3]; // this is a pointer to a list of the parents that are used to update that node (by indices)
  double *lambdas; // pointer to the lambdas used to update each node using their parents
  double *mus;
  int *current_states;
  eik_vals = malloc(triM_2D->nPoints*sizeof(double)); 
  current_states = malloc(triM_2D->nPoints*sizeof(int));
  eik_grad = malloc(2*triM_2D->nPoints*sizeof(double)); // each gradient has two coordinates (for each point)
  parents_path = malloc(3*triM_2D->nPoints*sizeof(int)); // each node has two parents by index
  lambdas = malloc(triM_2D->nPoints*sizeof(double)); // each node was updated with one lambda which is a double
  mus = malloc(triM_2D->nPoints*sizeof(double));
  for(int i = 0; i<triM_2D->nPoints; i++){
    eik_vals[i] = INFINITY; // set them all to infinity
    current_states[i] = 0; // set them all to far
    eik_grad[i][0] = 0; // initialize all the gradients to zero
    eik_grad[i][1] = 0;
    parents_path[i][0] = i; // initialize both parents as themselves
    parents_path[i][1] = i;
    parents_path[i][2] = i;
    lambdas[i] = 0; // initialize all lambdas to 0 (meaning that its parent is itself, useful for the starting points)
    mus[i] = 0;
  }
  eik_g->eik_vals = eik_vals;
  eik_g->current_states = current_states;
  eik_g->eik_grad = eik_grad;
  eik_g->parents_path = parents_path;
  eik_g->lambdas = lambdas;
  eik_g->mus = mus;

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
  triMesh_2Ds *triM_2D;
  triMesh_2Dalloc(&triM_2D);
  triMesh2_init_from_meshpy(triM_2D, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions, pathBoundaryTan, pathBoundaryChain);
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
    printf("Index   %d    ||  Coordinates:   (%fl, %fl)    ||  Eikonal:   %fl     ||  Current state:   %d ||  Eik gradient:   (%fl, %fl)||  Parents:   (%d, %d, %d)||  LamOpt:   %fl ||  MuOpt:   %fl \n", i, x[0], x[1]  , eik_g->eik_vals[i] , eik_g->current_states[i], eik_g->eik_grad[i][0], eik_g->eik_grad[i][1], eik_g->parents_path[i][0], eik_g->parents_path[i][1], eik_g->parents_path[i][2], eik_g->lambdas[i], eik_g->mus[i]);
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
  *That2 = eikApproxLin(T0, T1, *lambda, x0, x1, xHat, indexRef);

}

void twoStepUpdate(double x0[2], double x1[2], double x2[2], double xHat[2], double T0, double T1, double indexRef_01, double indexRef_02, double *That_step2, double *lambda, double *mu){
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
  *lambda = optimizers[0];
  *mu = optimizers[1];
}

void approximateEikonalGradient(double x0[2], double x1[2], double xHat[2], double parameterization, double indexRefraction, double grad[2]) {
  // given two points in a base with its corresponding parameterization we approximate the gradient of the Eikonal
  double MinParametrization, xBasePart1[2], xBasePart2[2], xBase[2], normBase, direction[2];
  MinParametrization = 1 - parameterization;
  scalar_times_2vec( MinParametrization, x0, xBasePart1 );
  scalar_times_2vec( parameterization, x1, xBasePart2 );
  vec2_addition( xBasePart1, xBasePart2, xBase );
  vec2_subtraction(xHat, xBase, direction);
  normBase = l2norm(direction);
  grad[0] = indexRefraction*direction[0]/normBase;
  grad[1] = indexRefraction*direction[1]/normBase;
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
    xStart[0] = eik_g->triM_2D->points->x[ indexStart ];
    xStart[1] = eik_g->triM_2D->points->y[ indexStart];
    initialIndexRefraction = s_function_threeSections(xStart, eik_g->triM_2D->indexRegions[ eik_g->triM_2D->incidentFaces[indexStart].neis_i[0] ]);
    for(int i = 0; i<eik_g->triM_2D->nPoints; i++){
      xCurrent[0] = eik_g->triM_2D->points->x[i];
      xCurrent[1] = eik_g->triM_2D->points->y[i];
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

void updateCurrentValues(eik_gridS *eik_g, int indexToBeUpdated, int parent1, int parent2, double param, double TFound, double indexRefraction) {
  // this void manages the updates, it doesn't compute them. Once computed the necessary updated (in another funcion)
  // updateCurrentValues changes everything inside (updates the current eikonal value, its parents, the paramenter used, etc
  double grad[2], x0[2], x1[2], xHat[2];
  // if it was previously far then we add this to the queue directly, if it was trial we just update the queue
  if(eik_g->current_states[indexToBeUpdated] == 0){
    eik_g->current_states[indexToBeUpdated] = 1;
    insert(eik_g->p_queueG, TFound, indexToBeUpdated);
  }
  else if (eik_g->current_states[indexToBeUpdated] == 1)
  {
    update(eik_g->p_queueG, TFound, indexToBeUpdated);
  }
  eik_g->current_states[indexToBeUpdated] = 1;
  eik_g->lambdas[indexToBeUpdated] = param;
  eik_g->mus[indexToBeUpdated] = 0.0;
  eik_g->eik_vals[indexToBeUpdated] = TFound;
  eik_g->parents_path[indexToBeUpdated][0] = parent1;
  eik_g->parents_path[indexToBeUpdated][1] = parent2;
  eik_g->parents_path[indexToBeUpdated][2] = parent1;
  // calculate the gradient
  x0[0] = eik_g->triM_2D->points->x[parent1];
  x0[1] = eik_g->triM_2D->points->y[parent1];
  x1[0] = eik_g->triM_2D->points->x[parent2];
  x1[1] = eik_g->triM_2D->points->y[parent2];
  xHat[0] = eik_g->triM_2D->points->x[indexToBeUpdated];
  xHat[1] = eik_g->triM_2D->points->y[indexToBeUpdated];
  approximateEikonalGradient(x0, x1, xHat, param, indexRefraction, grad);
  eik_g->eik_grad[indexToBeUpdated][0] = grad[0];
  eik_g->eik_grad[indexToBeUpdated][1] = grad[1];

}


void updateCurrentValues3(eik_gridS *eik_g, int indexToBeUpdated, int parent0, int parent1, int parent2, double lambda, double mu, double TFound, double indexRefraction) {
  // this void manages the updates, it doesn't compute them. Once computed the necessary updated (in another funcion)
  // updateCurrentValues changes everything inside (updates the current eikonal value, its parents, the paramenter used, etc
  double grad[2], x1[2], x2[2], xHat[2];
  // if it was previously far then we add this to the queue directly, if it was trial we just update the queue
  if(eik_g->current_states[indexToBeUpdated] == 0){
    eik_g->current_states[indexToBeUpdated] = 1;
    insert(eik_g->p_queueG, TFound, indexToBeUpdated);
  }
  else if (eik_g->current_states[indexToBeUpdated] == 1)
  {
    update(eik_g->p_queueG, TFound, indexToBeUpdated);
  }
  eik_g->current_states[indexToBeUpdated] = 1;
  eik_g->lambdas[indexToBeUpdated] = lambda;
  eik_g->mus[indexToBeUpdated] = mu;
  eik_g->eik_vals[indexToBeUpdated] = TFound;
  eik_g->parents_path[indexToBeUpdated][0] = parent0;
  eik_g->parents_path[indexToBeUpdated][1] = parent1;
  eik_g->parents_path[indexToBeUpdated][2] = parent2;
  // calculate the gradient
  x1[0] = eik_g->triM_2D->points->x[parent1];
  x1[1] = eik_g->triM_2D->points->y[parent1];
  x2[0] = eik_g->triM_2D->points->x[parent2];
  x2[1] = eik_g->triM_2D->points->y[parent2];
  xHat[0] = eik_g->triM_2D->points->x[indexToBeUpdated];
  xHat[1] = eik_g->triM_2D->points->y[indexToBeUpdated];
  approximateEikonalGradient(x1, x2, xHat, mu, indexRefraction, grad);
  eik_g->eik_grad[indexToBeUpdated][0] = grad[0];
  eik_g->eik_grad[indexToBeUpdated][1] = grad[1];

}


void addNeighbors_fromAccepted(eik_gridS *eik_g, int indexAccepted) {
  // from the point indexAccepted which was recently set to valid we update its neighbors that are not set to valid
  // one point updates are skipped (since we have initialized points near the source we don't need one point updates
  // anymore. This function manages all the cases (if its a simple update or a two step update,
  int nNeis, x1_ind, xHat_ind;
  double x0[2], x1[2], x2[2], xHat[2], pi, currentTHat, param, T0, T1, lambda, mu;
  pi = acos(-1.0);
  nNeis = eik_g->triM_2D->neighbors[indexAccepted].len; // to know the amount of neighbors we might update
  x0[0] = eik_g->triM_2D->points->x[indexAccepted];
  x0[1] = eik_g->triM_2D->points->y[indexAccepted];
  T0 = eik_g->eik_vals[indexAccepted];
  // printf("\n\n\nAdding neighbors from index accepted %d, with coordinates (%fl,   %fl)\n\n", indexAccepted, x0[0], x0[1]);
  for( int i = 0; i < nNeis; i++ ){
    // iterate on the possible xHats
    xHat_ind = eik_g->triM_2D->neighbors[indexAccepted].neis_i[i];
    if( eik_g->current_states[xHat_ind] != 2 ){
      // we can only update those neighbors that are not valid at the moment
      xHat[0] = eik_g->triM_2D->points->x[xHat_ind];
      xHat[1] = eik_g->triM_2D->points->y[xHat_ind];
      // printf("The coordinates of xHat: %d    are    (%fl, %fl)\n", xHat_ind, xHat[0], xHat[1]);
      for( int j = 0; j < nNeis; j++ ){
        x1_ind = eik_g->triM_2D->neighbors[indexAccepted].neis_i[j];
        // printf("\nFrom x0 (index accepted): %d we are considering x1: %d   to update xHat:  %d\n\n", indexAccepted, x1_ind, xHat_ind);
        T1 = eik_g->eik_vals[x1_ind];
        // iterate on the possible bases 
        if( i != j & eik_g->current_states[x1_ind] == 2 ){
          // the base should include another point different from both x0 and xHat but x1 SHOULD BE VALID
          x1[0] =  eik_g->triM_2D->points->x[x1_ind];
          x1[1] =  eik_g->triM_2D->points->y[x1_ind];
          // first we need to know if its going to be a simple update or a two step update
          infoTwoPartUpdate *infoOut;
          infoTwoPartUpdate_alloc(&infoOut);

          // FIRST DIRECTION
          pointWhereRegionChanges(eik_g->triM_2D, indexAccepted, x1_ind, xHat_ind, 0, infoOut);
          // two options, either it changes region or it doesn't but we much watch out for the angles here
          if( infoOut->xChange_ind == -1 & infoOut->angle_xHat <= pi ){
            // this means that there is no change in direction + we can build a straight line from x0x1 to xHat
            // printf("\n\nThere is no change in index of refraction, trying a simple update\n");
            simple_Update(x0, x1, xHat, T0, T1, infoOut->indexRef_01, &currentTHat, &lambda);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              // printf("We found a better value with lamdba: %fl\n", param);
              updateCurrentValues(eik_g, xHat_ind, indexAccepted, x1_ind, lambda, currentTHat, infoOut->indexRef_01);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }
          if( infoOut->xChange_ind != -1 & infoOut->angle_xChange <= pi & fabs(infoOut->angle_xChange - infoOut->angle_xHat) <= pi ){
            // this means that there is a change in direction + we can build two straight lines, one from x0x1 to x0x2 and the other one from x0x2 to x0xHat
            // printf("\n\nThere IS a change in the index of refraction, trying a two part update\n");
            x2[0] = eik_g->triM_2D->points->x[infoOut->xChange_ind];
            x2[1] = eik_g->triM_2D->points->y[infoOut->xChange_ind];
            twoStepUpdate(x0, x1, x2, xHat, T0, T1, infoOut->indexRef_01, infoOut->indexRef_02, &currentTHat, &lambda, &mu);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              // printf("We found a better value with mu: %fl\n", param);
              updateCurrentValues3(eik_g, xHat_ind, indexAccepted, x1_ind, infoOut->xChange_ind, lambda, mu, currentTHat, infoOut->indexRef_02);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }

          // GO ON THE OTHER DIRECTION
          pointWhereRegionChanges(eik_g->triM_2D, indexAccepted, x1_ind, xHat_ind, 1, infoOut);
          // two options, either it changes region or it doesn't but we much watch out for the angles here
          if( infoOut->xChange_ind == -1 & infoOut->angle_xHat <= pi ){
            // this means that there is no change in direction + we can build a straight line from x0x1 to xHat
            // printf("\n\nThere is no change in index of refraction, trying a simple update\n");
            simple_Update(x0, x1, xHat, T0, T1, infoOut->indexRef_01, &currentTHat, &lambda);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              // printf("We found a better value with lamdba: %fl\n", param);
              updateCurrentValues(eik_g, xHat_ind, indexAccepted, x1_ind, lambda, currentTHat, infoOut->indexRef_01);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }
          if( infoOut->xChange_ind != -1 & infoOut->angle_xChange <= pi & fabs(infoOut->angle_xChange - infoOut->angle_xHat) <= pi ){
            // this means that there is a change in direction + we can build two straight lines, one from x0x1 to x0x2 and the other one from x0x2 to x0xHat
            // printf("\n\nThere IS a change in the index of refraction, trying a two part update\n");
            x2[0] = eik_g->triM_2D->points->x[infoOut->xChange_ind];
            x2[1] = eik_g->triM_2D->points->y[infoOut->xChange_ind];
            twoStepUpdate(x0, x1, x2, xHat, T0, T1, infoOut->indexRef_01, infoOut->indexRef_02, &currentTHat, &lambda, &mu);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              // printf("We found a better value with mu: %fl\n", param);
              updateCurrentValues3(eik_g, xHat_ind, indexAccepted, x1_ind, infoOut->xChange_ind, lambda, mu, currentTHat, infoOut->indexRef_02);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }
          infoTwoPartUpdate_dalloc(&infoOut); // deallocate this struct
        }
      }
    }
  }
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
  // printf("The coordinates of the current minimum value just accepted are: (%fl,%fl)\n", eik_g->triM_2D->points->x[minIndex], eik_g->triM_2D->points->y[minIndex]);
  // printf("Its neighbors are: \n");
  // nNeighs = eik_g->triM_2D->neighbors[minIndex].len; // get the number of neighbors of the minimum index
  // for(int i=0; i<nNeighs; i++){
  //   printf("|     %d     |", eik_g->triM_2D->neighbors[minIndex].neis_i[i] );
  // }
  // printf("\n");
  // for(int i=0; i<nNeighs; i++){
  //   printf("|     (%fl, %fl)     |", eik_g->triM_2D->points->x[eik_g->triM_2D->neighbors[minIndex].neis_i[i]], eik_g->triM_2D->points->y[eik_g->triM_2D->neighbors[minIndex].neis_i[i]] );
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
  fwrite(eik_g->eik_vals, sizeof(double), eik_g->triM_2D->nPoints, fp);
  fclose(fp);
}

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile) {
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->triM_2D->nPoints; ++i){
      fwrite(eik_g->eik_grad[i], sizeof(double), 2, fp);
    }

    fclose(fp);
}

void saveComputedParents(eik_gridS *eik_g, const char *pathFile) {
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->triM_2D->nPoints; ++i){
      fwrite(eik_g->parents_path[i], sizeof(int), 2, fp);
    }

    fclose(fp);
}

void saveComputedLambdas(eik_gridS *eik_g, const char *pathFile) {
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->lambdas, sizeof(double), eik_g->triM_2D->nPoints, fp);
  fclose(fp);
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
///// INTERPOLATING T(LAMBDA) WITH CUBIC HERMITE INTERPOLATION

void simple_UpdateCubic(eik_gridS *eik_g, double x0[2], double x1[2], double xHat[2], double T0,
			double T1, double grad0[2], double grad1[2],
			int indexAccepted, int xHat_ind, double *indexRef, double *That2, double *lambda) {
  // this is a simple two point update, meaning that there are no change in regions considered
  // uses optimization to find xLam, which is on the line defined from x0 to x1 parametrized by lambda
  // T(xLam) is approximated using cubic hermite polynomial

  // first we find the optimal lambda
  double tol, lambda0, lambda1;
  int maxIter;
  tol = 0.00001;
  maxIter = 30;
  lambda0 = 0;
  lambda1 = 1;

  *lambda = secantCubic_2D(lambda0, lambda1, T0, T1, grad0, grad1, x0, x1, xHat, tol, maxIter, *indexRef); // optimal lambda found
  if(*lambda == 0){
    *indexRef = s_function_threeSections(x0, regionBetweenTwoPoints(eik_g->triM_2D, indexAccepted, xHat_ind));
  }
  // compute the approximation to That2
  *That2 = eikApproxCubic(T0, T1, grad0, grad1, *lambda, x0, x1, xHat, *indexRef);

}

void twoStepUpdateCubic(eik_gridS *eik_g, double x0[2], double x1[2], double x2[2], double xHat[2], double T0, double T1,
			double grad0[2], double grad1[2], int indexAccepted, int xHat_ind,
			double *indexRef_01, double *indexRef_02, double *That_step2, double *lambda, double *mu){
  // a two step update, there is a change of the index of refraction from indexRef_01 to indexRef_02 from the triangles
  // x2 x0 x1 to the triangle xHat x0 x2
  double tol, optimizers[2];
  int maxIter;
  tol = 0.00001;
  maxIter = 30;
  // we use the projected gradient descent method to get both lambda and mu optimal
  projectedGradientDescentCubic(optimizers, T0, T1, grad0, grad1, x0, x1, x2, xHat, tol, maxIter, *indexRef_01, *indexRef_02);
  // use them to get the minimum path and the approximated eikonal USING CUBIC INTERPOLATION
  // notice that if both mu and lambda are zero then it means that we have an update via the edge
  // in this particular case we need to consider the smallest index of refraction between the two areas
  if( optimizers[0] == 0 & optimizers[1] == 0){
    *indexRef_01 = s_function_threeSections(x0, regionBetweenTwoPoints(eik_g->triM_2D, indexAccepted, xHat_ind));
    *indexRef_02 = s_function_threeSections(x0, regionBetweenTwoPoints(eik_g->triM_2D, indexAccepted, xHat_ind));
  }
  *That_step2 = eikApproxCubic_2Regions(T0, T1, grad0, grad1, optimizers[0], optimizers[1], x0, x1, x2, xHat, *indexRef_01, *indexRef_02);
  // and extract the optimal mu from the optimizers found (we're going to use this to approximate the gradient)
  *lambda = optimizers[0];
  *mu = optimizers[1];
}



void addNeighbors_fromAcceptedCubic(eik_gridS *eik_g, int indexAccepted) {
  // from the point indexAccepted which was recently set to valid we update its neighbors that are not set to valid
  // one point updates are skipped (since we have initialized points near the source we don't need one point updates
  // anymore. This function manages all the cases (if its a simple update or a two step update,
  int nNeis, x1_ind, xHat_ind;
  double x0[2], x1[2], x2[2], xHat[2], pi, currentTHat, lambda, mu, T0, T1, grad0[2], grad1[2];
  double indexRef_01, indexRef_02;
  pi = acos(-1.0);
  nNeis = eik_g->triM_2D->neighbors[indexAccepted].len; // to know the amount of neighbors we might update
  x0[0] = eik_g->triM_2D->points->x[indexAccepted];
  x0[1] = eik_g->triM_2D->points->y[indexAccepted];
  // since we are interpolating T(xlambda) with a cubir polynomial we also need the
  // computed gradient of T at x0
  grad0[0] = eik_g->eik_grad[indexAccepted][0];
  grad0[1] = eik_g->eik_grad[indexAccepted][1];
  T0 = eik_g->eik_vals[indexAccepted];
  // printf("\n\n\nAdding neighbors from index accepted %d, with coordinates (%fl,   %fl)\n\n", indexAccepted, x0[0], x0[1]);
  for( int i = 0; i < nNeis; i++ ){
    // iterate on the possible xHats
    xHat_ind = eik_g->triM_2D->neighbors[indexAccepted].neis_i[i];
    if( eik_g->current_states[xHat_ind] != 2 ){
      // we can only update those neighbors that are not valid at the moment
      xHat[0] = eik_g->triM_2D->points->x[xHat_ind];
      xHat[1] = eik_g->triM_2D->points->y[xHat_ind];
      // printf("The coordinates of xHat: %d    are    (%fl, %fl)\n", xHat_ind, xHat[0], xHat[1]);
      for( int j = 0; j < nNeis; j++ ){
        x1_ind = eik_g->triM_2D->neighbors[indexAccepted].neis_i[j];
        // printf("\nFrom x0 (index accepted): %d we are considering x1: %d   to update xHat:  %d\n\n", indexAccepted, x1_ind, xHat_ind);
        T1 = eik_g->eik_vals[x1_ind];
        // iterate on the possible bases 
        if( i != j & eik_g->current_states[x1_ind] == 2 ){
          // the base should include another point different from both x0 and xHat but x1 SHOULD BE VALID
          x1[0] =  eik_g->triM_2D->points->x[x1_ind];
          x1[1] =  eik_g->triM_2D->points->y[x1_ind];
	  grad1[0] = eik_g->eik_grad[x1_ind][0];
	  grad1[1] = eik_g->eik_grad[x1_ind][1];
          // first we need to know if its going to be a simple update or a two step update
          infoTwoPartUpdate *infoOut;
          infoTwoPartUpdate_alloc(&infoOut);

          // FIRST DIRECTION
          pointWhereRegionChanges(eik_g->triM_2D, indexAccepted, x1_ind, xHat_ind, 0, infoOut);
          // two options, either it changes region or it doesn't but we much watch out for the angles here
          if( infoOut->xChange_ind == -1 & infoOut->angle_xHat <= pi ){
            // this means that there is no change in direction + we can build a straight line from x0x1 to xHat
            // printf("\n\nThere is no change in index of refraction, trying a simple update\n");
            simple_UpdateCubic(eik_g, x0, x1, xHat, T0, T1, grad0, grad1, indexAccepted, xHat_ind, &infoOut->indexRef_01, &currentTHat, &lambda);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              updateCurrentValues(eik_g, xHat_ind, indexAccepted, x1_ind, lambda, currentTHat, infoOut->indexRef_01);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }
          if( infoOut->xChange_ind != -1 & infoOut->angle_xChange <= pi & fabs(infoOut->angle_xChange - infoOut->angle_xHat) <= pi ){
            // this means that there is a change in direction + we can build two straight lines, one from x0x1 to x0x2 and the other one from x0x2 to x0xHat
            // printf("\n\nThere IS a change in the index of refraction, trying a two part update\n");
            x2[0] = eik_g->triM_2D->points->x[infoOut->xChange_ind];
            x2[1] = eik_g->triM_2D->points->y[infoOut->xChange_ind];
            twoStepUpdateCubic(eik_g, x0, x1, x2, xHat, T0, T1, grad0, grad1, indexAccepted, xHat_ind, &infoOut->indexRef_01, &infoOut->indexRef_02, &currentTHat, &lambda, &mu);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              // printf("We found a better value with mu: %fl\n", param);
              updateCurrentValues3(eik_g, xHat_ind, indexAccepted, x1_ind, infoOut->xChange_ind, lambda, mu, currentTHat, infoOut->indexRef_02);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }

          // GO ON THE OTHER DIRECTION
          pointWhereRegionChanges(eik_g->triM_2D, indexAccepted, x1_ind, xHat_ind, 1, infoOut);
          // two options, either it changes region or it doesn't but we much watch out for the angles here
          if( infoOut->xChange_ind == -1 & infoOut->angle_xHat <= pi ){
            // this means that there is no change in direction + we can build a straight line from x0x1 to xHat
            // printf("\n\nThere is no change in index of refraction, trying a simple update\n");
            simple_UpdateCubic(eik_g, x0, x1, xHat, T0, T1, grad0, grad1, indexAccepted, xHat_ind, &infoOut->indexRef_01, &currentTHat, &lambda);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              // printf("We found a better value with lamdba: %fl\n", param);
              updateCurrentValues(eik_g, xHat_ind, indexAccepted, x1_ind, lambda, currentTHat, infoOut->indexRef_01);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }
          if( infoOut->xChange_ind != -1 & infoOut->angle_xChange <= pi & fabs(infoOut->angle_xChange - infoOut->angle_xHat) <= pi ){
            // this means that there is a change in direction + we can build two straight lines, one from x0x1 to x0x2 and the other one from x0x2 to x0xHat
            // printf("\n\nThere IS a change in the index of refraction, trying a two part update\n");
            x2[0] = eik_g->triM_2D->points->x[infoOut->xChange_ind];
            x2[1] = eik_g->triM_2D->points->y[infoOut->xChange_ind];
            twoStepUpdateCubic(eik_g, x0, x1, x2, xHat, T0, T1, grad0, grad1, indexAccepted, xHat_ind, &infoOut->indexRef_01, &infoOut->indexRef_02, &currentTHat, &lambda, &mu);
            if(currentTHat < eik_g->eik_vals[xHat_ind]){
              // we've found a better value
              // printf("We found a better value with mu: %fl\n", param);
              updateCurrentValues3(eik_g, xHat_ind, indexAccepted, x1_ind, infoOut->xChange_ind, lambda, mu, currentTHat, infoOut->indexRef_02);
              // printf("The new value (lambda) is: %fl\n", eik_g->lambdas[xHat_ind]);
            }
          }
          infoTwoPartUpdate_dalloc(&infoOut); // deallocate this struct
        }
      }
    }
  }
}

void popAddNeighborsCubic(eik_gridS *eik_g){
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
  addNeighbors_fromAcceptedCubic(eik_g, minIndex); // add neighbors from the recently accepted index
  // printf("The coordinates of the current minimum value just accepted are: (%fl,%fl)\n", eik_g->triM_2D->points->x[minIndex], eik_g->triM_2D->points->y[minIndex]);
  // printf("Its neighbors are: \n");
  // nNeighs = eik_g->triM_2D->neighbors[minIndex].len; // get the number of neighbors of the minimum index
  // for(int i=0; i<nNeighs; i++){
  //   printf("|     %d     |", eik_g->triM_2D->neighbors[minIndex].neis_i[i] );
  // }
  // printf("\n");
  // for(int i=0; i<nNeighs; i++){
  //   printf("|     (%fl, %fl)     |", eik_g->triM_2D->points->x[eik_g->triM_2D->neighbors[minIndex].neis_i[i]], eik_g->triM_2D->points->y[eik_g->triM_2D->neighbors[minIndex].neis_i[i]] );
  // }
  // printf("\n");
}

