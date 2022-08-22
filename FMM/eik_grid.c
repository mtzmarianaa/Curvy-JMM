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

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathBoundaryPoints, char const *pathFacets, char const *pathFaces, char const *pathIndexRegions) {
  // the only difference between this and the previous method is that in here we do need to initialize the Mesh structure
  triMesh_2Ds *triM_2D;
  triMesh_2Dalloc(&triM_2D);
  triMesh2_init_from_meshpy(triM_2D, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);
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

// same principle as one point updates from square grid just that the length of the "arch" is not h, is not the same everywhere
void onePointUpdate_eikValue(eik_gridS *eik_g, int indexFrom, int indexTo, double *That1, int *regionIndex){
  double dist;
  double x1Minx0[2], x0[2], x1[0];
  // since we don't have a uniform square grid, we need to know how to handle these one point updates
  x0[0] = eik_g->triM_2D->points->x[indexFrom];
  x0[1] = eik_g->triM_2D->points->y[indexFrom];
  x1[0] = eik_g->triM_2D->points->x[indexTo];
  x1[1] = eik_g->triM_2D->points->y[indexTo];
  vec2_substraction(x0, x1, x1Minx0);
  dist = l2norm(x1Minx0);
  *regionIndex = regionBetweenTwoPoints(eik_g->triM_2D, indexFrom, indexTo);
  //printf("\n Region %d \n", *regionIndex);
  *That1 = eik_g->eik_vals[indexFrom] + s_function_threeSections(x0, *regionIndex)*dist; // 
}

void twoPointUpdate_eikValue(eik_gridS *eik_g, int x0_ind, int *x1_ind_original, int xHat_ind, double *lambda, double xlam[2], double *That2, int *regionIndex){
  // this is where we use the optimization problem (we find lambda and then update it)
  double lambda_opt, lambda0, lambda1, T0, T1, tol;
  double x0[2], x1[2], xHat[2], xlamNeigh[2];
  int maxIter, faceBetweenPoints, neighNeigh_ind, x1_ind;
  x1_ind = *x1_ind_original;
  lambda0 = 0.0;
  lambda1 = 1.0;
  maxIter = 25;
  tol = 0.0005;
  T0 = eik_g->eik_vals[x0_ind];
  T1 = eik_g->eik_vals[x1_ind];
  faceBetweenPoints = faceBetween3Points(eik_g->triM_2D, x0_ind, x1_ind, xHat_ind ); // get the face that is defined by these 3 points
  *regionIndex = eik_g->triM_2D->indexRegions[ faceBetweenPoints ]; // get the region where this face belongs to
  // get the coordinates of the points x0, x1, xHat
  x0[0] = eik_g->triM_2D->points->x[x0_ind];
  // printf("\n");
  // printf("\nx0 %fl\n", x0[0]);
  x0[1] = eik_g->triM_2D->points->y[x0_ind];
  // printf("\ny0 %fl\n", x0[1]);
  x1[0] = eik_g->triM_2D->points->x[x1_ind];
  // printf("\nx1 %fl\n", x1[0]);
  x1[1] = eik_g->triM_2D->points->y[x1_ind];
  // printf("\ny1 %fl\n", x1[1]);
  xHat[0] = eik_g->triM_2D->points->x[xHat_ind];
  // printf("\nxHat %fl\n", xHat[0]);
  xHat[1] = eik_g->triM_2D->points->y[xHat_ind];
  // printf("\nyHat %fl\n", xHat[1]);
  // compute the optimum lambda from  the linear model
  lambda_opt = secant_2D(lambda0, lambda1, T0, T1, x0, x1, xHat, tol, maxIter, *regionIndex);
  *lambda = lambda_opt;
  // printf("Lambda found %fl\n", lambda_opt);
  // save the xlam 
  xlam[0] = (1-lambda_opt)*x0[0] + lambda_opt*x1[0];
  xlam[1] = (1-lambda_opt)*x0[1] + lambda_opt*x1[1];
  // get the possible eikonal value for this two point update
  *That2 = eikApproxLin(T1, T0, lambda_opt, x0, x1, xHat, *regionIndex);
  // printf("Eikonal value before %fl\n", get_valueAtIndex(eik_g->p_queueG, xHat_ind));
  // printf("Eikonal value with two point update %fl\n", That2);
  
}

void lambdaWithArtificialTriangleAllowed(eik_gridS *eik_g, int x0_ind, int *x1_ind, int xHat_ind, double *lambda, double xlam[2], double *That2, int *regionIndex){
  // this is where it gets creative and allow two types of updates: with an artificial triangle (if both triangles involved 
  // belong to the same region or with a two piecewise linear update)
  // lambda0, lambda1 are the possible initial values for lambda T0 is the eikonal at x0, 
  // x0 is the newly accepted node, xhat is the node we want to update, x1_ind is the neighbor of x0 that is used for this weird two point update
  // we set the initial minimum value as the current value of xhat (which is just the one point update that was done from adding the neighbors of x0)
  
  // set the initial parameters
  double lambda0, lambda1, tol;
  int maxIter;
  lambda0 = 0;
  lambda1 = 1;
  tol = 0.0005;
  maxIter = 25;
  // printf("\n\n\n\nStarting to consider an artificial triangle update\n");
  int indexNeighNeigh, trA, trB; // indices of the neighbor of the neighbor (i.e. neighbor of x1)
  double xNeighNeigh[2], x0[2], xHat[2], TNeighNeigh, lambda_optNeighNeigh, ThatCurrent, T0;
  T0 = eik_g->eik_vals[x0_ind];
  *That2 = eik_g->eik_vals[xHat_ind]; // this is the value that we want to minimize using an artificial triangle
  *x1_ind = x0_ind; // and originally such value was just a one point update
  x0[0] = eik_g->triM_2D->points->x[x0_ind];
  x0[1] = eik_g->triM_2D->points->y[x0_ind]; 
  xlam[0] = x0[0];
  xlam[1] = x0[1]; 
  xHat[0] = eik_g->triM_2D->points->x[xHat_ind];
  xHat[1] = eik_g->triM_2D->points->y[xHat_ind]; 
  *lambda = 0; // because it was a one point update
  *regionIndex = regionBetweenTwoPoints(eik_g->triM_2D, x0_ind, xHat_ind); // the original one point update region is the region from this two points
  // first, we have to iterate through the neighbors of the newly accepted node that are set to valid
  for(int i = 0; i<eik_g->triM_2D->neighbors[x0_ind].len; i++ ){
    indexNeighNeigh = eik_g->triM_2D->neighbors[x0_ind].neis_i[i];
    if(eik_g->current_states[indexNeighNeigh] == 2){ // just look at the neighbors of x0 that are set to valid
      findTrATrB(eik_g->triM_2D, xHat_ind, x0_ind, indexNeighNeigh, &trA, &trB); // find the region(s) of the articial triangle
      if (trA != -1 & trB!= -1){
        // it means that we are considering triangles that are together, we don't want to consider triangles that don't share a side because that could imply that more than 2 regions
        // are involved and that gets messy ASK
        if(eik_g->triM_2D->indexRegions[trA] == eik_g->triM_2D->indexRegions[trB]){
          // if this happens then we are not changing regions, it is "easier"/less difficult because we can just use an artificial triangle
          // the corners of such artificial triangle are x0, indexNeighNeigh, xHat
          *regionIndex = eik_g->triM_2D->indexRegions[ trA ];
          // get the coordinates of the neighbor of the neighbor
          xNeighNeigh[0] = eik_g->triM_2D->points->x[indexNeighNeigh];
          xNeighNeigh[1] = eik_g->triM_2D->points->y[indexNeighNeigh];
          TNeighNeigh = eik_g->eik_vals[indexNeighNeigh];
          // compute the artificial triangle update
          lambda_optNeighNeigh = secant_2D(lambda0, lambda1, T0, TNeighNeigh, x0, xNeighNeigh, xHat, tol, maxIter, *regionIndex);
          ThatCurrent = eikApproxLin(TNeighNeigh, T0, lambda_optNeighNeigh, x0, xNeighNeigh, xHat, *regionIndex);
          // printf("\n Current ThatCurren: %fl for index %d,  comes from the parents (%d,%d)  that have valid state with eikonal (%fl, %fl)\n", ThatCurrent, xHat_ind, x0_ind, indexNeighNeigh, T0, TNeighNeigh);
          if( ThatCurrent < *That2){
            // printf("The current value %fl has been subsituted by this new value %fl\n", *That2, ThatCurrent);
            *That2 = ThatCurrent; // we found an update that is better
            *x1_ind = indexNeighNeigh; // we set this as the index neigh neigh (we need to save this to save the parents of xhat)
            *lambda = lambda_optNeighNeigh;
            // printf("\nThere was an artificial triangle used here\n");
            }
        }
        else{
          // if this happens then we ARE changing regions (horrible, we dont like this and we need a 2D optimization method gggg)
          // ASK SAM HOW TO SAVE THIS TWO PARAMETERS FOR THE PATH
          ThatCurrent = 6000;
          if( ThatCurrent < *That2){
            xlam[0] = 0;
            xlam[1] = 0;
            *x1_ind = -5;
            *That2 = 6000;
            }
          }
      }

    }
  }
}


// once accepted a node (i.e. we deleated the root from the priority queue)
// we need to add its un taged neighbors to the priority queue
// WE JUST ADD THEM IF THEY ARE CURRENTLY FAR
// ADD TO QUEUE != UPDATE (which might include a 2 node update)
void addNeighbors_fromAccepted(eik_gridS *eik_g, int index_accepted){
  // we iterate through its neighbors and those with current state = 0 we add them to the priority queue
  int neighborsIndex, nNei, regionIndex;
  double norm_div, That1, temp_substraction[2], xhat[2], xlam[2];
  temp_substraction[0] = 0;
  temp_substraction[1] = 0;
  xlam[0] = eik_g->triM_2D->points->x[index_accepted];
  xlam[1] = eik_g->triM_2D->points->y[index_accepted];
  nNei = eik_g->triM_2D->neighbors[index_accepted].len;
  for(int i = 0; i<nNei; i++){
    neighborsIndex = eik_g->triM_2D->neighbors[index_accepted].neis_i[i]; // a neighbor
    if(eik_g->current_states[neighborsIndex] == 0) { // meaning that we just add the neighbors which are currently set as far
      onePointUpdate_eikValue(eik_g, index_accepted, neighborsIndex, &That1, &regionIndex);
      insert(eik_g->p_queueG, That1 , neighborsIndex); // insert this one point update to the priority queue
      eik_g->current_states[neighborsIndex] = 1; // set this to TRIAL
      // now we add the approximated value of the gradient of the eikonal (this is kind of a trial of such gradient)
      xhat[0] = eik_g->triM_2D->points->x[neighborsIndex];
      xhat[1] = eik_g->triM_2D->points->y[neighborsIndex];
      vec2_substraction(xhat, xlam, temp_substraction);
      norm_div = l2norm(temp_substraction);
      eik_g->eik_grad[neighborsIndex][0] = s_function_threeSections(xhat, regionIndex)*temp_substraction[0]/norm_div;
      eik_g->eik_grad[neighborsIndex][1] = s_function_threeSections(xhat, regionIndex)*temp_substraction[1]/norm_div;
      eik_g->parents_path[neighborsIndex][0] = index_accepted; // it was updated directly from the newly accepted node
      // we don't update the lambdas because currently they're set to 0, i.e. the update is set directly from parents_path[neighborsIndex][0].
    }
  }
}

void updateWithArtificial(eik_gridS *eik_g, int index_accepted){
  // once we added the neighbors of index_accepted to the queue we can update them considering artificial triangles
  int neighborsIndex, nNei, regionIndex, x1_ind;
  double lambda, xlam[2], That2, xhat[2], temp_substraction[2], norm_div;
  nNei = eik_g->triM_2D->neighbors[index_accepted].len; // get the number of neighbors we might update
  for(int i = 0; i<nNei; i ++){
    neighborsIndex = eik_g->triM_2D->neighbors[index_accepted].neis_i[i]; // the current neighbor being considered
    if(eik_g->current_states[neighborsIndex] == 1) { // WE JUST UPDATE NODES THAT ARE CURRENTLY SET TO TRIAL
      lambdaWithArtificialTriangleAllowed(eik_g, index_accepted, &x1_ind, neighborsIndex, &lambda, xlam, &That2, &regionIndex);
      // we update the current eikonal values and add that to the queue
      update(eik_g->p_queueG, That2, neighborsIndex);
      if (That2 < eik_g->eik_vals[neighborsIndex]){
        eik_g->eik_vals[neighborsIndex] = That2;
        // we add the approximated value of the gradient of the eikonal
        xhat[0] = eik_g->triM_2D->points->x[neighborsIndex];
        xhat[1] = eik_g->triM_2D->points->y[neighborsIndex];
        vec2_substraction(xhat, xlam, temp_substraction);
        norm_div = l2norm(temp_substraction);
        eik_g->eik_grad[neighborsIndex][0] = s_function_threeSections(xhat, regionIndex)*temp_substraction[0]/norm_div;
        eik_g->eik_grad[neighborsIndex][1] = s_function_threeSections(xhat, regionIndex)*temp_substraction[1]/norm_div;
        eik_g->parents_path[neighborsIndex][0] = index_accepted;
        eik_g->parents_path[neighborsIndex][1] = x1_ind;
        eik_g->lambdas[neighborsIndex] = lambda;
      }
    }
  }
}

void update_step(eik_gridS *eik_g, int neighborValid, int neighborTrial, int index_accepted){
  // after accepting a node, if it has a valid neighbor and a trial neighbor we might perform an update in the 
  // current eikonal value of that trial neighbor using the newly accepted node and the valid neighbor (two point update)
  int regionIndex;
  double norm_div, twoPointVal, lambda_opt, temp_substraction[2], xlam[2], xhat[2];
  xhat[0] = eik_g->triM_2D->points->x[neighborTrial];
  xhat[1] = eik_g->triM_2D->points->y[neighborTrial];
  twoPointUpdate_eikValue(eik_g, index_accepted, &neighborValid, neighborTrial, &lambda_opt, xlam, &twoPointVal, &regionIndex);
  update(eik_g->p_queueG, twoPointVal, neighborTrial); // this function takes care, if its smaller it will update the new value
  if( twoPointVal < eik_g->eik_vals[neighborTrial] ){ // if we actually have a better update
      vec2_substraction(xhat, xlam, temp_substraction);
      norm_div = l2norm(temp_substraction);
      eik_g->eik_grad[neighborTrial][0] = s_function_threeSections(xhat, regionIndex)*temp_substraction[0]/norm_div;
      eik_g->eik_grad[neighborTrial][1] = s_function_threeSections(xhat, regionIndex)*temp_substraction[1]/norm_div;
      eik_g->parents_path[neighborTrial][0] = index_accepted;
      eik_g->parents_path[neighborTrial][1] = neighborValid;
      eik_g->lambdas[neighborTrial] = lambda_opt;
  }
}

// after accepting the neighboring nodes and putting them into the priority queue (their updates were all computed with
// either a one point update (if they were previously far) or stayed with their pre computed eikonal value)
// now we need to update those nodes that can be updated via a two point update
void update_afterAccepted(eik_gridS *eik_g, int index_accepted){
  // we need to first find the incident faces to the index_accepted point, then we iterate through those
  // faces to find the 2 other points that belong to the same face, if one of them is set as valid then we might consider performing
  // a two point update
  int nFaces, faceIndex, k, regionIndex, neighborValid, neighborTrial;
  int neis[2];
  double twoPointVal;
  nFaces = eik_g->triM_2D->incidentFaces[index_accepted].len; // get the number of incident faces to the newly accepted point
  // we iterate through those faces
  for (int i = 0; i<nFaces; i++){
    faceIndex = eik_g->triM_2D->incidentFaces[index_accepted].neis_i[i]; // i-th incident face 
    // now we get the two other points that belong to that face
    k = 0;
    for (int j = 0; j<3; j++) {
      if( eik_g->triM_2D->faces->points[faceIndex][j] != index_accepted ) { // we want the point in the face that are not index_accepted
        neis[k] = eik_g->triM_2D->faces->points[faceIndex][j];
        k++;
      }
    }
    // once we have the other two points that form part of that face (and should at least be trial)
    // if one of them is valid and the other one is trial we consider doing a two point update
    if( eik_g->current_states[neis[0]] == 2 &  eik_g->current_states[neis[1]] == 1 ) {
      // neis[0] is valid, neis[1] is trial
      // two point update considered
      neighborValid = neis[0];
      neighborTrial = neis[1];
      update_step(eik_g, neighborValid, neighborTrial, index_accepted);
    }
    if( eik_g->current_states[neis[0]] == 1 & eik_g->current_states[neis[1]] == 2  ) {
      // neis[0] is trial, neis[1] is valid
      // two point update considered
      neighborValid = neis[1];
      neighborTrial = neis[0];
      update_step(eik_g, neighborValid, neighborTrial, index_accepted);
    }
    // if both of them are trial it means that they were recently added as the neighbors of index_accepted, no two point update can be done in this case
  }
}

void popAddNeighbors(eik_gridS *eik_g){
  int minIndex = indexRoot(eik_g->p_queueG);
  double minEikVal = valueRoot(eik_g->p_queueG);
  deleteRoot(eik_g->p_queueG); // delete the root from the priority queue
  eik_g->current_states[minIndex] = 2; // set the newly accepted index to valid
  eik_g->eik_vals[minIndex] = minEikVal; // add the computed eikonal value to the list of eikonal values
  addNeighbors_fromAccepted(eik_g, minIndex); // add neighbors from the recently accepted index
}

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
