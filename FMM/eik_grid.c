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

struct eik_grid {
  int *start; // the index of the point that is the source (could be multiple, that's why its a pointer)
  int nStart; // number of points in start
  triMesh_2Ds *triM_2D; // mesh structure that includes coordinates, neighbors, incident faces, etc.
  double *eik_vals; // the current Eikonal values for each indexed point in the mesh
  double (*eik_grad)[2]; // this is a pointer to a list of the gradients of the eikonal
  p_queue *p_queueG; // priority queue struct
  int *current_states; // 0 far, 1 trial, 2 valid
  pathS *pathsTaken; // path taken from x0 to each of the nodes
  int currentUpdate[2]; // 0 is xlambda, 1 is xhat (currently)
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
  int *current_states;
  eik_vals = malloc(triM_2D->nPoints*sizeof(double)); 
  current_states = malloc(triM_2D->nPoints*sizeof(int));
  eik_grad = malloc(2*triM_2D->nPoints*sizeof(double)); // each gradient has two coordinates (for each point)
  for(int i = 0; i<triM_2D->nPoints; i++){
    eik_vals[i] = INFINITY; // set them all to infinity
    current_states[i] = 0; // set them all to far
    eik_grad[i][0] = 0; // initialize all the gradients to zero
    eik_grad[i][1] = 0;
  }
  eik_g->eik_vals = eik_vals;
  eik_g->current_states = current_states;
  eik_g->eik_grad = eik_grad;

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

  // We also need to initialize the paths
  pathS *pathsTaken;
  path_alloc_n(&pathsTaken, triM_2D->nPoints);
  path_init(pathsTaken, triM_2D->nPoints);
  eik_g->pathsTaken = pathsTaken;
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
    printf("Index   %d    ||  Coordinates:   (%fl, %fl)    ||  Eikonal value:   %fl    ||  L2:   %fl    ||  Current state:   %d ||  Coordinates eik gradient:   (%fl, %fl) \n", i, x[0], x[1]  , eik_g->eik_vals[i], l2norm(x) , eik_g->current_states[i], eik_g->eik_grad[i][0], eik_g->eik_grad[i][1]);
  }
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
  printf("\n Region %d \n", *regionIndex);
  *That1 = eik_g->eik_vals[indexFrom] + s_function_threeSections(x0, *regionIndex)*dist; // 
}

void twoPointUpdate_eikValue(eik_gridS *eik_g, int x0_ind, int x1_ind, int xHat_ind, double xlam[2], double *That2, int *regionIndex){
  // this is where we use the optimization problem (we find lambda and then update it)
  double lambda_opt, lambda0, lambda1, T0, T1, tol;
  double x0[2], x1[2], xHat[2];
  int maxIter, faceBetweenPoints;
  lambda0 = 0.0;
  lambda1 = 1.0;
  maxIter = 25;
  tol = 0.001; // ask if these parameters are ok
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
  // printf("Lambda found %fl\n", lambda_opt);
  // save the xlam 
  xlam[0] = (1-lambda_opt)*x0[0] + lambda_opt*x1[0];
  xlam[1] = (1-lambda_opt)*x0[1] + lambda_opt*x1[1];
  // get the possible eikonal value for this two point update
  *That2 = eikApproxLin(T1, T0, lambda_opt, x0, x1, xHat, *regionIndex);
  // printf("Eikonal value before %fl\n", get_valueAtIndex(eik_g->p_queueG, xHat_ind));
  // printf("Eikonal value with two point update %fl\n", That2);
}


// once accepted a node (i.e. we deleated the root from the priority queue)
// we need to add its un taged neighbors to the priority queue
// WE JUST ADD THEM IF THEY ARE CURRENTLY FAR
// ADD TO QUEUE != UPDATE (which might include a 2 node update)
void addNeighbors_fromAccepted(eik_gridS *eik_g, int index_accepted){
  // we iterate through its neighbors and those with current state = 0 we add them to the priority queue
  int neighborsIndex, nNei, regionIndex;
  int *listNeighbors;
  double norm_div, That1, temp_substraction[2], xhat[2], xlam[2];
  temp_substraction[0] = 0;
  temp_substraction[1] = 0;
  xlam[0] = eik_g->triM_2D->points->x[index_accepted];
  xlam[1] = eik_g->triM_2D->points->y[index_accepted];
  nNei = eik_g->triM_2D->neighbors[index_accepted].len;
  for(int i = 0; i<nNei; i++){
    neighborsIndex = eik_g->triM_2D->neighbors[index_accepted].neis_i[i]; // a neighbor
    if(eik_g->current_states[neighborsIndex] == 0) {
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
    }
  }
}

// after accepting the neighboring nodes and putting them into the priority queue (their updates were all computed with
// either a one point update (if they were previously far) or stayed with their pre computed eikonal value)
// now we need to update those nodes that can be updated via a two point update
void update_afterAccepted(eik_gridS *eik_g, int index_accepted){
  // we need to first find the incident faces to the index_accepted point, then we iterate through those
  // faces to find the 2 other points that belong to the same face, if one of them is set as valid then we might consider performing
  // a two point update
  int nFaces, faceIndex, k, regionIndex;
  int neis[2];
  double twoPointVal;
  double norm_div, temp_substraction[2], xhat[2], xlam[2];
  temp_substraction[0] = 0;
  temp_substraction[1] = 0;
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
      xhat[0] = eik_g->triM_2D->points->x[neis[1]];
      xhat[1] = eik_g->triM_2D->points->y[neis[1]];
      twoPointUpdate_eikValue(eik_g, index_accepted, neis[0], neis[1], xlam, &twoPointVal, &regionIndex);
      update(eik_g->p_queueG, twoPointVal, neis[1]); // this function takes care, if its smaller it will update the new value
      vec2_substraction(xhat, xlam, temp_substraction);
      norm_div = l2norm(temp_substraction);
      eik_g->eik_grad[neis[1]][0] = s_function_threeSections(xhat, regionIndex)*temp_substraction[0]/norm_div;
      eik_g->eik_grad[neis[1]][1] = s_function_threeSections(xhat, regionIndex)*temp_substraction[1]/norm_div;
    }
    if( eik_g->current_states[neis[0]] == 1 & eik_g->current_states[neis[1]] == 2  ) {
      // neis[0] is trial, neis[1] is valid
      // two point update considered
      xhat[0] = eik_g->triM_2D->points->x[neis[0]];
      xhat[1] = eik_g->triM_2D->points->y[neis[0]];
      twoPointUpdate_eikValue(eik_g, index_accepted, neis[1], neis[0], xlam, &twoPointVal, &regionIndex);
      printf("\nValue from the two point update %fl \n", twoPointVal);
      printf("\nCoordinates of xlam (%fl, %fl) \n", xlam[0], xlam[1]);
      update(eik_g->p_queueG, twoPointVal, neis[0]);
      vec2_substraction(xhat, xlam, temp_substraction);
      norm_div = l2norm(temp_substraction);
      eik_g->eik_grad[neis[0]][0] = s_function_threeSections(xhat, regionIndex)*temp_substraction[0]/norm_div;
      eik_g->eik_grad[neis[0]][1] = s_function_threeSections(xhat, regionIndex)*temp_substraction[1]/norm_div;
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


void savePathsTaken(eik_gridS *eik_g, const char *pathFile){
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->triM_2D->nPoints; ++i){
      fwrite(eik_g->pathsTaken[i].individual_path, sizeof(double), eik_g->pathsTaken[i].len, fp);
    }

    fclose(fp);
}