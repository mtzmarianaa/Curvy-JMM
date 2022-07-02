/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"
#include "SoSFunction.h"
#include "linAlg.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct eik_grid {
  int *start; // the index of the point that is the source (could be multiple, that's why its a pointer)
  int nStart; // number of points in start
  triMesh_2Ds *triM_2D; // mesh structure that includes coordinates, neighbors, incident faces, etc.
  double *eik_vals; // the current Eikonal values for each indexed point in the mesh
  p_queue *p_queueG; // priority queue struct
  int *current_states; // 0 far, 1 trial, 2 valid
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
  int *current_states;
  eik_vals = malloc(triM_2D->nPoints*sizeof(double)); 
  current_states = malloc(triM_2D->nPoints*sizeof(int));
  for(int i = 0; i<triM_2D->nPoints; i++){
    eik_vals[i] = INFINITY; // set them all to infinity
    current_states[i] = 0; // set them all to far
  }
  eik_g->eik_vals = eik_vals;
  eik_g->current_states = current_states;

  // we initialize the priority queue, all the elements in start are going to be inserted and their current states set to 1
  // notice that we need to add both the index of the starting points AND their eikonal value (which is 0) to the p_queue struct
  p_queue *p_queueG;
  priority_queue_alloc(&p_queueG); // allocate
  priority_queue_init(p_queueG); // initiate
  for(int i = 0; i<nStart; i++){
    insert(p_queueG, 0, start[i]); // insert all the starting points with eikonal value 0
    eik_vals[start[i]] = 0; // the final eikonal values of the starting points are zero
  }
  eik_g->p_queueG = p_queueG;
  assert(&eik_g != NULL); // eik_g should not be null
}

void eik_grid_initFromFile(eik_gridS *eik_g, int *start, int nStart, char const *pathPoints, char const *pathNeighbors, char const *pathIncidentFaces, char const *pathBoundaryPoints, char const *pathFacets, char const *pathFaces) {
  // the only difference between this and the previous method is that in here we do need to initialize the Mesh structure
  triMesh_2Ds *triM_2D;
  triMesh_2Dalloc(&triM_2D);
  triMesh2_init_from_meshpy(triM_2D, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces);
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
}

// same principle as one point updates from square grid just that the length of the "arch" is not h, is not the same everywhere
double onePointUpdate_eikValue(eik_gridS *eik_g, int indexFrom, int indexTo){
  double That1, dist;
  double x1Minx0[2], x0[2], x1[0];
  // since we don't have a uniform square grid, we need to know how to handle these one point updates
  x0[0] = eik_g->triM_2D->points->x[indexFrom];
  x0[1] = eik_g->triM_2D->points->y[indexFrom];
  x1[0] = eik_g->triM_2D->points->x[indexTo];
  x1[1] = eik_g->triM_2D->points->y[indexTo];
  vec2_substraction(x0, x1, x1Minx0);
  dist = l2norm(x1Minx0);
  That1 = s_function(x0)*dist;
  return That1;
}

double twoPointUpdate_eikValue(eik_gridS *eik_g, int x0_ind, int x1_ind, int xHat_ind){
  // this is where we use the optimization problem (we find lambda and then update it)
}


// once accepted a node (i.e. we deleated the root from the priority queue)
// we need to add its un taged neighbors to the priority queue
// WE JUST ADD THEM IF THEY ARE CURRENTLY FAR
// ADD TO QUEUE != UPDATE (which might include a 2 node update)
void addNeighbors_fromAccepted(eik_gridS *eik_g, int index_accepted){
  // we iterate through its neighbors and those with current state = 0 we add them to the priority queue
  int neighborsIndex, nNei;
  int *listNeighbors;
  nNei = eik_g->triM_2D->neighbors[index_accepted].len;
  for(int i = 0; i<nNei; i++){
    neighborsIndex = eik_g->triM_2D->neighbors[index_accepted].neis_i[i]; // a neighbor
    if(eik_g->current_states[neighborsIndex] == 0) {
      insert(eik_g->p_queueG, onePointUpdate_eikValue(eik_g, index_accepted, neighborsIndex) , neighborsIndex); // insert this one point update to the priority queue
      eik_g->current_states[neighborsIndex] = 1; // set this to TRIAL
    }
  }
}