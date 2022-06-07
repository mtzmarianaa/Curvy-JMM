/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"
#include "SoSFunction.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

struct eik_grid {
  double x0;
  double y0;
  int start[2];
  int M;
  int N;
  double h;
  double *x_linspace; // x coordinates of the nodes (length N)
  double *y_linspace; // y coordinates of the nodes (length M)
  double *eik_gridVals;
  p_queue *p_queueG; // priority queue struct
  int *current_states;
} ;

void eik_grid_alloc(eik_gridS **eik_g ) {
  *eik_g = malloc(sizeof(eik_gridS));
  assert(*eik_g != NULL);
}

void eik_grid_dealloc(eik_gridS **eik_g ) {
  free(*eik_g);
  *eik_g = NULL;
}

void eik_grid_init( eik_gridS *eik_g, double x_m, double y_m, int start[2], int m, int n, double H ) 
{
  int i;
  eik_g->start[0] = start[0];
  eik_g->start[1] = start[1];
  eik_g->M = m;
  eik_g->N = n;
  eik_g->h = H;
  // Initialize the values of the x and y coordinates
  double x_lin[m], y_lin[n];
  for(i = 0; i<m; i++)
  {
    x_lin[i] = x_m + i*H;
  }
  for(i = 0; i<m; i++)
  { 
    y_lin[i] = y_m + i*H;
  }
  eik_g->x_linspace = malloc(n*sizeof(double));
  eik_g->y_linspace = malloc(m*sizeof(double));
  eik_g->x_linspace = x_lin;
  eik_g->y_linspace = y_lin;
  eik_g->x0 = x_lin[start[0]];
  eik_g->y0 = y_lin[start[1]];
  // Initialize the current values of all the grid nodes as infinity
  eik_g->eik_gridVals = malloc(m*n*sizeof(double));
  // Initialize the priority queue struct
  p_queue *p_queueGs;
  priority_queue_alloc(&p_queueGs);
  eik_g->p_queueG = p_queueGs;
  // And the array of current states
  eik_g->current_states = malloc(m*n*sizeof(int));
  for(i = 0; i<m*n; i++){
    eik_g->current_states[i] = 0; // all nodes are classified as far
    eik_g->eik_gridVals[i] = INFINITY; // all nodes have eikonal value set to infinity
  }
  assert(&eik_g != NULL); // eik_g should not be null
}

void print_eikonal_grid(eik_gridS *eik_g)
{
  int i;
  int j;
  for(i = eik_g->N - 1; i>-1; i --)
  {
    for(j = 0; j<eik_g->M; j++)
    {
      printf("   %fl   ", eik_g->eik_gridVals[i*eik_g->N + j] );
    }
    printf("\n");
  }
  printf("\n");
}

void print_currentStates(eik_gridS *eik_g)
{
  int i;
  int j;
  for(i = eik_g->N - 1; i>-1; i --)
  {
    for(j = 0; j<eik_g->M; j++)
    {
      printf("   %d   ", eik_g->current_states[i*eik_g->N + j] );
    }
    printf("\n");
  }
  printf("\n");
}

void setState(eik_gridS *eik_g, int index, int new_state)
{
  eik_g->current_states[index] = new_state;
}

void setValue(eik_gridS *eik_g, int index, double new_value)
{
  eik_g->eik_gridVals[index] = new_value;
}

void add_toPriorityQueue(eik_gridS *eik_g, int index, double new_value)
{
  // NOT THE SAME AS UPDATE A NODE'S VALUE THAT IS CURRENTLY IN THE PRIORITY QUEUE
  // If we're goint to add something to the priority queue then this means that its current 
  // eikonal value is infinity, its current state is 0, and its proposed value is GREATER
  // than any other value inside the priority queue at the moment
  eik_g->current_states[index] = 1; // switch from far to close
  insert_end(eik_g->p_queueG, new_value, index); // hence we insert it at the end (faster)
}

double onePointUpdate( eik_gridS *eik_g, int index, int targetIndex )
{
  double x, y, update;
  x = eik_g->x_linspace( getXCoordFromIndex(targetIndex) );
  y = eik_g->y_linspace( getYCoordFromIndex(targetIndex)  );
  update = eik_g->eik_gridVals[index] + eik_g->h*s_function( x, y ) ;
  return update;
}

int getXCoordFromIndex(eik_gridS *eik_g, int index)
{
  return index%eik_g->M;
}

int getYCoordFromIndex(eik_gridS *eik_g, int index)
{
  return index/eik_g->M;
}

int getIndexFromCoordinates(eik_gridS *eik_g, int coord[2])
{
  return eik_g->M*coord[1] + coord[0];
}

int neighborSouthB(eik_gridS *eik_g, int index)
{
  // 1 if the current index has a southern neighbor, 0 if not
  int bool_south;
  if( index > eik_g->N )
  {
    bool_south = 1;
  }
  else
  {
    bool_south = 0;
  }
  return bool_south;
}

int indexNeighborSouth(eik_gridS *eik_g, int index)
{
  return index-eik_g->M;
}

int neighborWestB(eik_gridS *eik_g, int index)
{
  // 1 if the current index has a western neighbor, 0 if not
  int bool_west;
  if( index%eik_g->N != 0 )
  {
    bool_west = 1;
  }
  else
  {
    bool_west = 0;
  }
  return bool_west;
}

int indexNeighborWest(eik_gridS *eik_g, int index)
{
  return index-1;
}

int neighborEastB(eik_gridS *eik_g, int index)
{
  // 1 if the current index has a eastern neighbor, 0 if not
  int bool_east;
  if( index%eik_g->N != eik_g->N -1 )
  {
    bool_east = 1;
  }
  else
  {
    bool_east = 0;
  }
  return bool_east;
}

int indexNeighborWest(eik_gridS *eik_g, int index)
{
  return index+1;
}

int neighborNorthB(eik_gridS *eik_g, int index)
{
  // 1 if the current index has a northern neighbor, 0 if not
  int bool_north;
  if( index/eik_g->N < eik_g->M - 1 )
  {
    bool_north = 1;
  }
  else
  {
    bool_north = 0;
  }
  return bool_north;
}

int indexNeighborWest(eik_gridS *eik_g, int index)
{
  return index+eik_g->M;
}

void addNeighbors(eik_gridS *eik_g, int index)
{
  double update;
  // Given an index, add its (at most) 4 neighbors to the prioriry queue and update their current states
  // If it has a southern neighbor and is set as far, add it to the priority queue or update its current value
  if( neighborSouthB(eik_g, index) == 1 && eik_g->current_states[indexNeighborSouth(eik_g, index)] == 0 )
  {
    // southern neighbor wasn't in the priority queue
    update = s_function( getXCoordFromIndex(eik_g, index), getYCoordFromIndex(eik_g, index) );
    insert_end(eik_g->p_queueG, update, indexNeighborSouth(eik_g, index));
  }
  
  if( neighborWestB(eik_g, index) == 1 && eik_g->current_states[indexNeighborWest(eik_g, index)] == 0  )
  {
    // western neighbor wasn't in the priority queue
    update = 1;
  }

}
