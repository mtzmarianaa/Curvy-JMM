/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"

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

int *getCoordinatesFromIndex(eik_gridS *eik_g, int index)
{
  int coord[2];
  coord[0] = index%eik_g->M;
  coord[1] = index/eik_g->M;
  return coord;
}

int getIndexFromCoordinates(eik_gridS *eik_g, int coord[2])
{
  return eik_g->M*coord[1] + coord[0];
}

int *neighboursBool(eik_gridS *eik_g, int index)
{
  // returns an array of size 4 with 0 and 1. Position 0 is the southern neighbour, 1 is the western neighbour,
  // 2 is the eastern neighbour, 3 is the northern neighbour. If set to 0 then no neighbour found, if set to 1 neighbour found
  int neighbours_found[4];
  if( index > eik_g->N )
  {
    neighbours_found[0] = 1;
  }
  else
  {
    neighbours_found[0] = 0;
  }

  if( index%eik_g->N != 0 )
  {
    neighbours_found[1] = 1;
  }
  else
  {
    neighbours_found[1] = 0;
  }

  if( index%eik_g->N != eik_g->N -1 )
  {
    neighbours_found[2] = 1;
  }
  else
  {
    neighbours_found[2] = 0;
  }

  if( index/eik_g->N < eik_g->M - 1 )
  {
    neighbours_found[3] = 1;
  }
  else
  {
    neighbours_found[3] = 0;
  }

  return neighbours_found;
}



