/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.c"

#include <stdio.h>
#include <math.h>

typedef struct eik_grid {
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
} eik_gridS;

void eik_grid_alloc(eik_gridS *eik_g)
{
  eik_g = malloc(sizeof(eik_gridS));
}

void eik_grid_dealloc( eik_gridS **eik_g ) {
  free( *eik_g );
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
  eik_g->x_linspace = x_lin;
  eik_g->y_linspace = y_lin;
  eik_g->x0 = x_lin[start[0]];
  eik_g->y0 = y_lin[start[1]];
  // Initialize the current values of all the grid nodes as infinity
  eik_g->eik_gridVals = malloc(m*n*sizeof(double));
  // Initialize the priority queue struct
  p_queue *p_queueGs = malloc(sizeof(p_queue));
  eik_g->p_queueG = p_queueGs;
  priority_queue_init(eik_g->p_queueG);
  // And the array of current states
  eik_g->current_states = malloc(m*n*sizeof(int));
  for(i = 0; i<m*n; i++){
    eik_g->current_states[i] = 0; // all nodes are classified as far
    eik_g->eik_gridVals[i] = INFINITY; // all nodes have eikonal value set to infinity
  }
  assert(&eik_g != NULL); // eik_g should not be null
}

static void print_eikonal_grid(eik_gridS *eik_g)
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

static void print_currentStates(eik_gridS *eik_g)
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



