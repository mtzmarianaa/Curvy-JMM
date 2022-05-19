/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include <stdio.h>

typedef struct eik_grid {
  double x_min;
  double y_min;
  int start[2];
  int M;
  int N;
  double h;
  double *grid;
  double *eik_queue;
  int *index_queue;
  int *current_states;
} eik_gridS;

void eik_queue_alloc(double *eik_queue, int *index_queue)
{
    eik_queue = malloc(16*sizeof(double));
    index_queue = malloc(16*sizeof(int));
}

void eik_grid_alloc( eik_gridS **eik_g ) {
  *eik_g = malloc( sizeof(eik_gridS)  );
  assert(*eik_g != NULL); // eik_g should not be null
}

void eik_grid_dealloc( eik_gridS **eik_g ) {
  free( *eik_g );
  *eik_g = NULL;
}


