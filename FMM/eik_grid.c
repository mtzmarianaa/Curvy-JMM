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

static void initialize_queue(double *eik_queue, int *index_queue, int M, int N)
{
    eik_queue = malloc(M*N*sizeof(double));
    index_queue = malloc(M*N*sizeof(int));
}

