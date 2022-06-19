/*
This are the things we can do with coordinates (idk if this is the way to go)
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "coord.h"

struct coord {
  double *x;
  double *y;
  int nPoints;
} ;

void coord_alloc(coordS **coordinates ) {
  *coordinates = malloc(sizeof(coordS));
  assert(*coordinates != NULL);
}

void coord_dealloc(coordS **coordinates) {
    free(coordinates);
    *coordinates = NULL;
}

void coord_init(coordS *coordinates, double *x, double *y, int nPoints) {
    coordinates->x = x;
    coordinates->y = y;
    coordinates->nPoints = nPoints;
}

void print_coord(coordS *coordinates) {
    printf("Number of coordinates: %d \n", coordinates->nPoints);
    for(int i = 0; i< coordinates->nPoints; i ++) {
        printf("x: %f,  y: %f \n", coordinates->x[i], coordinates->y[i]);
    }
}
