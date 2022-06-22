#pragma once

typedef struct {
  double *x;
  double *y;
  int nPoints;
} coordS;

void coord_alloc(coordS **coordinates );

void coord_dealloc(coordS **coordinates);

void coord_init(coordS *coordinates, double *x, double *y, int nPoints);

void print_coord(coordS *coordinates);

void coord_initFromFile(coordS *coordinates, char const *pathPoints);