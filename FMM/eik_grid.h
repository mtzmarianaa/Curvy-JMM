#pragma once

typedef struct eik_grid eik_gridS;

void eik_grid_alloc(eik_gridS *eik_g);

void eik_grid_dealloc( eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, double x_min, double y_min, int start[2], int M, int N, double h );

static void print_eikonal_grid(eik_gridS *eik_g);

static void print_currentStates(eik_gridS *eik_g);