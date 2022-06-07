#pragma once

typedef struct eik_grid eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, double x_m, double y_m, int start[2], int m, int n, double H ) ;

void print_eikonal_grid(eik_gridS *eik_g);

void print_currentStates(eik_gridS *eik_g);

void setState(eik_gridS *eik_g, int index, int new_state);

void setValue(eik_gridS *eik_g, int index, double new_value);

void add_toPriorityQueue(eik_gridS *eik_g, int index, double new_value);

double onePointUpdate( eik_gridS *eik_g, int index, int targetIndex );

int getXCoordFromIndex(eik_gridS *eik_g, int index);

int getYCoordFromIndex(eik_gridS *eik_g, int index);

int getIndexFromCoordinates(eik_gridS *eik_g, int coord[2]);

int neighborSouthB(eik_gridS *eik_g, int index);

int neighborWestB(eik_gridS *eik_g, int index);

int neighborEastB(eik_gridS *eik_g, int index);

int neighborNorthB(eik_gridS *eik_g, int index);