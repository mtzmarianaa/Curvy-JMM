#pragma once

typedef struct eik_grid eik_gridS;

void eik_grid_alloc(eik_gridS **eik_g );

void eik_grid_dealloc(eik_gridS **eik_g );

void eik_grid_init( eik_gridS *eik_g, int *start, int nStart, triMesh_2Ds *triM_2D) ;

void printGeneralInfo(eik_gridS *eik_g);
