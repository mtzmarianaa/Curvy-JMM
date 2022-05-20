/*
TESTS FOR THE GRID
*/
#include <math.h>
#include "eik_grid.c"

int main()
{
  
  eik_gridS *eikonal_g = malloc(sizeof(eik_gridS));
  int start[2];
  start[0] = 1;
  start[1] = 2;
  eik_grid_init( eikonal_g, 0.0, 0.0, start, 8, 8, 0.1 );
  print_eikonal_grid(eikonal_g);
  print_currentStates(eikonal_g);
  eikonal_g->current_states[ eikonal_g->start[0]*8 + eikonal_g->start[1] ] = 2;
  insert( eikonal_g->p_queueG, 0.0,  eikonal_g->start[0]*8 + eikonal_g->start[1]  ); // set the starting point's Eikonal value to 0
  deleteRoot(eikonal_g->p_queueG); // delete the root because we have set this point to valid
  eikonal_g->eik_gridVals[ eikonal_g->start[0]*8 + eikonal_g->start[1] ] = 0.0;
  print_eikonal_grid(eikonal_g);
  print_currentStates(eikonal_g);


  

}