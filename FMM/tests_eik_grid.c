/*
TESTS FOR THE GRID
*/
#include <math.h>
#include "eik_grid.h"
#include "priority_queue.h"

int main()
{
  
  eik_gridS *eikonal_g ;
  eik_grid_alloc(&eikonal_g );
  int start[2];
  start[0] = 1;
  start[1] = 2;
  eik_grid_init( eikonal_g, 0.0, 0.0, start, 8, 8, 0.1 );
  print_eikonal_grid(eikonal_g);
  print_currentStates(eikonal_g);
  setState(eikonal_g, start[0]*8 + start[1], 2);
  setValue(eikonal_g, start[0]*8 + start[1], 0);
  print_eikonal_grid(eikonal_g);
  print_currentStates(eikonal_g);


  

}