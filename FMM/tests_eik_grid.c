/*
TESTS FOR THE GRID
*/
#include <math.h>
#include <stdio.h>
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
  // testing all the functions
  int neighbours_found[4];
  neighbours_found = neighboursBool(eikonal_g, 5);
  printf( "%d /n", neighbours_found[0] );
  printf( "%d /n", neighbours_found[1] );
  printf( "%d /n", neighbours_found[2] );
  printf( "%d /n", neighbours_found[3] );


  

}