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
  int southN, westN, eastN, northN, index;
  // no southern or western neighbors
  index = 0;
  southN = neighborSouthB(eikonal_g, index);
  westN = neighborWestB(eikonal_g, index);
  eastN = neighborEastB(eikonal_g, index);
  northN = neighborNorthB(eikonal_g, index);
  printf(" no southern or western neighbors \n");
  printf( "%d \n", southN );
  printf( "%d \n", westN );
  printf( "%d \n", eastN );
  printf( "%d \n", northN  );
  // no northern or eastern
  index = 63;
  southN = neighborSouthB(eikonal_g, index);
  westN = neighborWestB(eikonal_g, index);
  eastN = neighborEastB(eikonal_g, index);
  northN = neighborNorthB(eikonal_g, index);
  printf(" no northern or eastern \n");
  printf( "%d \n", southN );
  printf( "%d \n", westN );
  printf( "%d \n", eastN );
  printf( "%d \n", northN  );
  // all neighbours
  index = 10;
  southN = neighborSouthB(eikonal_g, index);
  westN = neighborWestB(eikonal_g, index);
  eastN = neighborEastB(eikonal_g, index);
  northN = neighborNorthB(eikonal_g, index);
  printf(" all neighbours \n");
  printf( "%d \n", southN );
  printf( "%d \n", westN );
  printf( "%d \n", eastN );
  printf( "%d \n", northN  );


  

}