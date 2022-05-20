/*
TESTS FOR THE 2d FMM
*/
#include <math.h>
#include "fmm_2d.c"


int main()
{
  
  eik_gridS *eikonal_g = malloc(sizeof(eik_gridS));
  int start[2];
  start[0] = 1;
  start[1] = 2;
  eik_grid_init( eikonal_g, 0.0, 0.0, start, 8, 8, 1 );
  //print_eikonal_grid(eikonal_g);
  //print_currentStates(eikonal_g);
  FMM_2D( eikonal_g );
  printf("OK");
  print_eikonal_grid(eikonal_g);
  print_currentStates(eikonal_g);


  

}