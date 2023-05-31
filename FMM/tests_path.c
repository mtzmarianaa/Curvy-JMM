/*
TESTS FOR THE PRIORITY QUEUE
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "path.h"

int main()
{
  
   pathS *path_test;
   double *coord_x, *coord_y, xlam[2];
   int nPoints = 4;
   double x[nPoints], y[nPoints];
   x[0] = 0.0;
   y[0] = 0.0;
   x[1] = 1.0;
   y[1] = 0.0;
   x[2] = 0.0;
   y[2] = 1.0;
   x[3] = 1.0;
   y[3] = 1.0;

   coord_x = x;
   coord_y = y;

   xlam[0] = -1.0;
   xlam[1] = -1.0;

   path_alloc_n(&path_test, nPoints);

   path_init(path_test, nPoints, x, y);

   printf("Initializing paths\n");

   printPaths(path_test, nPoints);

   insertOneToPath(path_test, 0, xlam);

   printf("\n\nAdding (-1,-1) to the path of the origin\n");

   printPaths(path_test, nPoints);

   xlam[0] = -sqrt(2)/2;
   xlam[1] = -sqrt(2)/2;

   updatePath(path_test, 0, xlam);

   printf("\n\nAdding updating the path of the origin with sqr(2)/2\n");

   printPaths(path_test, nPoints);

   addAllPathFromOneNode(path_test, 1, 0);

   printf("\n\nFrom the origin we add (1, 0)\n");

   printPaths(path_test, nPoints);
   

}