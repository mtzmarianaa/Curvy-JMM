/*
TESTS FOR THE PRIORITY QUEUE
*/
#include <stdio.h>
#include <stdlib.h>

#include "path.h"

int main()
{
  
   pathS *path_test;
   int nPoints = 15;

   path_alloc_n(&path_test, nPoints);

   path_init(path_test, nPoints);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 1 from node 0\n");

   insertAfterAcceptedToPath(path_test, 1, 0);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 2 from node 1\n");

   insertAfterAcceptedToPath(path_test, 2, 1);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 3 from node 2\n");

   insertAfterAcceptedToPath(path_test, 3, 2);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 4 from node 3\n");

   insertAfterAcceptedToPath(path_test, 4, 3);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 5 from node 4\n");

   insertAfterAcceptedToPath(path_test, 5, 4);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 6 from node 5\n");

   insertAfterAcceptedToPath(path_test, 6, 5);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 7 from node 6\n");

   insertAfterAcceptedToPath(path_test, 7, 6);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 8 from node 7\n");

   insertAfterAcceptedToPath(path_test, 8, 7);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 9 from node 8\n");

   insertAfterAcceptedToPath(path_test, 9, 8);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 10 from node 9\n");

   insertAfterAcceptedToPath(path_test, 10, 9);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 11 from node 8\n");

   insertAfterAcceptedToPath(path_test, 11, 8);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 12 from node 3\n");

   insertAfterAcceptedToPath(path_test, 12, 3);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 13 from node 0\n");

   insertAfterAcceptedToPath(path_test, 13, 0);

   printPaths(path_test, nPoints);

   printf("\n\n\n We accept node 14 from node 6\n");

   insertAfterAcceptedToPath(path_test, 14, 6);

   printPaths(path_test, nPoints);
   

}