/*
This is the path struct, for each node in the mesh this records the path
taken from x0 to the nodes
*/

#include "path.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


void path_alloc_n(pathS **path, int nPoints) 
{
    // memory allocates the number of points in the mesh
    *path = malloc(nPoints*sizeof(pathS));
    assert(*path != NULL);
}

void path_dealloc(pathS **path)
{
    free(*path);
    *path = NULL;
}

void grow_This_path( pathS *path, int indexHat ) 
{
  path[indexHat].maxSize *= 2;
  path[indexHat].individual_path_x = realloc( path[indexHat].individual_path_x, path[indexHat].maxSize*sizeof(double)  );
  path[indexHat].individual_path_y = realloc( path[indexHat].individual_path_y, path[indexHat].maxSize*sizeof(double)  );
}

void path_init(pathS *path, int nPoints, double *coord_x, double *coord_y) 
{
    for(int i = 0; i<nPoints; i++){
        path[i].maxSize = 10;
        path[i].len = 1;
        path[i].individual_path_x = malloc(10*sizeof(double));
        path[i].individual_path_x[0] = coord_x[i];
        path[i].individual_path_y = malloc(10*sizeof(double));
        path[i].individual_path_y[0] = coord_y[i];
    }
}

void insertOneToPath(pathS *path, int indexHat, double xlam[2])
{
    // we want to add a new point to the path but in the second-last index becuse the last index should hold the coordinates of 
    // the indexed point in question (because the path should end at such point)
    int last;
    double this_x, this_y;
    if( path[indexHat].len == path[indexHat].maxSize ){
        grow_This_path( path, indexHat );
    }
    last = path[indexHat].len - 1 ; // get the length of the current path computed
    this_x = path[indexHat].individual_path_x[last]; // x coordinate of the point in question
    this_y = path[indexHat].individual_path_y[last]; // y coordinate of the point in question
    path[indexHat].individual_path_x[last] = xlam[0]; // add the new xlambda[0] to the path
    path[indexHat].individual_path_y[last] = xlam[1]; // add the new xlambda[1] to the path
    path[indexHat].individual_path_x[last + 1] = this_x; // we have to have the coordinates of the point in question at the last index
    path[indexHat].individual_path_y[last + 1] = this_y;
    path[indexHat].len ++;
}

void updatePath(pathS *path, int indexHat, double xlam[2])
{
    // if we have a new xlambda such that the eikonal is minimized in this iteration then we must update the path taken
    int last;
    last = path[indexHat].len -1;
    path[indexHat].individual_path_x[last -1] = xlam[0];
    path[indexHat].individual_path_y[last -1] = xlam[1];
}

void addAllPathFromOneNode(pathS *path, int indexHat, int xNewAccepted) 
{
    // After adding the neighbors of a newly accepted node we must also add 
    // the path that goes from x0 all the way up to the newly accepted node
    int nInPath;
    double xlam[2];
    nInPath = path[xNewAccepted].len;
    for (int i=0; i<nInPath; i++){
        // get the coordinates of the ith point in the path from x0 to the newly accepted node
        xlam[0] = path[xNewAccepted].individual_path_x[i];
        xlam[1] = path[xNewAccepted].individual_path_y[i];
        // add those coordinates to the path of the neighbor
        insertOneToPath(path, indexHat, xlam);
    }
}

void printPaths(pathS *path, int nPoints) 
{
    for(int i = 0; i<nPoints; i++) {
        printf("\nFor index: %d. Length of this path: %d. Maximum length allowed in this path: %d\n", i, path[i].len, path[i].maxSize);
        printf("The path taken was:");
        for(int j=0; j<path[i].len; j++){
            printf("     (%fl, %fl)     |", path[i].individual_path_x[j], path[i].individual_path_y[j]);
        }
    }
}