/*
This are the things we can do with coordinates (idk if this is the way to go)
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "coord.h"



void coord_alloc(coordS **coordinates ) {
  *coordinates = malloc(sizeof(coordS));
  assert(*coordinates != NULL);
}

void coord_dealloc(coordS **coordinates) {
    free(coordinates);
    *coordinates = NULL;
}

void coord_init(coordS *coordinates, double *x, double *y, int nPoints) {
    coordinates->x = x;
    coordinates->y = y;
    coordinates->nPoints = nPoints;
}

void print_coord(coordS *coordinates) {
    printf("Number of coordinates: %d \n", coordinates->nPoints);
    for(int i = 0; i< coordinates->nPoints; i ++) {
        printf("x: %f,  y: %f \n", coordinates->x[i], coordinates->y[i]);
    }
}

void coord_initFromFile(coordS *coordinates, char const *pathPoints){
  // This method is going to store what is found un the pathPoints text file 
      FILE *fp; // the file were going to read
    char * line = NULL; // each line of the file
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathPoints, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }

    // scan such file

        while ((read = getline(&line, &len, fp)) != -1) {
        // printf("Retrieved line of length %d:\n", nCharInLine);
        // printf("Iteration %d \n", i);
        // printf("K %d \n", k);
        // printf("%s", line);
        // printf("Memory allocation of line: %p \n", &line);
        i ++;
        // printf("Temporary line: %s", line);
        // printf("Number of indices found in curreny line: %d :\n", nNei);
        // printf("\n");
        neighbors[k].neis_i = malloc(nNei*sizeof(int));
        separateARow(line, nNei, neighbors[k].neis_i);
        // printThisLinesNeighbors(neighbors[k].neis_i, nNei);
        // printf("\n\n");
        k++;
    }
    fclose(fp);
    if (line)
        free(line);
}
