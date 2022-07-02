/*
This are the things we can do with coordinates (idk if this is the way to go)
*/
#include "coord.h"
#include "files_methods.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

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

void print_generalCoord(coordS *coordinates){
    printf("Number of coordinates: %d \n", coordinates->nPoints);
}

void coord_initFromFile(coordS *coordinates, char const *pathPoints) {
  // This method is going to store what is found un the pathPoints text file 
    FILE *fp; // the file were going to read
    char * line = NULL; // each line of the file
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    double row[2]; // temporary row
    row[0] = 0;
    row[1] = 0;
    int nPoints;

    nPoints = numLinesInFile(pathPoints); // we have the number of points in the mesh
    coordinates->nPoints = nPoints;
    coordinates->x = malloc(nPoints*sizeof(double));
    coordinates->y = malloc(nPoints*sizeof(double));

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
    fp = fopen(pathPoints, "r");
    while ((read = getline(&line, &len, fp)) != -1) {
    //printf("Iteration %d \n", i);
    //printf("%s", line);
    //printf("\n");
    separateARowDb(line, 2, row);
    //printf("First element found: %lf\n", row[0]);
    //printf("Second element found: %lf\n", row[1]);
    coordinates->x[i] = row[0];
    coordinates->y[i] = row[1];
    //printf("\n\n");
    i ++;
    }
    fclose(fp);
    if (line)
        free(line);
}
