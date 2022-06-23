/*
This is a Nx2 int array type of structue. It is useful for "connecting indexed dots"
*/

#include "facets.h"
#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


void facets_alloc(facetsS **facets ) {
  *facets = malloc(sizeof(facetsS));
  assert(*facets != NULL);
}

void facets_dealloc(facetsS **facets) {
    free(facets);
    *facets = NULL;
}

void facets_init(facetsS *facets, int *from, int *to, int nFacets) {
    facets->from = from;
    facets->to = to;
    facets->nFacets = nFacets;
}

void print_facets(facetsS *facets) {
    printf("Number of facets: %d \n", facets->nFacets);
    for(int i = 0; i< facets->nFacets; i ++) {
        printf("From index: %d,  To index: %d \n", facets->from[i], facets->to[i]);
    }
}

void facets_initFromFile(facetsS *facets, char const *pathFacets) {
  // This method is going to store what is found un the pathFacets text file 
    FILE *fp; // the file were going to read
    char * line = NULL; // each line of the file
    int i = 0;
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    int *row; // temporary row of the file we're reading
    row = malloc(2*sizeof(int));
    row[0] = 0;
    row[1] = 0;
    int nFacets;

    nFacets = numLinesInFile(pathFacets); // we have the number of Facets in the file
    facets->nFacets = nFacets;
    facets->from = malloc(nFacets*sizeof(int));
    facets->to = malloc(nFacets*sizeof(int));

    printf("\ntrying to open to file\n");
    // Check if the file exists under that path
    fp = fopen(pathFacets, "r");
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
        // printf("Iteration %d \n", i);
        // printf("%s", line);
        i ++;
        // printf("\n");
        separateARowInt(line, 2, row);
        // printf("First element found: %d\n", row[0]);
        // printf("Second element found: %d\n", row[1]);
        facets->from[i] = row[0];
        facets->to[i] = row[1];
        // printf("\n\n");
    }
    fclose(fp);
    if (line)
        free(line);
}

