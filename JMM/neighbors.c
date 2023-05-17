#include "neighbors.h"

/*
This is the neighbor struct. Since not all nodes have the same amount of neighbors,
we have to build a struct so that it mimics this "2d array with varying sizes of rows"
*/

#include "feature_test.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <assert.h>

#include "files_methods.h"

void neighborsRSalloc(neighborsRS **neighbors) {
    // memory allocates one thing of type neighborsRS
    *neighbors = malloc(sizeof(neighborsRS));
    assert(*neighbors != NULL);
}

void neighborsRSalloc_n(neighborsRS **neighbors, int N) {
    // memory allocates N things of type neighborsRS
    *neighbors = malloc(N*sizeof(neighborsRS));
    assert(*neighbors != NULL);
}

void neighborsRSdealloc(neighborsRS **neighbors){
    free(*neighbors);
    *neighbors = NULL;
}

void neighbors_init(neighborsRS *neighbors, char const *pathNeighbors, int N) {
    // here we read the txt file, which is the output from meshpy and save this to
    // our struct neighborsRS called neighbors in this particular case N are the
    // number of lines that the file has
    FILE *fp; // the file were going to read
    char * line = NULL; // each line of the file
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    int nNei;
    int i = 0;
    int k = 0;

    printf("\ntrying to open to file\n");

    fp = fopen(pathNeighbors, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file");
        exit(EXIT_FAILURE);
    }
    else{
        printf("\nFile successfully found\n");
    }
    // If we can indeed open such file we read it line by line because each line has different amount of elements
    // we need the number of lines in the file first thing
    // nPoints = numLinesInFile(pathNeighbors);
    // assert(nPoints == N); //they should match, if not then either the file is wrong or something is off


    fp = fopen(pathNeighbors, "r");
    // now we read each line and get the indices for each neighbor
    while ((read = getline(&line, &len, fp)) != -1) {
        int nCharInLine = (int) read -1;
        // printf("Retrieved line of length %d:\n", nCharInLine);
        // printf("Iteration %d \n", i);
        // printf("K %d \n", k);
        // printf("%s", line);
        // printf("Memory allocation of line: %p \n", &line);
        i ++;
        // printf("Temporary line: %s", line);
        nNei = numberElementsInRow(line, nCharInLine);
        neighbors[k].len = nNei;
        // printf("Number of indices found in curreny line: %d :\n", nNei);
        // printf("\n");
        neighbors[k].neis_i = malloc(nNei*sizeof(int));
        separateARowInt(line, nNei, neighbors[k].neis_i);
        // printThisLinesNeighbors(neighbors[k].neis_i, nNei);
        // printf("\n\n");
        k++;
    }
    fclose(fp);
    if (line)
        free(line);

}


void printThisLinesNeighbors(int *neighborsRow, int SizeRow) {
    printf("Number of neighbors: %d \n ", SizeRow);
    for (int i = 0; i<SizeRow; i++){
        printf("|  %d  |", neighborsRow[i]);
    }
    printf("\n \n");
}

void printAllNeighbors(neighborsRS *neighbors, int N) {
    for (int p=0; p<N; p++){
        int currentN;
        int *Neighs;
        currentN = neighbors[p].len;
        Neighs = neighbors[p].neis_i;
        printf("Neighbors for point indexed: %d \n", p);
        printf("Has %d neighbors: \n", currentN);
        for (int i= 0;i<currentN; i++){
            printf("|          %d          |", Neighs[i]);
        }
        printf("\n\n");
    }
}
