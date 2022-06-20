/*
This is the neighbor struct. Since not all nodes have the same amount of neighbors,
we have to build a struct so that it mimics this "2d array with varying sizes of rows"
*/

#include "neighbors.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

struct neighborsR {
    int len;
    int* neis_i;
};

void neigborsRSalloc(neighborsRS **neighbors){
    *neighbors = malloc(sizeof(neighborsRS));
    assert(*neighbors != NULL);
}

void neighborsRSdealloc(neighborsRS **neighbors){
    free(*neighbors);
    *neighbors = NULL;
}

void neighbors_init(neighborsRS *neighbors, char const *pathNeighbors) {
    // here we read the txt file, which is the output from meshpy and save this to 
    // our struct neighborsRS called neighbors in this particular case
    FILE *fp; // the file were going to read
    char * line = NULL; // each line of the file
    size_t len = 0; // length of each line of the file
    ssize_t read; // reading each line in the file
    int nNei, nPoints;
    int i = 0;
    char c;
    int k = 0;

    fp = fopen(pathNeighbors, "r");
    // if the file doesnt exist or for some reason it can't be opened:
    if (fp == NULL) {
        printf("No such file");
        exit(EXIT_FAILURE);
    }
    // If we can indeed open such file we read it line by line because each line has different amount of elements
    // we need the number of lines in the file first thing
    nPoints = 0;
    for (c = getc(fp); c != EOF; c = getc(fp)){
        if (c == '\n'){
            nPoints ++;
        }
    }
    // we allocate memory for the amount of neighbors well need
    neigborsRSalloc(&neighbors);
    // now we read each line and get the indices for each neighbor
    while ((read = getline(&line, &len, fp)) != -1) {
        int nCharInLine = (int) read -1;
        printf("Retrieved line of length %d:\n", nCharInLine);
        printf("Iteration %d \n", i);
        printf("K %d \n", k);
        printf("%s", line);
        printf("Memory allocation of line: %p \n", &line);
        i ++;
        printf("Temporary line: %s", line);
        nNei = numberNeighborsFound(line, nCharInLine);
        neighbors[k].len = nNei;
        printf("Number of indices found in curreny line: %d :\n", nNei);
        printf("\n");
        int *neighborsRow = malloc(nNei*sizeof(int));
        separateARow(line, nNei, neighborsRow);
        neighbors[k].neis_i = neighborsRow;
        printThisLinesNeighbors(neighborsRow, nNei);
        printf("\n\n");
        k++;
    }
    fclose(fp);
    if (line)
        free(line);
    exit(EXIT_SUCCESS);
}

int numberNeighborsFound(char *line, int nCharInLine) {
    int i, count;
    count = 0;
    printf("%s  \n", line);
    for (i=0; i<nCharInLine; i++){
        if ((char) line[i] == ','){
            count+= 1;
        }
    }
    count += 1; // plus 1 because in theory we're counting the ,
    return count;
}

void separateARow(char *line, int nNei, int *neighborsRow) {
    // In this method given a line of text read from a file 
    // indices of neighbors we separate them and put it in neighborsRow.
    char *char_part;
    int int_part;
    int_part = strtol(line, &char_part, 10);
    neighborsRow[0] = int_part;
    char *end = char_part;
    char_part = char_part + 1;
    for (int i = 1; i<nNei; i++){
        int_part = strtol(char_part++, &end , 10);
        neighborsRow[i] = int_part;
        char_part = end + 1;
    }
}

void printThisLinesNeighbors(int *neighborsRow, int SizeRow) {
    for (int i = 0; i<SizeRow; i++){
        printf(" %d ", neighborsRow[i]);
    }
}

void printAllNeighbors(neighborsRS *neighbors) {
    int nPoints;
    nPoints = 100;
    for (int p=0; p<nPoints; p++){
        int currentN;
        int *Neighs;
        currentN = neighbors[p].len;
        Neighs = neighbors[p].neis_i;
        for (int i= 0;i<currentN; i++){
            printf("%d  ", Neighs[i]);
        }
        printf("\n");
    }
}