/*
This is the neighbor struct. Since not all nodes have the same amount of neighbors,
we have to build a struct so that it mimics this "2d array with varying sizes of rows"
*/

#include "neighbors.h"

#include <stdio.h>
#include <stdlib.h>

struct neighborsR {
    int len;
    int* neis_i;
};

void neighbors_init(neighborsRS *neighbors, char const *pathNeighbors) {
    // here we read the txt file, which is the output from meshpy and save this to 
    // our struct neighborsRS called neighbors in this particular case

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
    // In this method given a line of text read from a file with nNeighRow 
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