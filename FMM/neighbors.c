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