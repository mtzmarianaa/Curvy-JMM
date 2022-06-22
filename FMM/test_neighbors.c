#include "neighbors.h"
#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(){

    int N;
    char const *pathNeighbors;

    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Neigh.txt";
    N = numLinesInFile(pathNeighbors);

    neighborsRS *neighbors;

    neighborsRSalloc_n(&neighbors, N); // allocate those N lines of organized information


    printf("\n Size 1 %p", &neighbors);
    printf("\n Size 2%lu", sizeof(neighbors));
    

    neighbors_init(neighbors, pathNeighbors, N);

    printf("\n Size 1 %p", &neighbors);
    printf("\n Size 2%lu", sizeof(neighbors));
    printf("\n Size 3%lu", sizeof(*neighbors));

    printf("\n Printing all the neighbors found: \n");
    printAllNeighbors(neighbors, N);

}