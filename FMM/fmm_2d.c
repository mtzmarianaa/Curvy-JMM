#include "eik_grid.h"

#include <stdio.h>

void FMM_2D( eik_gridS *eik_g ){
    int currentMinInd;
    printf("\nStarting iteration\n");
    while(nStillInQueue(eik_g) != 0){
        // while we still have something in the queue
        currentMinInd = currentMinIndex(eik_g);
        printf("There are still %d in the queue \n", nStillInQueue(eik_g));
        printf("The minimum index is at %d\n", currentMinInd);
        popAddNeighbors(eik_g); // first find minimum and add its neighbors if classified before as FAR
        update_afterAccepted(eik_g, currentMinInd);
    }
    printGeneralInfo(eik_g);
}