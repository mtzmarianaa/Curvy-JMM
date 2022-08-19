#include "eik_grid.h"

#include <stdio.h>

void FMM_2D( eik_gridS *eik_g, int L){
    // L is going to determine weather we skip updates with nonzero Lagrange multiplier or not, it's a flag
    int currentMinInd;
    printf("\n\nStarting iteration\n\n");
    while(nStillInQueue(eik_g) != 0){
        // while we still have something in the queue
        currentMinInd = currentMinIndex(eik_g);
        //printf("There are still %d in the queue \n", nStillInQueue(eik_g));
        //printf("The minimum index is at %d\n", currentMinInd);
        popAddNeighbors(eik_g); // first find minimum and add its neighbors if classified before as FAR
        if (L == 1){
            // if this happens then after adding the neighbors we do a round of artificial triangle update
            updateWithArtificial(eik_g, currentMinInd);
        }
        else{
            update_afterAccepted(eik_g, currentMinInd);
        }
        printf("\n\n\n\n\n New ITERATION \n\n");
        printGeneralInfo(eik_g);
    }
    //printGeneralInfo(eik_g);
}