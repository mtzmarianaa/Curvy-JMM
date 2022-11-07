#include "eik_grid.h"

#include <stdio.h>

void FMM_2D( eik_gridS *eik_g, double rBall){
    // we first add directly the points that are close
    initializePointsNear(eik_g, rBall);
    // then we can start marching
    //int currentMinInd;
    while(nStillInQueue(eik_g) != 0){
        // while we still have something in the queue
        //currentMinInd = currentMinIndex(eik_g);
        // printf("\n\n\n New iteration \n\n\n");
        // printf("There are still %d in the queue \n", nStillInQueue(eik_g));
        // printf("The minimum index is at %d\n", currentMinInd);
        popAddNeighbors(eik_g);
        // printGeneralInfo(eik_g);
    }
    // printGeneralInfo(eik_g);
}
