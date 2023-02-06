#include "triMesh_2D.h"
#include "eik_grid.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){
    double pi;
    pi = acos(-1.0);

    triMesh_2Ds *triM_2D;
    triMesh_2Dalloc(&triM_2D);
    eik_gridS *eik_g;
    eik_grid_alloc(&eik_g);
    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions, *path_BoundaryTan, *path_BoundaryChain;
    int trA, trB;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Neigh.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_IncidentFaces.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_FacesLabel.txt";
    path_BoundaryTan = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Boundary_tan.txt";
    path_BoundaryChain = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Boundary_chain.txt";

    int *start;
    int nstart, s;
    s = 0;
    start = &s;
    nstart = 1;
    eik_grid_initFromFile(eik_g, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions, path_BoundaryTan, path_BoundaryChain);


    
    triMesh2_init_from_meshpy(triM_2D, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions, path_BoundaryTan, path_BoundaryChain);

    printf("GENERAL INFO \n\n");
    printGeneralInfoMesh(triM_2D);

    printf("\n\n\n\nALL INFO \n\n");
    printEverythingInMesh(triM_2D);

    printf("\nTesting the marching along triangles thing\n\n");

    infoTriangleFan *infoOut;

    infoTriangleFan_alloc(&infoOut);

    pointWhereRegionChanges(triM_2D, eik_g->current_states, 115, 116, 114, 1, infoOut);

    printf("\n\n\n");

    // printf("If there was a change in region %d \n", infoOut->xChange_ind);

    printf("\n Number of changes in the index of refraction: %d\n", infoOut->nChanges);
    printf("If this update should be considered or not %d\n", infoOut->updatable);

    printf("\n  %fl  \n", infoOut->indexRef_01);
    printf("\n  %fl  \n", infoOut->indexRef_02);

    printf("The index of refraction changed from   %fl   to  %fl\n\n", infoOut->indexRef_01, infoOut->indexRef_02);

    printf("The angle from x0x1 to xHat is: %fl\n", infoOut->angle_xHat/pi);

    printf("The angle from x0x1 to xChange is: %fl\n", infoOut->angle_xChange/pi);

    printf("\n\n\n\n\nTrying the other way\n\n");

    pointWhereRegionChanges(triM_2D, eik_g->current_states, 115, 116, 114, 0, infoOut);

    printf("\n\n\n");

    // printf("If there was a change in region %d \n", infoOut->xChange_ind);

    printf("\n Number of changes in the index of refraction: %d\n", infoOut->nChanges);
    printf("If this update should be considered or not %d\n", infoOut->updatable);

    printf("\n  %fl  \n", infoOut->indexRef_01);
    printf("\n  %fl  \n", infoOut->indexRef_02);

    printf("The index of refraction changed from   %fl   to  %fl\n\n", infoOut->indexRef_01, infoOut->indexRef_02);

    printf("The angle from x0x1 to xHat is: %fl\n", infoOut->angle_xHat/pi);

    printf("The angle from x0x1 to xChange is: %fl\n", infoOut->angle_xChange/pi);


}
