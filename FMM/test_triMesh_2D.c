#include "coord.h"
#include "triMesh_2D.h"
#include "faces.h"
#include "neighbors.h"
#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){
    double pi;
    pi = acos(-1.0);

    triMesh_2Ds *triMesh;
    triMesh_2Dalloc(&triMesh);
    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions;
    int trA, trB;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_IncidentFaces.txt";
    pathFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_Faces.txt";
    pathIndexRegions =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_FacesLabel.txt";
    triMesh2_init_from_meshpy(triMesh, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions);

    printf("GENERAL INFO \n\n");
    printGeneralInfoMesh(triMesh);

    printf("\n\n\n\nALL INFO \n\n");
    printEverythingInMesh(triMesh);

    printf("\nTesting the marching along triangles thing\n\n");

    infoTwoPartUpdate *infoOut;

    infoTwoPartUpdate_alloc(&infoOut);

    pointWhereRegionChanges(triMesh, 5, 33, 16, 1, infoOut);

    printf("\n\n\n");

    // printf("If there was a change in region %d \n", infoOut->xChange_ind);

    printf("\n  %fl  \n", infoOut->indexRef_01);
    printf("\n  %fl  \n", infoOut->indexRef_02);

    printf("The index of refraction changed from   %fl   to  %fl\n\n", infoOut->indexRef_01, infoOut->indexRef_02);

    printf("The angle from x0x1 to xHat is: %fl\n", infoOut->angle_xHat/pi);

    printf("The angle from x0x1 to xChange is: %fl\n", infoOut->angle_xChange/pi);

    printf("\n\n\n\n\nTrying the other way\n\n");

    pointWhereRegionChanges(triMesh, 5, 33, 16, 0, infoOut);

    printf("\n\n\n");

    // printf("If there was a change in region %d \n", infoOut->xChange_ind);

    printf("\n  %fl  \n", infoOut->indexRef_01);
    printf("\n  %fl  \n", infoOut->indexRef_02);

    printf("The index of refraction changed from   %fl   to  %fl\n\n", infoOut->indexRef_01, infoOut->indexRef_02);

    printf("The angle from x0x1 to xHat is: %fl\n", infoOut->angle_xHat/pi);

    printf("The angle from x0x1 to xChange is: %fl\n", infoOut->angle_xChange/pi);


}