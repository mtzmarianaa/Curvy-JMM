#include "coord.h"
#include "triMesh_2D.h"
#include "facets.h"
#include "faces.h"
#include "neighbors.h"
#include "files_methods.h"

#include <stdio.h>
#include <stdlib.h>


int main(){

    triMesh_2Ds *triMesh;
    triMesh_2Dalloc(&triMesh);
    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions;
    int trA, trB;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_IncidentFaces.txt";
    pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_BoundaryPoints.txt";
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_Facets.txt";
    pathFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_Faces.txt";
    pathIndexRegions =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/H1_FacesLabel.txt";
    triMesh2_init_from_meshpy(triMesh, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    printf("GENERAL INFO \n\n");
    printGeneralInfoMesh(triMesh);

    printf("\n\n\n\nALL INFO \n\n");
    printEverythingInMesh(triMesh);

    printf("\n\n Testing the face index stuff \n\n");

    printf( "\n Should be 1 %d \n", regionBetweenTwoPoints(triMesh, 4, 28)  );

    printf("\n Should be 6 %d \n", faceBetween3Points(triMesh, 5, 10, 12));

    findTrATrB(triMesh, 5, 10, 14, &trA, &trB);

    printf("\n Triangle 1 should be 6 %d \n", trA);

    printf("\n Should be 10, 12, 5: %d | %d | %d", triMesh->faces->points[trA][0], triMesh->faces->points[trA][1], triMesh->faces->points[trA][2]);

    printf("\n Triangle 2 should be 13 %d \n", trB);

    printf("\n Should be 10, 12, 14: %d | %d | %d", triMesh->faces->points[trB][0], triMesh->faces->points[trB][1], triMesh->faces->points[trB][2]);

}