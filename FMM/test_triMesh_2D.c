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
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/IncidentFaces.txt";
    pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/BoundaryPoints.txt";
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Facets.txt";
    pathFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Faces.txt";
    pathIndexRegions =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/FacesLabel.txt";
    triMesh2_init_from_meshpy(triMesh, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    printf("GENERAL INFO \n\n");
    printGeneralInfoMesh(triMesh);

    printf("\n\n\n\nALL INFO \n\n");
    printEverythingInMesh(triMesh);

    printf("\n\n Testing the face index stuff \n\n");

    printf( "\n Should be 1 %d \n", regionBetweenTwoPoints(triMesh, 4, 28)  );

    printf("\n Should be 6 %d \n", faceBetween3Points(triMesh, 5, 10, 12));



}