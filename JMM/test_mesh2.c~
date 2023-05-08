#include "triMesh_2D.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){
    double pi;
    pi = acos(-1.0);

    triMesh_2Ds *triMesh;
    triMesh_2Dalloc(&triMesh);
    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions;
    const char *pathBoundaryTan, *pathBoundaryChain;
    int trA, trB;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_IncidentFaces.txt";
    pathFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Faces.txt";
    pathIndexRegions =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_FacesLabel.txt";
    pathBoundaryTan = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Boundary_tan.txt"; 
    pathBoundaryChain = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/TestBaseSnow/H0/H0_Boundary_chain.txt";
    triMesh2_init_from_meshpy(triMesh, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions, pathBoundaryTan, pathBoundaryChain);

    printf("GENERAL INFO \n\n");
    printGeneralInfoMesh(triMesh);

    printf("\n\n\n\nALL INFO \n\n");
    printEverythingInMesh(triMesh);


}
