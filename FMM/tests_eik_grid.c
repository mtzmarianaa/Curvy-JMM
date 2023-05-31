/*
TESTS FOR THE GRID
*/
#include "eik_grid.h"
#include "triMesh_2D.h"
#include "priority_queue.h"

#include <math.h>
#include <stdio.h>


int main()
{
    // we first test the ugly way of initializing the eik_grid struct, by building everything
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING FROM A TRIANGULAR MESH FOR A SIMPLE SQUARE \n\n\n\n");
    // triMesh_2Ds *triM_2D;
    // triMesh_2Dalloc(&triM_2D);
    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_IncidentFaces.txt";
    pathFaces =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_Faces.txt";
    pathIndexRegions =  "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H0_1/H0_1_FacesLabel.txt";

    int *start;
    int nStart, s;

    s = 0;
    start = &s;
    nStart = 1;


    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nStart, pathPoints, pathNeighbors, pathIncidentFaces, pathFaces, pathIndexRegions);
    printGeneralInfo(eik_g1);


    printf("\n\nInitializing the points within a radious of 3.2 from the starting point (0,-2)");

    double rBall = 3.2;

    initializePointsNear(eik_g1, rBall);

    printGeneralInfo(eik_g1);
    
    

    eik_grid_dealloc(&eik_g1);
}