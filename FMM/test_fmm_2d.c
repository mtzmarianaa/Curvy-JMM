/*
TESTS FOR THE GRID
*/
#include "eik_grid.h"
#include "fmm_2d.h"

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
    const char *pathPoints_sq, *pathNeighbors_sq, *pathIncidentFaces_sq, *pathBoundaryPoints_sq, *pathFacets_sq, *pathFaces_sq;
    pathPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/MeshPoints.txt";
    pathNeighbors_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Neigh.txt";
    pathIncidentFaces_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/IncidentFaces.txt";
    pathBoundaryPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/BoundaryPoints.txt";
    pathFacets_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Facets.txt";
    pathFaces_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Faces.txt";
 
    int *start;
    int nStart, s;

    s = 4;
    start = &s;
    nStart = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SIMPLE SQUARE \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nStart, pathPoints_sq, pathNeighbors_sq, pathIncidentFaces_sq, pathBoundaryPoints_sq, pathFacets_sq, pathFaces_sq);
    printGeneralInfo(eik_g1);

    // we use the fmm olim method
    FMM_2D( eik_g1 );

    eik_grid_dealloc(&eik_g1);
}