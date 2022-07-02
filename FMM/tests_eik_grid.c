/*
TESTS FOR THE GRID
*/
#include "eik_grid.h"
#include "triMesh_2D.h"

#include <math.h>
#include <stdio.h>


int main()
{
    // we first test the ugly way of initializing the eik_grid struct, by building everything
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM A TRIANGULAR MESH FOR A SIMPLE SQUARE \n\n\n\n");
    triMesh_2Ds *triM_2D;
    triMesh_2Dalloc(&triM_2D);
    const char *pathPoints_sq, *pathNeighbors_sq, *pathIncidentFaces_sq, *pathBoundaryPoints_sq, *pathFacets_sq, *pathFaces_sq;
    pathPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/MeshPoints.txt";
    pathNeighbors_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Neigh.txt";
    pathIncidentFaces_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/IncidentFaces.txt";
    pathBoundaryPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/BoundaryPoints.txt";
    pathFacets_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Facets.txt";
    pathFaces_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Faces.txt";
    triMesh2_init_from_meshpy(triM_2D, pathPoints_sq, pathNeighbors_sq, pathIncidentFaces_sq, pathBoundaryPoints_sq, pathFacets_sq, pathFaces_sq);

    int *start;
    int nStart, s;

    s = 4;
    start = &s;
    nStart = 1;

    eik_gridS *eik_g;
    eik_grid_alloc(&eik_g);

    eik_grid_init( eik_g, start, nStart, triM_2D);

    printGeneralInfo(eik_g);

    triMesh_2Ddalloc(&triM_2D);
    eik_grid_dealloc(&eik_g);

    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SIMPLE SQUARE \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nStart, pathPoints_sq, pathNeighbors_sq, pathIncidentFaces_sq, pathBoundaryPoints_sq, pathFacets_sq, pathFaces_sq);
    printGeneralInfo(eik_g1);
    eik_grid_dealloc(&eik_g1);
    

}