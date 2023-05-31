/*
TESTS FOR THE GRID
*/
#include "eik_grid.h"
#include "fmm_2d.h"

#include <math.h>
#include <stdio.h>


int main()
{
    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions, *pathToSave;
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/MeshPoints_Sq.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Neigh_Sq.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/IncidentFaces_Sq.txt";
    pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/BoundaryPoints_Sq.txt";
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Facets_Sq.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Faces_Sq.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/FacesLabel_Sq.txt";

    pathToSave = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/ComputedValues.bin";
 
    int *start;
    int nStart, s;

    s = 562;
    start = &s;
    nStart = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TEST GEOMETRY (SNOWMAN) \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nStart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);
    printGeneralInfo(eik_g1);

    // we use the fmm olim method
    FMM_2D( eik_g1 );

    saveComputedValues(eik_g1, pathToSave);

    eik_grid_dealloc(&eik_g1);




}