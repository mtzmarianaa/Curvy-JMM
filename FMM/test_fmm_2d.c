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
    const char *pathPoints_sq, *pathNeighbors_sq, *pathIncidentFaces_sq, *pathBoundaryPoints_sq, *pathFacets_sq, *pathFaces_sq, *pathIndexRegions, *pathToSave, *pathSaveGradients;
    pathPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/MeshPoints.txt";
    pathNeighbors_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Neigh.txt";
    pathIncidentFaces_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/IncidentFaces.txt";
    pathBoundaryPoints_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/BoundaryPoints.txt";
    pathFacets_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Facets.txt";
    pathFaces_sq = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/FacesLabel.txt";

    pathToSave = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/ComputedValues.bin";
    pathSaveGradients = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquare/ComputedGradients.bin";
 
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
    eik_grid_initFromFile(eik_g1, start, nStart, pathPoints_sq, pathNeighbors_sq, pathIncidentFaces_sq, pathBoundaryPoints_sq, pathFacets_sq, pathFaces_sq, pathIndexRegions);
    printGeneralInfo(eik_g1);

    // we use the fmm olim method
    FMM_2D( eik_g1 );

    saveComputedValues(eik_g1, pathToSave);

    saveComputedGradients(eik_g1, pathSaveGradients);

    eik_grid_dealloc(&eik_g1);



    // // TEST SQUARE WITH INVERTED TRIANGLE - TWO SECTIONS

    // const char *pathPoints_sqTr, *pathNeighbors_sqTr, *pathIncidentFaces_sqTr, *pathBoundaryPoints_sqTr, *pathFacets_sqTr, *pathFaces_sqTr, *pathIndexRegionsTr, *pathToSaveTr;
    // pathPoints_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/MeshPoints.txt";
    // pathNeighbors_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/Neigh.txt";
    // pathIncidentFaces_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/IncidentFaces.txt";
    // pathBoundaryPoints_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/BoundaryPoints.txt";
    // pathFacets_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/Facets.txt";
    // pathFaces_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/Faces.txt";
    // pathIndexRegionsTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/FacesLabel.txt";

    // pathToSaveTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/ComputedValues.bin";

    // int *start2;
    // int nStart2, s2;

    // s2 = 5;
    // start2 = &s2;
    // nStart2 = 1;
    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN INVERTED TRIANGLE - TWO SECTIONS \n\n\n\n");
    // eik_gridS *eik_g2;
    // eik_grid_alloc(&eik_g2);
    // eik_grid_initFromFile(eik_g2, start2, nStart2, pathPoints_sqTr, pathNeighbors_sqTr, pathIncidentFaces_sqTr, pathBoundaryPoints_sqTr, pathFacets_sqTr, pathFaces_sqTr, pathIndexRegionsTr);
    // printGeneralInfo(eik_g2);

    // FMM_2D( eik_g2 );

    // saveComputedValues(eik_g2, pathToSaveTr);

    // eik_grid_dealloc(&eik_g2);

    // const char *pathToSaveTr2;

    // pathToSaveTr2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/ComputedValues2.bin";


    // s2 = 432;
    // start2 = &s2;
    // nStart2 = 1;
    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN INVERTED TRIANGLE - TWO SECTIONS OTHER STARTING POINT\n\n\n\n");
    // eik_gridS *eik_g3;
    // eik_grid_alloc(&eik_g3);
    // eik_grid_initFromFile(eik_g3, start2, nStart2, pathPoints_sqTr, pathNeighbors_sqTr, pathIncidentFaces_sqTr, pathBoundaryPoints_sqTr, pathFacets_sqTr, pathFaces_sqTr, pathIndexRegionsTr);
    // printGeneralInfo(eik_g3);

    // FMM_2D( eik_g3 );

    // saveComputedValues(eik_g3, pathToSaveTr2);

    // eik_grid_dealloc(&eik_g3);


    // // TEST SQUARE WITH ARCH - TWO SECTIONS

    // const char *pathPoints_sqArch, *pathNeighbors_sqArch, *pathIncidentFaces_sqArch, *pathBoundaryPoints_sqArch, *pathFacets_sqArch, *pathFaces_sqArch, *pathIndexRegionsArch, *pathToSaveArch;
    // pathPoints_sqArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/MeshPoints.txt";
    // pathNeighbors_sqArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/Neigh.txt";
    // pathIncidentFaces_sqArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/IncidentFaces.txt";
    // pathBoundaryPoints_sqArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/BoundaryPoints.txt";
    // pathFacets_sqArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/Facets.txt";
    // pathFaces_sqArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/Faces.txt";
    // pathIndexRegionsArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/FacesLabel.txt";

    // pathToSaveArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/ComputedValues.bin";

    // int *start3;
    // int nStart3, s3;

    // s3 = 217;
    // start3 = &s3;
    // nStart3 = 1;
    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN ARCH - TWO SECTIONS \n\n\n\n");
    // eik_gridS *eik_g4;
    // eik_grid_alloc(&eik_g4);
    // eik_grid_initFromFile(eik_g4, start3, nStart3, pathPoints_sqArch, pathNeighbors_sqArch, pathIncidentFaces_sqArch, pathBoundaryPoints_sqArch, pathFacets_sqArch, pathFaces_sqArch, pathIndexRegionsArch);
    // printGeneralInfo(eik_g4);

    // FMM_2D( eik_g4 );

    // saveComputedValues(eik_g4, pathToSaveArch);

    // eik_grid_dealloc(&eik_g4);

    // int s3_2[2];
    // s3_2[0] = 25;
    // s3_2[1] = 10;
    // start3 = &s3;
    // nStart3 = 2;

    // pathToSaveArch = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestArch/ComputedValues2.bin";

    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN ARCH - TWO SECTIONS AND TWO STARTING POINTS \n\n\n\n");
    // eik_gridS *eik_g5;
    // eik_grid_alloc(&eik_g5);
    // eik_grid_initFromFile(eik_g5, start3, nStart3, pathPoints_sqArch, pathNeighbors_sqArch, pathIncidentFaces_sqArch, pathBoundaryPoints_sqArch, pathFacets_sqArch, pathFaces_sqArch, pathIndexRegionsArch);
    // printGeneralInfo(eik_g5);

    // FMM_2D( eik_g5 );

    // saveComputedValues(eik_g5, pathToSaveArch);

    // eik_grid_dealloc(&eik_g5);


    





}