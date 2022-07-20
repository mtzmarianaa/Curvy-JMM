/*
TESTS FOR THE GRID
*/
#include "eik_grid.h"
#include "fmm_2d.h"

#include <math.h>
#include <stdio.h>
#include <time.h>


int main()
{
    double times[5];
    clock_t start_t, end_t;
    double time;
    // // // TEST SQUARE WITH INVERTED TRIANGLE - TWO SECTIONS


    //   H1

    const char *pathPoints_H1_sqTr, *pathNeighbors_H1_sqTr, *pathIncidentFaces_H1_sqTr, *pathBoundaryPoints_H1_sqTr, *pathFacets_H1_sqTr, *pathFaces_H1_sqTr, *pathIndexRegions_H1_sqTr, *pathToSaveTr_H1_, *pathSaveGradientsTr_H1_;
    pathPoints_H1_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_MeshPoints.txt";
    pathNeighbors_H1_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_Neigh.txt";
    pathIncidentFaces_H1_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_IncidentFaces.txt";
    pathBoundaryPoints_H1_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_BoundaryPoints.txt";
    pathFacets_H1_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_Facets.txt";
    pathFaces_H1_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_Faces.txt";
    pathIndexRegions_H1_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_FacesLabel.txt";

    pathToSaveTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_ComputedValues.bin";
    pathSaveGradientsTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_ComputedGradients.bin";

    int *start_H1_2;
    int nstart_H1_2, s_H1;

    s_H1 = 5;
    start_H1_2 = &s_H1;
    nstart_H1_2 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN INVERTED TRIANGLE - TWO SECTIONS \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start_H1_2, nstart_H1_2, pathPoints_H1_sqTr, pathNeighbors_H1_sqTr, pathIncidentFaces_H1_sqTr, pathBoundaryPoints_H1_sqTr, pathFacets_H1_sqTr, pathFaces_H1_sqTr, pathIndexRegions_H1_sqTr);
    printGeneralInfo(eik_g1);

    start_t = clock();
    FMM_2D( eik_g1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[0] = time;
    printGeneralInfo(eik_g1);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g1, pathToSaveTr_H1_);

    saveComputedGradients(eik_g1, pathSaveGradientsTr_H1_);

    eik_grid_dealloc(&eik_g1);


    //   H2

    const char *pathPoints_H2_sqTr, *pathNeighbors_H2_sqTr, *pathIncidentFaces_H2_sqTr, *pathBoundaryPoints_H2_sqTr, *pathFacets_H2_sqTr, *pathFaces_H2_sqTr, *pathIndexRegions_H2_sqTr, *pathToSaveTr_H2_, *pathSaveGradientsTr_H2_;
    pathPoints_H2_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_MeshPoints.txt";
    pathNeighbors_H2_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_Neigh.txt";
    pathIncidentFaces_H2_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_IncidentFaces.txt";
    pathBoundaryPoints_H2_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_BoundaryPoints.txt";
    pathFacets_H2_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_Facets.txt";
    pathFaces_H2_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_Faces.txt";
    pathIndexRegions_H2_sqTr= "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_FacesLabel.txt";

    pathToSaveTr_H2_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_ComputedValues.bin";
    pathSaveGradientsTr_H2_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H2_ComputedGradients.bin";

    int *start_H2_2;
    int nstart_H2_2, s_H2;

    s_H2 = 5;
    start_H2_2 = &s_H2;
    nstart_H2_2 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN INVERTED TRIANGLE - TWO SECTIONS \n\n\n\n");
    eik_gridS *eik_g2;
    eik_grid_alloc(&eik_g2);
    eik_grid_initFromFile(eik_g2, start_H2_2, nstart_H2_2, pathPoints_H2_sqTr, pathNeighbors_H2_sqTr, pathIncidentFaces_H2_sqTr, pathBoundaryPoints_H2_sqTr, pathFacets_H2_sqTr, pathFaces_H2_sqTr, pathIndexRegions_H2_sqTr);
    printGeneralInfo(eik_g2);

    start_t = clock();
    FMM_2D( eik_g2 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[1] = time;
    printGeneralInfo(eik_g2);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g2, pathToSaveTr_H2_);

    saveComputedGradients(eik_g2, pathSaveGradientsTr_H2_);

    eik_grid_dealloc(&eik_g2);


    //   H3

    const char *pathPoints_H3_sqTr, *pathNeighbors_H3_sqTr, *pathIncidentFaces_H3_sqTr, *pathBoundaryPoints_H3_sqTr, *pathFacets_H3_sqTr, *pathFaces_H3_sqTr, *pathIndexRegions_H3_sqTr, *pathToSaveTr_H3_, *pathSaveGradientsTr_H3_;
    pathPoints_H3_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_MeshPoints.txt";
    pathNeighbors_H3_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_Neigh.txt";
    pathIncidentFaces_H3_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_IncidentFaces.txt";
    pathBoundaryPoints_H3_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_BoundaryPoints.txt";
    pathFacets_H3_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_Facets.txt";
    pathFaces_H3_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_Faces.txt";
    pathIndexRegions_H3_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_FacesLabel.txt";

    pathToSaveTr_H3_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_ComputedValues.bin";
    pathSaveGradientsTr_H3_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H3_ComputedGradients.bin";

    int *start_H3_2;
    int nstart_H3_2, s_H3;

    s_H3 = 5;
    start_H3_2 = &s_H3;
    nstart_H3_2 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN INVERTED TRIANGLE - TWO SECTIONS \n\n\n\n");
    eik_gridS *eik_g3;
    eik_grid_alloc(&eik_g3);
    eik_grid_initFromFile(eik_g3, start_H3_2, nstart_H3_2, pathPoints_H3_sqTr, pathNeighbors_H3_sqTr, pathIncidentFaces_H3_sqTr, pathBoundaryPoints_H3_sqTr, pathFacets_H3_sqTr, pathFaces_H3_sqTr, pathIndexRegions_H3_sqTr);
    printGeneralInfo(eik_g3);

    start_t = clock();
    FMM_2D( eik_g3 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[2] = time;
    printGeneralInfo(eik_g3);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g3, pathToSaveTr_H3_);

    saveComputedGradients(eik_g3, pathSaveGradientsTr_H3_);

    eik_grid_dealloc(&eik_g3);



    //   H4

    const char *pathPoints_H4_sqTr, *pathNeighbors_H4_sqTr, *pathIncidentFaces_H4_sqTr, *pathBoundaryPoints_H4_sqTr, *pathFacets_H4_sqTr, *pathFaces_H4_sqTr, *pathIndexRegions_H4_sqTr, *pathToSaveTr_H4_, *pathSaveGradientsTr_H4_;
    pathPoints_H4_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_MeshPoints.txt";
    pathNeighbors_H4_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_Neigh.txt";
    pathIncidentFaces_H4_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_IncidentFaces.txt";
    pathBoundaryPoints_H4_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_BoundaryPoints.txt";
    pathFacets_H4_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_Facets.txt";
    pathFaces_H4_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_Faces.txt";
    pathIndexRegions_H4_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_FacesLabel.txt";

    pathToSaveTr_H4_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_ComputedValues.bin";
    pathSaveGradientsTr_H4_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H4_ComputedGradients.bin";

    int *start_H4_2;
    int nstart_H4_2, s_H4;

    s_H4 = 5;
    start_H4_2 = &s_H4;
    nstart_H4_2 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN INVERTED TRIANGLE - TWO SECTIONS \n\n\n\n");
    eik_gridS *eik_g4;
    eik_grid_alloc(&eik_g4);
    eik_grid_initFromFile(eik_g4, start_H4_2, nstart_H4_2, pathPoints_H4_sqTr, pathNeighbors_H4_sqTr, pathIncidentFaces_H4_sqTr, pathBoundaryPoints_H4_sqTr, pathFacets_H4_sqTr, pathFaces_H4_sqTr, pathIndexRegions_H4_sqTr);
    printGeneralInfo(eik_g4);

    start_t = clock();
    FMM_2D( eik_g4 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[3] = time;
    printGeneralInfo(eik_g4);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g4, pathToSaveTr_H4_);

    saveComputedGradients(eik_g4, pathSaveGradientsTr_H4_);

    eik_grid_dealloc(&eik_g4);



    //   H6

    const char *pathPoints_H6_sqTr, *pathNeighbors_H6_sqTr, *pathIncidentFaces_H6_sqTr, *pathBoundaryPoints_H6_sqTr, *pathFacets_H6_sqTr, *pathFaces_H6_sqTr, *pathIndexRegions_H6_sqTr, *pathToSaveTr_H6_, *pathSaveGradientsTr_H6_;
    pathPoints_H6_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_MeshPoints.txt";
    pathNeighbors_H6_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_Neigh.txt";
    pathIncidentFaces_H6_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_IncidentFaces.txt";
    pathBoundaryPoints_H6_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_BoundaryPoints.txt";
    pathFacets_H6_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_Facets.txt";
    pathFaces_H6_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_Faces.txt";
    pathIndexRegions_H6_sqTr = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_FacesLabel.txt";

    pathToSaveTr_H6_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_ComputedValues.bin";
    pathSaveGradientsTr_H6_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H6_ComputedGradients.bin";

    int *start_H6_2;
    int nstart_H6_2, s_H6;

    s_H6 = 5;
    start_H6_2 = &s_H6;
    nstart_H6_2 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR A SQUARE WITH AN INVERTED TRIANGLE - TWO SECTIONS \n\n\n\n");
    eik_gridS *eik_g6;
    eik_grid_alloc(&eik_g6);
    eik_grid_initFromFile(eik_g6, start_H6_2, nstart_H6_2, pathPoints_H6_sqTr, pathNeighbors_H6_sqTr, pathIncidentFaces_H6_sqTr, pathBoundaryPoints_H6_sqTr, pathFacets_H6_sqTr, pathFaces_H6_sqTr, pathIndexRegions_H6_sqTr);
    printGeneralInfo(eik_g6);

    start_t = clock();
    FMM_2D( eik_g6 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[4] = time;
    printGeneralInfo(eik_g6);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g6, pathToSaveTr_H6_);

    saveComputedGradients(eik_g6, pathSaveGradientsTr_H6_);

    eik_grid_dealloc(&eik_g6);

    const char *pathTimes;

    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/Times.bin";

    saveTimes(times, pathTimes);


}