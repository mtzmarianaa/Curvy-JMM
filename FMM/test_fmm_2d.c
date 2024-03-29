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
    double times[7];
    clock_t start_t, end_t;
    double time;
    // // // TEST SQUARE WITH INVERTED TRIANGLE - TWO SECTIONS


    //   H1

    const char *pathPoints_H1, *pathNeighbors_H1, *pathIncidentFaces_H1, *pathBoundaryPoints_H1, *pathFacets_H1, *pathFaces_H1, *pathIndexRegions_H1, *pathToSaveTr_H1_, *pathSaveGradientsTr_H1_;
    const char *pathSavePath_H1, *pathSaveLambdas_H1;
    pathPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_MeshPoints.txt";
    pathNeighbors_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Neigh.txt";
    pathIncidentFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_IncidentFaces.txt";
    pathBoundaryPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_BoundaryPoints.txt";
    pathFacets_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Facets.txt";
    pathFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Faces.txt";
    pathIndexRegions_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_FacesLabel.txt";

    pathToSaveTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_Parents_ARTIFICIAL.bin";
    pathSaveLambdas_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H1/H1_LambdasOpt_ARTIFICIAL.bin";

    int *start_H1;
    int nstart_H1, s_H1;

    s_H1 = 0;
    start_H1 = &s_H1;
    nstart_H1 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H1 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start_H1, nstart_H1, pathPoints_H1, pathNeighbors_H1, pathIncidentFaces_H1, pathBoundaryPoints_H1, pathFacets_H1, pathFaces_H1, pathIndexRegions_H1);
    printGeneralInfo(eik_g1);

    start_t = clock();
    FMM_2D( eik_g1 , 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[0] = time;
    printGeneralInfo(eik_g1);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g1, pathToSaveTr_H1_);

    saveComputedGradients(eik_g1, pathSaveGradientsTr_H1_);

    saveComputedParents(eik_g1, pathSavePath_H1);

    saveComputedLambdas(eik_g1, pathSaveLambdas_H1);

    eik_grid_dealloc(&eik_g1);


    //   H2

    const char *pathPoints_H2, *pathNeighbors_H2, *pathIncidentFaces_H2, *pathBoundaryPoints_H2, *pathFacets_H2, *pathFaces_H2, *pathIndexRegions_H2, *pathToSaveTr_H2_, *pathSaveGradientsTr_H2_;
    const char *pathSavePath_H2, *pathSaveLambdas_H2;
    pathPoints_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_MeshPoints.txt";
    pathNeighbors_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Neigh.txt";
    pathIncidentFaces_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_IncidentFaces.txt";
    pathBoundaryPoints_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_BoundaryPoints.txt";
    pathFacets_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Facets.txt";
    pathFaces_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Faces.txt";
    pathIndexRegions_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_FacesLabel.txt";

    pathToSaveTr_H2_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_H2_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_Parents_ARTIFICIAL.bin";
    pathSaveLambdas_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H2/H2_LambdasOpt_ARTIFICIAL.bin";

    int *start_H2;
    int nstart_H2, s_H2;

    s_H2 = 0;
    start_H2 = &s_H2;
    nstart_H2 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H2 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g2;
    eik_grid_alloc(&eik_g2);
    eik_grid_initFromFile(eik_g2, start_H2, nstart_H2, pathPoints_H2, pathNeighbors_H2, pathIncidentFaces_H2, pathBoundaryPoints_H2, pathFacets_H2, pathFaces_H2, pathIndexRegions_H2);
    printGeneralInfo(eik_g2);

    start_t = clock();
    FMM_2D( eik_g2 , 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[1] = time;
    printGeneralInfo(eik_g2);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g2, pathToSaveTr_H2_);

    saveComputedGradients(eik_g2, pathSaveGradientsTr_H2_);

    saveComputedParents(eik_g2, pathSavePath_H2);

    saveComputedLambdas(eik_g2, pathSaveLambdas_H2);

    eik_grid_dealloc(&eik_g2);



    //   H3

    const char *pathPoints_H3, *pathNeighbors_H3, *pathIncidentFaces_H3, *pathBoundaryPoints_H3, *pathFacets_H3, *pathFaces_H3, *pathIndexRegions_H3, *pathToSaveTr_H3_, *pathSaveGradientsTr_H3_;
    const char *pathSavePath_H3, *pathSaveLambdas_H3;
    pathPoints_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_MeshPoints.txt";
    pathNeighbors_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_Neigh.txt";
    pathIncidentFaces_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_IncidentFaces.txt";
    pathBoundaryPoints_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_BoundaryPoints.txt";
    pathFacets_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_Facets.txt";
    pathFaces_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_Faces.txt";
    pathIndexRegions_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_FacesLabel.txt";

    pathToSaveTr_H3_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_H3_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_Parents_ARTIFICIAL.bin";
    pathSaveLambdas_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H3/H3_LambdasOpt_ARTIFICIAL.bin";

    int *start_H3;
    int nstart_H3, s_H3;

    s_H3 = 0;
    start_H3 = &s_H3;
    nstart_H3 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H3 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g3;
    eik_grid_alloc(&eik_g3);
    eik_grid_initFromFile(eik_g3, start_H3, nstart_H3, pathPoints_H3, pathNeighbors_H3, pathIncidentFaces_H3, pathBoundaryPoints_H3, pathFacets_H3, pathFaces_H3, pathIndexRegions_H3);
    printGeneralInfo(eik_g3);

    start_t = clock();
    FMM_2D( eik_g3 , 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[2] = time;
    printGeneralInfo(eik_g3);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g3, pathToSaveTr_H3_);

    saveComputedGradients(eik_g3, pathSaveGradientsTr_H3_);

    saveComputedParents(eik_g3, pathSavePath_H3);

    saveComputedLambdas(eik_g3, pathSaveLambdas_H3);

    eik_grid_dealloc(&eik_g3);


    //   H4

    const char *pathPoints_H4, *pathNeighbors_H4, *pathIncidentFaces_H4, *pathBoundaryPoints_H4, *pathFacets_H4, *pathFaces_H4, *pathIndexRegions_H4, *pathToSaveTr_H4_, *pathSaveGradientsTr_H4_;
    const char *pathSavePath_H4, *pathSaveLambdas_H4;
    pathPoints_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_MeshPoints.txt";
    pathNeighbors_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_Neigh.txt";
    pathIncidentFaces_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_IncidentFaces.txt";
    pathBoundaryPoints_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_BoundaryPoints.txt";
    pathFacets_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_Facets.txt";
    pathFaces_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_Faces.txt";
    pathIndexRegions_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_FacesLabel.txt";

    pathToSaveTr_H4_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_H4_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_Parents_ARTIFICIAL.bin";
    pathSaveLambdas_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H4/H4_LambdasOpt_ARTIFICIAL.bin";

    int *start_H4;
    int nstart_H4, s_H4;

    s_H4 = 0;
    start_H4 = &s_H4;
    nstart_H4 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H4 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g4;
    eik_grid_alloc(&eik_g4);
    eik_grid_initFromFile(eik_g4, start_H4, nstart_H4, pathPoints_H4, pathNeighbors_H4, pathIncidentFaces_H4, pathBoundaryPoints_H4, pathFacets_H4, pathFaces_H4, pathIndexRegions_H4);
    printGeneralInfo(eik_g4);

    start_t = clock();
    FMM_2D( eik_g4 , 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[3] = time;
    printGeneralInfo(eik_g4);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g4, pathToSaveTr_H4_);

    saveComputedGradients(eik_g4, pathSaveGradientsTr_H4_);

    saveComputedParents(eik_g4, pathSavePath_H4);

    saveComputedLambdas(eik_g4, pathSaveLambdas_H4);

    eik_grid_dealloc(&eik_g4);


    //   H5

    const char *pathPoints_H5, *pathNeighbors_H5, *pathIncidentFaces_H5, *pathBoundaryPoints_H5, *pathFacets_H5, *pathFaces_H5, *pathIndexRegions_H5, *pathToSaveTr_H5_, *pathSaveGradientsTr_H5_;
    const char *pathSavePath_H5, *pathSaveLambdas_H5;
    pathPoints_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_MeshPoints.txt";
    pathNeighbors_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_Neigh.txt";
    pathIncidentFaces_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_IncidentFaces.txt";
    pathBoundaryPoints_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_BoundaryPoints.txt";
    pathFacets_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_Facets.txt";
    pathFaces_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_Faces.txt";
    pathIndexRegions_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_FacesLabel.txt";

    pathToSaveTr_H5_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_H5_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_Parents_ARTIFICIAL.bin";
    pathSaveLambdas_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H5/H5_LambdasOpt_ARTIFICIAL.bin";

    int *start_H5;
    int nstart_H5, s_H5;

    s_H5 = 0;
    start_H5 = &s_H5;
    nstart_H5 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H5 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g5;
    eik_grid_alloc(&eik_g5);
    eik_grid_initFromFile(eik_g5, start_H5, nstart_H5, pathPoints_H5, pathNeighbors_H5, pathIncidentFaces_H5, pathBoundaryPoints_H5, pathFacets_H5, pathFaces_H5, pathIndexRegions_H5);
    printGeneralInfo(eik_g5);

    start_t = clock();
    FMM_2D( eik_g5, 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[4] = time;
    printGeneralInfo(eik_g5);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g5, pathToSaveTr_H5_);

    saveComputedGradients(eik_g5, pathSaveGradientsTr_H5_);

    saveComputedParents(eik_g5, pathSavePath_H5);

    saveComputedLambdas(eik_g5, pathSaveLambdas_H5);

    eik_grid_dealloc(&eik_g5);


    //   H6

    const char *pathPoints_H6, *pathNeighbors_H6, *pathIncidentFaces_H6, *pathBoundaryPoints_H6, *pathFacets_H6, *pathFaces_H6, *pathIndexRegions_H6, *pathToSaveTr_H6_, *pathSaveGradientsTr_H6_;
    const char *pathSavePath_H6, *pathSaveLambdas_H6;
    pathPoints_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_MeshPoints.txt";
    pathNeighbors_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_Neigh.txt";
    pathIncidentFaces_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_IncidentFaces.txt";
    pathBoundaryPoints_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_BoundaryPoints.txt";
    pathFacets_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_Facets.txt";
    pathFaces_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_Faces.txt";
    pathIndexRegions_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_FacesLabel.txt";

    pathToSaveTr_H6_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_H6_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_Parents_ARTIFICIAL.bin";
    pathSaveLambdas_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H6/H6_LambdasOpt_ARTIFICIAL.bin";

    int *start_H6;
    int nstart_H6, s_H6;

    s_H6 = 0;
    start_H6 = &s_H6;
    nstart_H6 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H6 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g6;
    eik_grid_alloc(&eik_g6);
    eik_grid_initFromFile(eik_g6, start_H6, nstart_H6, pathPoints_H6, pathNeighbors_H6, pathIncidentFaces_H6, pathBoundaryPoints_H6, pathFacets_H6, pathFaces_H6, pathIndexRegions_H6);
    printGeneralInfo(eik_g6);

    start_t = clock();
    FMM_2D( eik_g6, 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[5] = time;
    printGeneralInfo(eik_g6);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g6, pathToSaveTr_H6_);

    saveComputedGradients(eik_g6, pathSaveGradientsTr_H6_);

    saveComputedParents(eik_g6, pathSavePath_H6);

    saveComputedLambdas(eik_g6, pathSaveLambdas_H6);

    eik_grid_dealloc(&eik_g6);


    //   H0

    const char *pathPoints_H0, *pathNeighbors_H0, *pathIncidentFaces_H0, *pathBoundaryPoints_H0, *pathFacets_H0, *pathFaces_H0, *pathIndexRegions_H0, *pathToSaveTr_H0_, *pathSaveGradientsTr_H0_;
    const char *pathSavePath_H0, *pathSaveLambdas_H0;
    pathPoints_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_MeshPoints.txt";
    pathNeighbors_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_Neigh.txt";
    pathIncidentFaces_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_IncidentFaces.txt";
    pathBoundaryPoints_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_BoundaryPoints.txt";
    pathFacets_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_Facets.txt";
    pathFaces_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_Faces.txt";
    pathIndexRegions_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_FacesLabel.txt";

    pathToSaveTr_H0_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_H0_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_Parents_ARTIFICIAL.bin";
    pathSaveLambdas_H0 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/H-2/H-2_LambdasOpt_ARTIFICIAL.bin";

    int *start_H0;
    int nstart_H0, s_H0;

    s_H0 = 0;
    start_H0 = &s_H0;
    nstart_H0 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H0 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g7;
    eik_grid_alloc(&eik_g7);
    eik_grid_initFromFile(eik_g7, start_H0, nstart_H0, pathPoints_H0, pathNeighbors_H0, pathIncidentFaces_H0, pathBoundaryPoints_H0, pathFacets_H0, pathFaces_H0, pathIndexRegions_H0);
    printGeneralInfo(eik_g7);

    start_t = clock();
    FMM_2D( eik_g7, 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[6] = time;
    printGeneralInfo(eik_g7);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g7, pathToSaveTr_H0_);

    saveComputedGradients(eik_g7, pathSaveGradientsTr_H0_);

    saveComputedParents(eik_g7, pathSavePath_H0);

    saveComputedLambdas(eik_g7, pathSaveLambdas_H0);

    eik_grid_dealloc(&eik_g7);


    const char *pathTimes;

    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestBaseSnow/Times_ARTIFICIAL.bin";

    saveTimes(times, pathTimes);


}