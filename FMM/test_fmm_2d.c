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
    pathPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_MeshPoints.txt";
    pathNeighbors_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_Neigh.txt";
    pathIncidentFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_IncidentFaces.txt";
    pathBoundaryPoints_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_BoundaryPoints.txt";
    pathFacets_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_Facets.txt";
    pathFaces_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_Faces.txt";
    pathIndexRegions_H1 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_FacesLabel.txt";

    pathToSaveTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_Computed1-3Values.bin";
    pathSaveGradientsTr_H1_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H1/H1_Computed1-3Gradients.bin";

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

    const char *pathPoints_H2, *pathNeighbors_H2, *pathIncidentFaces_H2, *pathBoundaryPoints_H2, *pathFacets_H2, *pathFaces_H2, *pathIndexRegions_H2, *pathToSaveTr_H2_, *pathSaveGradientsTr_H2_;
    pathPoints_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_MeshPoints.txt";
    pathNeighbors_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_Neigh.txt";
    pathIncidentFaces_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_IncidentFaces.txt";
    pathBoundaryPoints_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_BoundaryPoints.txt";
    pathFacets_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_Facets.txt";
    pathFaces_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_Faces.txt";
    pathIndexRegions_H2 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_FacesLabel.txt";

    pathToSaveTr_H2_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_Computed1-3Values.bin";
    pathSaveGradientsTr_H2_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H2/H2_Computed1-3Gradients.bin";

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

    const char *pathPoints_H3, *pathNeighbors_H3, *pathIncidentFaces_H3, *pathBoundaryPoints_H3, *pathFacets_H3, *pathFaces_H3, *pathIndexRegions_H3, *pathToSaveTr_H3_, *pathSaveGradientsTr_H3_;
    pathPoints_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_MeshPoints.txt";
    pathNeighbors_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_Neigh.txt";
    pathIncidentFaces_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_IncidentFaces.txt";
    pathBoundaryPoints_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_BoundaryPoints.txt";
    pathFacets_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_Facets.txt";
    pathFaces_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_Faces.txt";
    pathIndexRegions_H3 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_FacesLabel.txt";

    pathToSaveTr_H3_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_Computed1-3Values.bin";
    pathSaveGradientsTr_H3_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H3/H3_Computed1-3Gradients.bin";

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

    const char *pathPoints_H4, *pathNeighbors_H4, *pathIncidentFaces_H4, *pathBoundaryPoints_H4, *pathFacets_H4, *pathFaces_H4, *pathIndexRegions_H4, *pathToSaveTr_H4_, *pathSaveGradientsTr_H4_;
    pathPoints_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_MeshPoints.txt";
    pathNeighbors_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_Neigh.txt";
    pathIncidentFaces_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_IncidentFaces.txt";
    pathBoundaryPoints_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_BoundaryPoints.txt";
    pathFacets_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_Facets.txt";
    pathFaces_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_Faces.txt";
    pathIndexRegions_H4 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_FacesLabel.txt";

    pathToSaveTr_H4_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_Computed1-3Values.bin";
    pathSaveGradientsTr_H4_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H4/H4_Computed1-3Gradients.bin";

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
    FMM_2D( eik_g4 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[3] = time;
    printGeneralInfo(eik_g4);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g4, pathToSaveTr_H4_);

    saveComputedGradients(eik_g4, pathSaveGradientsTr_H4_);

    eik_grid_dealloc(&eik_g4);


    //   H5

    const char *pathPoints_H5, *pathNeighbors_H5, *pathIncidentFaces_H5, *pathBoundaryPoints_H5, *pathFacets_H5, *pathFaces_H5, *pathIndexRegions_H5, *pathToSaveTr_H5_, *pathSaveGradientsTr_H5_;
    pathPoints_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_MeshPoints.txt";
    pathNeighbors_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_Neigh.txt";
    pathIncidentFaces_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_IncidentFaces.txt";
    pathBoundaryPoints_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_BoundaryPoints.txt";
    pathFacets_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_Facets.txt";
    pathFaces_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_Faces.txt";
    pathIndexRegions_H5 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_FacesLabel.txt";

    pathToSaveTr_H5_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_Computed1-3Values.bin";
    pathSaveGradientsTr_H5_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H5/H5_Computed1-3Gradients.bin";

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
    FMM_2D( eik_g5 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[4] = time;
    printGeneralInfo(eik_g5);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g5, pathToSaveTr_H5_);

    saveComputedGradients(eik_g5, pathSaveGradientsTr_H5_);

    eik_grid_dealloc(&eik_g5);


    //   H6

    const char *pathPoints_H6, *pathNeighbors_H6, *pathIncidentFaces_H6, *pathBoundaryPoints_H6, *pathFacets_H6, *pathFaces_H6, *pathIndexRegions_H6, *pathToSaveTr_H6_, *pathSaveGradientsTr_H6_;
    pathPoints_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_MeshPoints.txt";
    pathNeighbors_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_Neigh.txt";
    pathIncidentFaces_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_IncidentFaces.txt";
    pathBoundaryPoints_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_BoundaryPoints.txt";
    pathFacets_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_Facets.txt";
    pathFaces_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_Faces.txt";
    pathIndexRegions_H6 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_FacesLabel.txt";

    pathToSaveTr_H6_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_Computed1-3Values.bin";
    pathSaveGradientsTr_H6_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H6/H6_Computed1-3Gradients.bin";

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
    FMM_2D( eik_g6 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[5] = time;
    printGeneralInfo(eik_g6);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g6, pathToSaveTr_H6_);

    saveComputedGradients(eik_g6, pathSaveGradientsTr_H6_);

    eik_grid_dealloc(&eik_g6);


    //   H7

    const char *pathPoints_H7, *pathNeighbors_H7, *pathIncidentFaces_H7, *pathBoundaryPoints_H7, *pathFacets_H7, *pathFaces_H7, *pathIndexRegions_H7, *pathToSaveTr_H7_, *pathSaveGradientsTr_H7_;
    pathPoints_H7 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_MeshPoints.txt";
    pathNeighbors_H7 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_Neigh.txt";
    pathIncidentFaces_H7 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_IncidentFaces.txt";
    pathBoundaryPoints_H7 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_BoundaryPoints.txt";
    pathFacets_H7 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_Facets.txt";
    pathFaces_H7 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_Faces.txt";
    pathIndexRegions_H7 = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_FacesLabel.txt";

    pathToSaveTr_H7_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_Computed1-3Values.bin";
    pathSaveGradientsTr_H7_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/H7/H7_Computed1-3Gradients.bin";

    int *start_H7;
    int nstart_H7, s_H7;

    s_H7 = 0;
    start_H7 = &s_H7;
    nstart_H7 = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING FROM FILES OF A TRIANGULAR MESH FOR TEST GEOMERTY H7 - THREE SECTIONS \n\n\n\n");
    eik_gridS *eik_g7;
    eik_grid_alloc(&eik_g7);
    eik_grid_initFromFile(eik_g7, start_H7, nstart_H7, pathPoints_H7, pathNeighbors_H7, pathIncidentFaces_H7, pathBoundaryPoints_H7, pathFacets_H7, pathFaces_H7, pathIndexRegions_H7);
    printGeneralInfo(eik_g7);

    start_t = clock();
    FMM_2D( eik_g7 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    times[6] = time;
    printGeneralInfo(eik_g7);
    printf("Time taken %fl\n", time);
    

    saveComputedValues(eik_g7, pathToSaveTr_H7_);

    saveComputedGradients(eik_g7, pathSaveGradientsTr_H7_);

    eik_grid_dealloc(&eik_g7);


    const char *pathTimes;

    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/MeshInfo/Times.bin";

    saveTimes(times, pathTimes);


}