// testing the artificial triangle update

#include "eik_grid.h"
#include "fmm_2d.h"

#include <math.h>
#include <stdio.h>
#include <time.h>


int main()
{
    // TEST GEOMETRY FOR DIFFERENT INDICES OF REFRACTION - SOME PATHS GO AROUND, NOT INTO REG3



    const char *pathPoints, *pathNeighbors, *pathIncidentFaces, *pathBoundaryPoints, *pathFacets, *pathFaces, *pathIndexRegions, *pathToSaveTr_, *pathSaveGradientsTr_;
    const char *pathSavePath, *pathSaveLambdas, *pathTimes;
    int *start;
    int nstart, s;
    double time_taken[1], time;
    clock_t start_t, end_t;



    // pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_MeshPoints.txt";
    // pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_Neigh.txt";
    // pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_IncidentFaces.txt";
    // pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_BoundaryPoints.txt";
    // pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_Facets.txt";
    // pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_Faces.txt";
    // pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_FacesLabel.txt";

    // pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_ComputedValues.bin";
    // pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_ComputedGradients.bin";
    // pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_Parents.bin";
    // pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_LambdasOpt.bin";
    // pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_14/H1_14_Times.bin";


    // s = 0;
    // start = &s;
    // nstart = 1;
    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES H1_14 \n\n\n\n");
    // eik_gridS *eik_g1;
    // eik_grid_alloc(&eik_g1);
    // eik_grid_initFromFile(eik_g1, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    // // printAllInfoMesh(eik_g1);

    // start_t = clock();
    // FMM_2D( eik_g1, 0 );
    // end_t = clock();
    // time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    // time_taken[0] = time;
    

    // saveComputedValues(eik_g1, pathToSaveTr_);

    // saveComputedGradients(eik_g1, pathSaveGradientsTr_);

    // saveComputedParents(eik_g1, pathSavePath);

    // saveComputedLambdas(eik_g1, pathSaveLambdas);

    // saveTimes(time_taken, pathTimes);

    // printGeneralInfo(eik_g1);

    // eik_grid_dealloc(&eik_g1);


    // // SECOND

    // pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_MeshPoints.txt";
    // pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_Neigh.txt";
    // pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_IncidentFaces.txt";
    // pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_BoundaryPoints.txt";
    // pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_Facets.txt";
    // pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_Faces.txt";
    // pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_FacesLabel.txt";

    // pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_ComputedValues.bin";
    // pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_ComputedGradients.bin";
    // pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_Parents.bin";
    // pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_LambdasOpt.bin";
    // pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_13/H1_13_Times.bin";


    // s = 0;
    // start = &s;
    // nstart = 1;
    // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES H1_13 \n\n\n\n");
    // eik_gridS *eik_g2;
    // eik_grid_alloc(&eik_g2);
    // eik_grid_initFromFile(eik_g2, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    // printAllInfoMesh(eik_g2);

    // start_t = clock();
    // FMM_2D( eik_g2, 0 );
    // end_t = clock();
    // time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    // time_taken[0] = time;
    

    // saveComputedValues(eik_g2, pathToSaveTr_);

    // saveComputedGradients(eik_g2, pathSaveGradientsTr_);

    // saveComputedParents(eik_g2, pathSavePath);

    // saveComputedLambdas(eik_g2, pathSaveLambdas);

    // saveTimes(time_taken, pathTimes);

    // printGeneralInfo(eik_g2);

    // eik_grid_dealloc(&eik_g2);



    // THIRD

    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_Neigh.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_IncidentFaces.txt";
    pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_BoundaryPoints.txt";
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_Facets.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_FacesLabel.txt";

    pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_ComputedValues.bin";
    pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_ComputedGradients.bin";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_Parents.bin";
    pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_LambdasOpt.bin";
    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_1/H1_1_Times.bin";


    s = 0;
    start = &s;
    nstart = 1;
    //now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES H1_1 \n\n\n\n");
    eik_gridS *eik_g3;
    eik_grid_alloc(&eik_g3);
    eik_grid_initFromFile(eik_g3, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    printAllInfoMesh(eik_g3);

    start_t = clock();
    FMM_2D( eik_g3, 0 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g3, pathToSaveTr_);

    saveComputedGradients(eik_g3, pathSaveGradientsTr_);

    saveComputedParents(eik_g3, pathSavePath);

    saveComputedLambdas(eik_g3, pathSaveLambdas);

    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g3);

    eik_grid_dealloc(&eik_g3);



    // FOURTH

    // pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_MeshPoints.txt";
    // pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_Neigh.txt";
    // pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_IncidentFaces.txt";
    // pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_BoundaryPoints.txt";
    // pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_Facets.txt";
    // pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_Faces.txt";
    // pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_FacesLabel.txt";

    // pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_ComputedValues.bin";
    // pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_ComputedGradients.bin";
    // pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_Parents.bin";
    // pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_LambdasOpt.bin";
    // pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestTriangleSquare/H1_11/H1_11_Times.bin";


    // s = 0;
    // start = &s;
    // nstart = 1;
    // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES H1_11 \n\n\n\n");
    // eik_gridS *eik_g4;
    // eik_grid_alloc(&eik_g4);
    // eik_grid_initFromFile(eik_g4, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    // printAllInfoMesh(eik_g4);

    // start_t = clock();
    // FMM_2D( eik_g4, 0 );
    // end_t = clock();
    // time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    // time_taken[0] = time;
    

    // saveComputedValues(eik_g4, pathToSaveTr_);

    // saveComputedGradients(eik_g4, pathSaveGradientsTr_);

    // saveComputedParents(eik_g4, pathSavePath);

    // saveComputedLambdas(eik_g4, pathSaveLambdas);

    // saveTimes(time_taken, pathTimes);

    // printGeneralInfo(eik_g4);

    // eik_grid_dealloc(&eik_g4);



}