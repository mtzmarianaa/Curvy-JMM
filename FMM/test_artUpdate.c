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
    pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_MeshPoints.txt";
    pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_Neigh.txt";
    pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_IncidentFaces.txt";
    pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_BoundaryPoints.txt";
    pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_Facets.txt";
    pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_Faces.txt";
    pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_FacesLabel.txt";

    pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_ComputedValues_ARTIFICIAL.bin";
    pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_ComputedGradients_ARTIFICIAL.bin";
    pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_Parents_ARTIFICIAL.bin";
    pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_LambdasOpt_ARTIFICIAL.bin";
    pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H1/H1_Times_ARTIFICIAL.bin";

    int *start;
    int nstart, s;
    double time_taken[1], time;
    clock_t start_t, end_t;

    s = 0;
    start = &s;
    nstart = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES \n\n\n\n");
    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);
    eik_grid_initFromFile(eik_g1, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    // printAllInfoMesh(eik_g1);

    start_t = clock();
    FMM_2D( eik_g1, 1 );
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g1, pathToSaveTr_);

    saveComputedGradients(eik_g1, pathSaveGradientsTr_);

    saveComputedParents(eik_g1, pathSavePath);

    saveComputedLambdas(eik_g1, pathSaveLambdas);

    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g1);

    eik_grid_dealloc(&eik_g1);


    // // SECOND

    // pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_MeshPoints.txt";
    // pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_Neigh.txt";
    // pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_IncidentFaces.txt";
    // pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_BoundaryPoints.txt";
    // pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_Facets.txt";
    // pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_Faces.txt";
    // pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_FacesLabel.txt";

    // pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_ComputedValues_ARTIFICIAL.bin";
    // pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_ComputedGradients_ARTIFICIAL.bin";
    // pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_Parents_ARTIFICIAL.bin";
    // pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_LambdasOpt_ARTIFICIAL.bin";
    // pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H103/H103_Times_ARTIFICIAL.bin";


    // s = 0;
    // start = &s;
    // nstart = 1;
    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES \n\n\n\n");
    // eik_gridS *eik_g2;
    // eik_grid_alloc(&eik_g2);
    // eik_grid_initFromFile(eik_g2, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    // // printAllInfoMesh(eik_g2);

    // start_t = clock();
    // FMM_2D( eik_g2, 1 );
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



    // // THIRD

    // pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_MeshPoints.txt";
    // pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_Neigh.txt";
    // pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_IncidentFaces.txt";
    // pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_BoundaryPoints.txt";
    // pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_Facets.txt";
    // pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_Faces.txt";
    // pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_FacesLabel.txt";

    // pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_ComputedValues_ARTIFICIAL.bin";
    // pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_ComputedGradients_ARTIFICIAL.bin";
    // pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_Parents_ARTIFICIAL.bin";
    // pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_LambdasOpt_ARTIFICIAL.bin";
    // pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H102/H102_Times_ARTIFICIAL.bin";


    // s = 0;
    // start = &s;
    // nstart = 1;
    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES \n\n\n\n");
    // eik_gridS *eik_g3;
    // eik_grid_alloc(&eik_g3);
    // eik_grid_initFromFile(eik_g3, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    // // printAllInfoMesh(eik_g3);

    // start_t = clock();
    // FMM_2D( eik_g3, 1 );
    // end_t = clock();
    // time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    // time_taken[0] = time;
    

    // saveComputedValues(eik_g3, pathToSaveTr_);

    // saveComputedGradients(eik_g3, pathSaveGradientsTr_);

    // saveComputedParents(eik_g3, pathSavePath);

    // saveComputedLambdas(eik_g3, pathSaveLambdas);

    // saveTimes(time_taken, pathTimes);

    // printGeneralInfo(eik_g3);

    // eik_grid_dealloc(&eik_g3);



    // // FOURTH

    // pathPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_MeshPoints.txt";
    // pathNeighbors = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_Neigh.txt";
    // pathIncidentFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_IncidentFaces.txt";
    // pathBoundaryPoints = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_BoundaryPoints.txt";
    // pathFacets = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_Facets.txt";
    // pathFaces = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_Faces.txt";
    // pathIndexRegions = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_FacesLabel.txt";

    // pathToSaveTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_ComputedValues_ARTIFICIAL.bin";
    // pathSaveGradientsTr_ = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_ComputedGradients_ARTIFICIAL.bin";
    // pathSavePath = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_Parents_ARTIFICIAL.bin";
    // pathSaveLambdas = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_LambdasOpt_ARTIFICIAL.bin";
    // pathTimes = "/Users/marianamartinez/Documents/NYU-Courant/FMM-Project/FMM/TestSquareTriangle/H101/H101_Times_ARTIFICIAL.bin";


    // s = 0;
    // start = &s;
    // nstart = 1;
    // // now we test the init with just the path to the files

    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n------------------------------------");
    // printf("\n\n\n TESTING THE UPDATES WITH ARTIFICIAL TRIANGLES \n\n\n\n");
    // eik_gridS *eik_g4;
    // eik_grid_alloc(&eik_g4);
    // eik_grid_initFromFile(eik_g4, start, nstart, pathPoints, pathNeighbors, pathIncidentFaces, pathBoundaryPoints, pathFacets, pathFaces, pathIndexRegions);

    // // printAllInfoMesh(eik_g4);

    // start_t = clock();
    // FMM_2D( eik_g4, 1 );
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