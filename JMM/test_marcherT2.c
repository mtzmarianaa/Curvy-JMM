// testing the marcher with hermite interpolation for T(xLambda)

#include "eik_grid.h"
#include "marcher_T2.h"

#include <math.h>
#include <stdio.h>
#include <time.h>


int main()
{
    // TEST GEOMETRY FOR DIFFERENT INDICES OF REFRACTION - SOME PATHS GO AROUND, NOT INTO REG3


    const char *pathPoints, *pathFaces, *pathEdges, *pathEdgesInFace, *pathNeighbors;
    const char *pathIncidentFaces, *pathIndices, *pathBoundary;
    const char *pathSaveEiks, *pathSaveGrads, *pathTimes;
    double rBall = 2.0;

    pathPoints = "./H0/H0_MeshPoints.txt";
    pathFaces = "./H0/H0_Faces.txt";
    pathEdges = "./H0/H0_Edges.txt";
    pathEdgesInFace = "./H0/H0_EdgesInFace.txt";
    pathNeighbors = "./H0/H0_Neigh.txt";
    pathIncidentFaces = "./H0/H0_IncidentFaces.txt";
    pathIndices = "./H0/H0_Indices.txt";
    pathIndices = "./H0/H0_Indices.txt";
    pathBoundary = "./H0/H0_BoundaryCurve.txt";

    pathSaveEiks = "./H0/H0_ComputedValues.bin";
    pathSaveGrads = "./H0/H0_ComputedGradients.bin";
    pathTimes = "./H0/H0_Times.bin";

    double time_taken[1], time;
    clock_t start_t, end_t;

    size_t *start, start_int, nStart;

    start_int = 0;
    start = &start_int;
    nStart = 1;
    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING THE UPDATES WITH CUBIC HERMITE INTERPOLATION \n\n\n\n");

    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);

    eik_grid_initFromFile(eik_g1, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    // printAllInfoMesh(eik_g1);

    start_t = clock();
    marcher_T2(eik_g1, rBall);;
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g1, pathSaveEiks);
    saveComputedGradients(eik_g1, pathSaveGrads);
    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g1);

    eik_grid_dealloc(&eik_g1);





    pathPoints = "./H1/H1_MeshPoints.txt";
    pathFaces = "./H1/H1_Faces.txt";
    pathEdges = "./H1/H1_Edges.txt";
    pathEdgesInFace = "./H1/H1_EdgesInFace.txt";
    pathNeighbors = "./H1/H1_Neigh.txt";
    pathIncidentFaces = "./H1/H1_IncidentFaces.txt";
    pathIndices = "./H1/H1_Indices.txt";
    pathIndices = "./H1/H1_Indices.txt";
    pathBoundary = "./H1/H1_BoundaryCurve.txt";

    pathSaveEiks = "./H1/H1_ComputedValues.bin";
    pathSaveGrads = "./H1/H1_ComputedGradients.bin";
    pathTimes = "./H1/H1_Times.bin";

    eik_gridS *eik_g2;
    eik_grid_alloc(&eik_g2);

    eik_grid_initFromFile(eik_g2, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    // printAllInfoMesh(eik_g1);

    start_t = clock();
    marcher_T2(eik_g2, rBall);;
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g2, pathSaveEiks);
    saveComputedGradients(eik_g2, pathSaveGrads);
    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g2);

    eik_grid_dealloc(&eik_g2);



}
