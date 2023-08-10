#include "mesh2D.h"
#include "eik_grid.h"
#include "marcher_T2.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


int main(){

    const char *pathPoints, *pathFaces, *pathEdges, *pathEdgesInFace, *pathNeighbors;
    const char *pathIncidentFaces, *pathIndices, *pathBoundary;
    const char *pathSaveEiks, *pathSaveGrads, *pathTimes;
    double rBall = 0.5;
    double time_taken[1], time;
    clock_t start_t, end_t;

    size_t *start, start_int, nStart;

    start_int = 0;
    start = &start_int;
    nStart = 1;

    /* pathPoints = "./EasyGeom/H00/H00_MeshPoints.txt"; */
    /* pathFaces = "./EasyGeom/H00/H00_Faces.txt"; */
    /* pathEdges = "./EasyGeom/H00/H00_Edges.txt"; */
    /* pathEdgesInFace = "./EasyGeom/H00/H00_EdgesInFace.txt"; */
    /* pathNeighbors = "./EasyGeom/H00/H00_Neigh.txt"; */
    /* pathIncidentFaces = "./EasyGeom/H00/H00_IncidentFaces.txt"; */
    /* pathIndices = "./EasyGeom/H00/H00_Indices.txt"; */
    /* pathIndices = "./EasyGeom/H00/H00_Indices.txt"; */
    /* pathBoundary = "./EasyGeom/H00/H00_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./EasyGeom/H00/H00_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./EasyGeom/H00/H00_ComputedGradientsFast.bin"; */
    /* pathTimes = "./EasyGeom/H00/H00_TimesFast.bin"; */

    /* // now we test the init with just the path to the files */

    /* printf("\n------------------------------------"); */
    /* printf("\n------------------------------------"); */
    /* printf("\n------------------------------------"); */
    /* printf("\n\n\n TESTING TH0E UPDATES WITH0 CUBIC H0ERMITE INTERPOLATION \n\n\n\n"); */

    /* eik_gridS *eik_g1; */
    /* eik_grid_alloc(&eik_g1); */

    /* eik_grid_initFromFile(eik_g1, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace, */
    /* 			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary); */

    /* printEverythingInMesh(eik_g1->mesh2); */

    /* start_t = clock(); */
    /* marcher_T2(eik_g1, rBall);; */
    /* end_t = clock(); */
    /* time = (double)(end_t - start_t)/ CLOCKS_PER_SEC; */
    /* time_taken[0] = time; */
    

    /* saveComputedValues(eik_g1, pathSaveEiks); */
    /* saveComputedGradients(eik_g1, pathSaveGrads); */
    /* saveTimes(time_taken, pathTimes); */

    /* printGeneralInfo(eik_g1); */

    /* eik_grid_dealloc(&eik_g1); */



    pathPoints = "./EasyGeom/H10/H10_MeshPoints.txt";
    pathFaces = "./EasyGeom/H10/H10_Faces.txt";
    pathEdges = "./EasyGeom/H10/H10_Edges.txt";
    pathEdgesInFace = "./EasyGeom/H10/H10_EdgesInFace.txt";
    pathNeighbors = "./EasyGeom/H10/H10_Neigh.txt";
    pathIncidentFaces = "./EasyGeom/H10/H10_IncidentFaces.txt";
    pathIndices = "./EasyGeom/H10/H10_Indices.txt";
    pathIndices = "./EasyGeom/H10/H10_Indices.txt";
    pathBoundary = "./EasyGeom/H10/H10_BoundaryCurve.txt";

    pathSaveEiks = "./EasyGeom/H10/H10_ComputedValuesFast.bin";
    pathSaveGrads = "./EasyGeom/H10/H10_ComputedGradientsFast.bin";
    pathTimes = "./EasyGeom/H10/H10_TimesFast.bin";

    eik_gridS *eik_g4;
    eik_grid_alloc(&eik_g4);

    eik_grid_initFromFile(eik_g4, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    printEverythingInMesh(eik_g4->mesh2);

    start_t = clock();
    marcher_T2(eik_g4, rBall);;
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g4, pathSaveEiks);
    saveComputedGradients(eik_g4, pathSaveGrads);
    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g4);

    eik_grid_dealloc(&eik_g4);




    /* pathPoints = "./EasyGeom/H01/H01_MeshPoints.txt"; */
    /* pathFaces = "./EasyGeom/H01/H01_Faces.txt"; */
    /* pathEdges = "./EasyGeom/H01/H01_Edges.txt"; */
    /* pathEdgesInFace = "./EasyGeom/H01/H01_EdgesInFace.txt"; */
    /* pathNeighbors = "./EasyGeom/H01/H01_Neigh.txt"; */
    /* pathIncidentFaces = "./EasyGeom/H01/H01_IncidentFaces.txt"; */
    /* pathIndices = "./EasyGeom/H01/H01_Indices.txt"; */
    /* pathIndices = "./EasyGeom/H01/H01_Indices.txt"; */
    /* pathBoundary = "./EasyGeom/H01/H01_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./EasyGeom/H01/H01_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./EasyGeom/H01/H01_ComputedGradientsFast.bin"; */
    /* pathTimes = "./EasyGeom/H01/H01_TimesFast.bin"; */

    /* eik_gridS *eik_g2; */
    /* eik_grid_alloc(&eik_g2); */

    /* eik_grid_initFromFile(eik_g2, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace, */
    /* 			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary); */

    /* printEverythingInMesh(eik_g2->mesh2); */

    /* start_t = clock(); */
    /* marcher_T2(eik_g2, rBall);; */
    /* end_t = clock(); */
    /* time = (double)(end_t - start_t)/CLOCKS_PER_SEC; */
    /* time_taken[0] = time; */
    

    /* saveComputedValues(eik_g2, pathSaveEiks); */
    /* saveComputedGradients(eik_g2, pathSaveGrads); */
    /* saveTimes(time_taken, pathTimes); */

    /* printGeneralInfo(eik_g2); */

    /* eik_grid_dealloc(&eik_g2); */



    /* pathPoints = "./EasyGeom/H11/H11_MeshPoints.txt"; */
    /* pathFaces = "./EasyGeom/H11/H11_Faces.txt"; */
    /* pathEdges = "./EasyGeom/H11/H11_Edges.txt"; */
    /* pathEdgesInFace = "./EasyGeom/H11/H11_EdgesInFace.txt"; */
    /* pathNeighbors = "./EasyGeom/H11/H11_Neigh.txt"; */
    /* pathIncidentFaces = "./EasyGeom/H11/H11_IncidentFaces.txt"; */
    /* pathIndices = "./EasyGeom/H11/H11_Indices.txt"; */
    /* pathIndices = "./EasyGeom/H11/H11_Indices.txt"; */
    /* pathBoundary = "./EasyGeom/H11/H11_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./EasyGeom/H11/H11_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./EasyGeom/H11/H11_ComputedGradientsFast.bin"; */
    /* pathTimes = "./EasyGeom/H11/H11_TimesFast.bin"; */

    /* eik_gridS *eik_g5; */
    /* eik_grid_alloc(&eik_g5); */

    /* eik_grid_initFromFile(eik_g5, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace, */
    /* 			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary); */

    /* printEverythingInMesh(eik_g5->mesh2); */


    /* start_t = clock(); */
    /* marcher_T2(eik_g5, rBall);; */
    /* end_t = clock(); */
    /* time = (double)(end_t - start_t)/ CLOCKS_PER_SEC; */
    /* time_taken[0] = time; */
    

    /* saveComputedValues(eik_g5, pathSaveEiks); */
    /* saveComputedGradients(eik_g5, pathSaveGrads); */
    /* saveTimes(time_taken, pathTimes); */

    /* printGeneralInfo(eik_g5); */

    /* eik_grid_dealloc(&eik_g5); */





    /* pathPoints = "./EasyGeom/H02/H02_MeshPoints.txt"; */
    /* pathFaces = "./EasyGeom/H02/H02_Faces.txt"; */
    /* pathEdges = "./EasyGeom/H02/H02_Edges.txt"; */
    /* pathEdgesInFace = "./EasyGeom/H02/H02_EdgesInFace.txt"; */
    /* pathNeighbors = "./EasyGeom/H02/H02_Neigh.txt"; */
    /* pathIncidentFaces = "./EasyGeom/H02/H02_IncidentFaces.txt"; */
    /* pathIndices = "./EasyGeom/H02/H02_Indices.txt"; */
    /* pathIndices = "./EasyGeom/H02/H02_Indices.txt"; */
    /* pathBoundary = "./EasyGeom/H02/H02_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./EasyGeom/H02/H02_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./EasyGeom/H02/H02_ComputedGradientsFast.bin"; */
    /* pathTimes = "./EasyGeom/H02/H02_TimesFast.bin"; */

    /* eik_gridS *eik_g3; */
    /* eik_grid_alloc(&eik_g3); */

    /* eik_grid_initFromFile(eik_g3, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace, */
    /* 			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary); */

    /* printEverythingInMesh(eik_g3->mesh2); */

    /* start_t = clock(); */
    /* marcher_T2(eik_g3, rBall);; */
    /* end_t = clock(); */
    /* time = (double)(end_t - start_t)/ CLOCKS_PER_SEC; */
    /* time_taken[0] = time; */
    

    /* saveComputedValues(eik_g3, pathSaveEiks); */
    /* saveComputedGradients(eik_g3, pathSaveGrads); */
    /* saveTimes(time_taken, pathTimes); */

    /* printGeneralInfo(eik_g3); */

    /* eik_grid_dealloc(&eik_g3); */


    




  

}
