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
    double rBall = 2.0;
    double time_taken[1], time;
    clock_t start_t, end_t;

    size_t *start, start_int, nStart;

    start_int = 0;
    start = &start_int;
    nStart = 1;

    pathPoints = "./H0/H0_MeshPoints.txt";
    pathFaces = "./H0/H0_Faces.txt";
    pathEdges = "./H0/H0_Edges.txt";
    pathEdgesInFace = "./H0/H0_EdgesInFace.txt";
    pathNeighbors = "./H0/H0_Neigh.txt";
    pathIncidentFaces = "./H0/H0_IncidentFaces.txt";
    pathIndices = "./H0/H0_Indices.txt";
    pathIndices = "./H0/H0_Indices.txt";
    pathBoundary = "./H0/H0_BoundaryCurve.txt";

    pathSaveEiks = "./H0/H0_ComputedValuesFast.bin";
    pathSaveGrads = "./H0/H0_ComputedGradientsFast.bin";
    pathTimes = "./H0/H0_TimesFast.bin";

    // now we test the init with just the path to the files

    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n------------------------------------");
    printf("\n\n\n TESTING THE UPDATES WITH CUBIC HERMITE INTERPOLATION \n\n\n\n");

    eik_gridS *eik_g1;
    eik_grid_alloc(&eik_g1);

    eik_grid_initFromFile(eik_g1, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    printEverythingInMesh(eik_g1->mesh2);

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

    pathSaveEiks = "./H1/H1_ComputedValuesFast.bin";
    pathSaveGrads = "./H1/H1_ComputedGradientsFast.bin";
    pathTimes = "./H1/H1_TimesFast.bin";

    eik_gridS *eik_g2;
    eik_grid_alloc(&eik_g2);

    eik_grid_initFromFile(eik_g2, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    printEverythingInMesh(eik_g2->mesh2);

    start_t = clock();
    marcher_T2(eik_g2, rBall);;
    end_t = clock();
    time = (double)(end_t - start_t)/CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g2, pathSaveEiks);
    saveComputedGradients(eik_g2, pathSaveGrads);
    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g2);

    eik_grid_dealloc(&eik_g2);


   





    pathPoints = "./H2/H2_MeshPoints.txt";
    pathFaces = "./H2/H2_Faces.txt";
    pathEdges = "./H2/H2_Edges.txt";
    pathEdgesInFace = "./H2/H2_EdgesInFace.txt";
    pathNeighbors = "./H2/H2_Neigh.txt";
    pathIncidentFaces = "./H2/H2_IncidentFaces.txt";
    pathIndices = "./H2/H2_Indices.txt";
    pathIndices = "./H2/H2_Indices.txt";
    pathBoundary = "./H2/H2_BoundaryCurve.txt";

    pathSaveEiks = "./H2/H2_ComputedValuesFast.bin";
    pathSaveGrads = "./H2/H2_ComputedGradientsFast.bin";
    pathTimes = "./H2/H2_TimesFast.bin";

    eik_gridS *eik_g3;
    eik_grid_alloc(&eik_g3);

    eik_grid_initFromFile(eik_g3, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    printEverythingInMesh(eik_g3->mesh2);

    start_t = clock();
    marcher_T2(eik_g3, rBall);;
    end_t = clock();
    time = (double)(end_t - start_t)/ CLOCKS_PER_SEC;
    time_taken[0] = time;
    

    saveComputedValues(eik_g3, pathSaveEiks);
    saveComputedGradients(eik_g3, pathSaveGrads);
    saveTimes(time_taken, pathTimes);

    printGeneralInfo(eik_g3);

    eik_grid_dealloc(&eik_g3);




    /* pathPoints = "./H3/H3_MeshPoints.txt"; */
    /* pathFaces = "./H3/H3_Faces.txt"; */
    /* pathEdges = "./H3/H3_Edges.txt"; */
    /* pathEdgesInFace = "./H3/H3_EdgesInFace.txt"; */
    /* pathNeighbors = "./H3/H3_Neigh.txt"; */
    /* pathIncidentFaces = "./H3/H3_IncidentFaces.txt"; */
    /* pathIndices = "./H3/H3_Indices.txt"; */
    /* pathIndices = "./H3/H3_Indices.txt"; */
    /* pathBoundary = "./H3/H3_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./H3/H3_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./H3/H3_ComputedGradientsFast.bin"; */
    /* pathTimes = "./H3/H3_TimesFast.bin"; */

    /* eik_gridS *eik_g4; */
    /* eik_grid_alloc(&eik_g4); */

    /* eik_grid_initFromFile(eik_g4, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace, */
    /* 			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary); */

    /* printEverythingInMesh(eik_g4->mesh2); */

    /* start_t = clock(); */
    /* marcher_T2(eik_g4, rBall);; */
    /* end_t = clock(); */
    /* time = (double)(end_t - start_t)/ CLOCKS_PER_SEC; */
    /* time_taken[0] = time; */
    

    /* saveComputedValues(eik_g4, pathSaveEiks); */
    /* saveComputedGradients(eik_g4, pathSaveGrads); */
    /* saveTimes(time_taken, pathTimes); */

    /* printGeneralInfo(eik_g4); */

    /* eik_grid_dealloc(&eik_g4); */


    

    /* pathPoints = "./H4/H4_MeshPoints.txt"; */
    /* pathFaces = "./H4/H4_Faces.txt"; */
    /* pathEdges = "./H4/H4_Edges.txt"; */
    /* pathEdgesInFace = "./H4/H4_EdgesInFace.txt"; */
    /* pathNeighbors = "./H4/H4_Neigh.txt"; */
    /* pathIncidentFaces = "./H4/H4_IncidentFaces.txt"; */
    /* pathIndices = "./H4/H4_Indices.txt"; */
    /* pathIndices = "./H4/H4_Indices.txt"; */
    /* pathBoundary = "./H4/H4_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./H4/H4_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./H4/H4_ComputedGradientsFast.bin"; */
    /* pathTimes = "./H4/H4_TimesFast.bin"; */

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

    


    /* pathPoints = "./H5/H5_MeshPoints.txt"; */
    /* pathFaces = "./H5/H5_Faces.txt"; */
    /* pathEdges = "./H5/H5_Edges.txt"; */
    /* pathEdgesInFace = "./H5/H5_EdgesInFace.txt"; */
    /* pathNeighbors = "./H5/H5_Neigh.txt"; */
    /* pathIncidentFaces = "./H5/H5_IncidentFaces.txt"; */
    /* pathIndices = "./H5/H5_Indices.txt"; */
    /* pathIndices = "./H5/H5_Indices.txt"; */
    /* pathBoundary = "./H5/H5_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./H5/H5_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./H5/H5_ComputedGradientsFast.bin"; */
    /* pathTimes = "./H5/H5_TimesFast.bin"; */

    /* eik_gridS *eik_g6; */
    /* eik_grid_alloc(&eik_g6); */

    /* eik_grid_initFromFile(eik_g6, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace, */
    /* 			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary); */

    /* printEverythingInMesh(eik_g6->mesh2); */

    /* start_t = clock(); */
    /* marcher_T2(eik_g6, rBall);; */
    /* end_t = clock(); */
    /* time = (double)(end_t - start_t)/ CLOCKS_PER_SEC; */
    /* time_taken[0] = time; */
    

    /* saveComputedValues(eik_g6, pathSaveEiks); */
    /* saveComputedGradients(eik_g6, pathSaveGrads); */
    /* saveTimes(time_taken, pathTimes); */

    /* printGeneralInfo(eik_g6); */

    /* eik_grid_dealloc(&eik_g6); */




    /* pathPoints = "./H6/H6_MeshPoints.txt"; */
    /* pathFaces = "./H6/H6_Faces.txt"; */
    /* pathEdges = "./H6/H6_Edges.txt"; */
    /* pathEdgesInFace = "./H6/H6_EdgesInFace.txt"; */
    /* pathNeighbors = "./H6/H6_Neigh.txt"; */
    /* pathIncidentFaces = "./H6/H6_IncidentFaces.txt"; */
    /* pathIndices = "./H6/H6_Indices.txt"; */
    /* pathIndices = "./H6/H6_Indices.txt"; */
    /* pathBoundary = "./H6/H6_BoundaryCurve.txt"; */

    /* pathSaveEiks = "./H6/H6_ComputedValuesFast.bin"; */
    /* pathSaveGrads = "./H6/H6_ComputedGradientsFast.bin"; */
    /* pathTimes = "./H6/H6_TimesFast.bin"; */

    /* eik_gridS *eik_g7; */
    /* eik_grid_alloc(&eik_g7); */

    /* eik_grid_initFromFile(eik_g7, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace, */
    /* 			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary); */

    /* // printAllInfoMesh(eik_g1); */

    /* start_t = clock(); */
    /* marcher_T2(eik_g7, rBall);; */
    /* end_t = clock(); */
    /* time = (double)(end_t - start_t)/ CLOCKS_PER_SEC; */
    /* time_taken[0] = time; */
    

    /* saveComputedValues(eik_g7, pathSaveEiks); */
    /* saveComputedGradients(eik_g7, pathSaveGrads); */
    /* saveTimes(time_taken, pathTimes); */

    /* printGeneralInfo(eik_g7); */

    /* eik_grid_dealloc(&eik_g7); */


  

}
