#include "mesh2D.h"
#include "eik_grid.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(){

    eik_gridS *eik_g;
    eik_grid_alloc(&eik_g);
    const char *pathPoints, *pathFaces, *pathEdges, *pathEdgesInFace, *pathNeighbors;
    const char *pathIncidentFaces, *pathIndices, *pathBoundary;
    
    pathPoints = "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow//H0/H0_MeshPoints.txt";
    pathFaces =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Faces.txt";
    pathEdges =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Edges.txt";
    pathEdgesInFace = "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_EdgesInFace.txt";
    pathNeighbors =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Neigh.txt";
    pathIncidentFaces =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_IncidentFaces.txt";
    pathIndices =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Indices.txt";
    pathIndices =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_Indices.txt";
    pathBoundary =  "/Users/marianamartinez/Documents/Curvy-JMM/TestBaseSnow/H0/H0_BoundaryCurve.txt";


    size_t *start, start_int, nStart;

    start_int = 0;
    start = &start_int;
    nStart = 1;


    eik_grid_initFromFile(eik_g, start, nStart, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			  pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);

    printf("Points 5: %f %f\n", eik_g->mesh2->points[5][0], eik_g->mesh2->points[5][1]);

    printAllInfoMesh(eik_g);


    //printGeneralInfo(eik_g);


    // test the initializing points near

    printGeneralInfo(eik_g);

    double rBall = 2.0;

    printf("Initialize Points Near \n\n\n");

    initializePointsNear(eik_g, rBall);  // there is an error here

    printGeneralInfo(eik_g);

    size_t testIndex;
    testIndex = 2521;

    // test findEdgesOnValidFront
    int indices1[2], indices2[2], firstTriangles[2];

    //findEdgesOnValidFront(eik_g, testIndex, indices1, indices2, firstTriangles);


    eik_grid_dealloc(&eik_g);


    mesh2S *mesh2;
    mesh2_alloc(&mesh2);

    
    mesh2_init_from_meshpy(mesh2, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			   pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);
    printf("Points 5 just mesh2  %f  %f\n", mesh2->points[5][0], mesh2->points[5][1]);
    

    printf("GENERAL INFO \n\n");
    printGeneralInfoMesh(mesh2);

    printf("\n\n\n\nALL INFO \n\n");
    printEverythingInMesh(mesh2);


    mesh2_dealloc(&mesh2);


    /* // Now we test the triangleFan */

    /* triangleFanS *triFan; */
    /* triangleFan_alloc(&triFan); */
    /* size_t nRegions, index0, index1, indexHat, *listIndicesNodes; */

    /* // going to use: */
    /* //      triangle 247 with points 89, 350, 513 */
    /* //      triangle 1017 with points 89, 90, 513 */

    /* nRegions = 2; */
    /* index0 = 89; // triangle 247 with points 89, 350, 513 */
    /* index1 = 350; */
    /* indexHat = 90; */

    /* listIndicesNodes = malloc(4*sizeof(size_t)); */
    /* listIndicesNodes[0] = 89; */
    /* listIndicesNodes[1] = 350; */
    /* listIndicesNodes[2] = 513; */
    /* listIndicesNodes[3] = 90; */

    /* triangleFan_initFromIndices(triFan, mesh2, nRegions, index0, index1, indexHat, listIndicesNodes); */

    /* printEverythingTriFan(triFan); */

    /* mesh2_dealloc(&mesh2); */

  // init a triangle fan

  /* triangleFanS *triFan; */
  /* size_t nRegions; */
  /* double x0[2], x1[2], xHat[2], *listIndices; */
  /* double (*listxk)[2], (*listBk)[2], (*listB0k)[2], (*listBkBk1)[2]; */
  /* int *listFaces, *listEdges; */
  
  /* triangleFan_alloc(&triFan); */

  /* nRegions = 2; */
  /* x0[0] = 9.01006524; */
  /* x0[1] = -4.739905; */
  /* x1[0] = 8.91006524; */
  /* x1[1] = -4.539905; */
  /* xHat[0] = 9.53879533; */
  /* xHat[1] = -3.97683432; */

  /* listIndices = malloc(5*sizeof(double)); */
  /* listIndices[0] = 1.0; */
  /* listIndices[1] = 1.0; */
  /* listIndices[2] = 1.0; */
  /* listIndices[3] = 1.452; */
  /* listIndices[4] = 1.0; */

  /* listxk = malloc(8*sizeof(double)); */
  /* listxk[0][0] =  9.01006524; */
  /* listxk[0][1] = -4.739905; */
  /* listxk[1][0] = 8.91006524; */
  /* listxk[1][1] = -4.539905; */
  /* listxk[2][0] = 9.23879533; */
  /* listxk[2][1] = -3.82683432; */
  /* listxk[3][0] = 9.53879533; */
  /* listxk[3][1] = -3.97683432; */

  /* listB0k = malloc(6*sizeof(double)); */
  /* listB0k[0][0] = -0.1; */
  /* listB0k[0][1] = 0.2; */
  /* listB0k[1][0] = 0.22873008; */
  /* listB0k[1][1] = 0.91307067; */
  /* listB0k[2][0] = 0.52873008; */
  /* listB0k[2][1] = 0.76307067; */

  /* listBk = malloc(6*sizeof(double)); */
  /* listBk[0][0] = -0.1; */
  /* listBk[0][1] = 0.2; */
  /* listBk[1][0] = 0.22873008; */
  /* listBk[1][1] = 0.91307067; */
  /* listBk[2][0] = 0.52873008; */
  /* listBk[2][1] = 0.76307067; */


  /* listBkBk1 = malloc(8*sizeof(double)); */
  /* listBkBk1[0][0] = 0.4022869; */
  /* listBkBk1[0][1] = 0.7895325; */
  /* listBkBk1[1][0] = 0.33910078; */
  /* listBkBk1[1][1] = 0.8186617; */
  /* listBkBk1[2][0] = 0.3; */
  /* listBkBk1[2][1] = -0.15; */
  /* listBkBk1[3][0] = 0.3; */
  /* listBkBk1[3][1] = -0.15; */

  /* listEdges = malloc(5*sizeof(int)); */
  /* listEdges[0] = -1; */
  /* listEdges[1] = -1; */
  /* listEdges[2] = -1; */
  /* listEdges[3] = -1; */
  /* listEdges[4] = -1; */

  /* listFaces = malloc(2*sizeof(int)); */
  /* listFaces[0] = -1; */
  /* listFaces[1] = -1; */

  


  /* triangleFan_init(triFan, nRegions, x0, x1, xHat, */
  /* 		   listFaces, listIndices, listEdges, */
  /* 		   listxk, listB0k, listBk, listBkBk1); */

  /* // print */

  /* printEverythingTriFan(triFan); */

  /* // build the triangle update */
  /* fanUpdateS *fanUpdate; */
  /* fanUpdate_alloc(&fanUpdate); */

  /* double T0, grad0[2], T1, grad1[2]; */

  /* T0 = 25.887855886616833; */
  /* grad0[0] = 0.62334049; */
  /* grad0[1] = 0.78195053; */
  /* T1 = 25.995574287564278; */
  /* grad1[0] = 0.4539905; */
  /* grad1[1] = 0.89100652; */

  /* fanUpdate_initPreOpti(fanUpdate, triFan, T0, grad0, T1, grad1); */

  /* printf("\nStart opti\n"); */

  /* //createJSONFile(fanUpdate, "/Users/marianamartinez/Documents/Curvy-JMM/JMM/update.json"); */

  /* char buffer[5000]; */
  /* int fd = open("/Users/marianamartinez/Documents/Curvy-JMM/JMM/update.json", O_RDONLY); */
  /* ssize_t num_read = read(fd, buffer, sizeof(buffer)); */
  /* json_object *output_obj = json_tokener_parse(buffer); */
  /* deserializeJSONoutput(fanUpdate, output_obj); */
  
  

  /* optimizeTriangleFan_wPython(fanUpdate); */

  /* printInfoFanUpdate(fanUpdate); */

  /* fanUpdate_dalloc(&fanUpdate); */
  /* triangleFan_dalloc(&triFan); */

}
