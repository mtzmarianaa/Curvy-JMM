/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"
// #include "opti_method.h" // currently using Python for the optimization
#include "linAlg.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <json-c/json.h> // used for reading the json type string from python


// things to call out python optimizer
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/wait.h>

void eik_grid_alloc(eik_gridS **eik_g ) {
  *eik_g = malloc(sizeof(eik_gridS));
  assert(*eik_g != NULL);
}

void eik_grid_dealloc(eik_gridS **eik_g ) {
  free(*eik_g);
  *eik_g = NULL;
}

void fanUpdate_alloc(fanUpdateS **fanUpdate) {
  *fanUpdate = malloc(sizeof(fanUpdateS));
  assert(*fanUpdate != NULL);
}

void fanUpdate_dalloc(fanUpdateS **fanUpdate) {
  free(*fanUpdate);
  *fanUpdate = NULL;
}

void eik_grid_init( eik_gridS *eik_g, size_t *start, size_t nStart, mesh2S *mesh2) {
  // the rest of the parameters, eik_vals, p_queueG, current_states are going to be assigned inside
  eik_g->start = start;
  eik_g->nStart = nStart;
  eik_g->mesh2 = mesh2;

  // we first set all the current eik_vals to infinity, set all the current_states to 0 (far)
  double *eik_vals;
  double (*eik_grad)[2]; // this is a pointer to a list of the gradients of the eikonal
  size_t *current_states;
  fanUpdateS *fanUpdates;
  
  eik_vals = malloc(mesh2->nPoints*sizeof(double)); 
  current_states = malloc(mesh2->nPoints*sizeof(int));
  eik_grad = malloc(2*mesh2->nPoints*sizeof(double)); // each gradient has two coordinates (for each point)
  
  for(int i = 0; i<mesh2->nPoints; i++){
    eik_vals[i] = INFINITY; // set them all to infinity
    current_states[i] = 0; // set them all to far
    eik_grad[i][0] = 0;
    eik_grad[i][1] = 0;
  }
  
  eik_g->eik_vals = eik_vals;
  eik_g->current_states = current_states;
  eik_g->eik_grad = eik_grad;

  // we initialize the priority queue, all the elements in start are going to be inserted and their current states set to 1
  // notice that we need to add both the index of the starting points AND their eikonal value (which is 0) to the p_queue struct
  p_queue *p_queueG;
  priority_queue_alloc(&p_queueG); // allocate
  priority_queue_init(p_queueG); // initiate
  for(int i = 0; i<nStart; i++){
    insert(p_queueG, 0, start[i]); // insert all the starting points with eikonal value 0
  }
  eik_g->p_queueG = p_queueG;
  assert(&eik_g != NULL); // eik_g should not be null
}


void fanUpdate_init(fanUpdateS *fanUpdate, triangleFanS *triFan, double *params,
		    double T0, double grad0[2], double T1, double grad1[2],
		    size_t nIndCrTop, size_t *indCrTop, double *paramsCrTop,
		    size_t nIndStTop, size_t *indStTop, double *paramsStTop,
		    double THat, double (*grads)[2], double (*path)[2],
		    double gradHat[2]) {
  // "manually" initialize a triangle fan update
  fanUpdate->triFan = triFan;
  fanUpdate->params = params;
  fanUpdate->T0 = T0;
  fanUpdate->grad0[0] = grad0[0];
  fanUpdate->grad0[1] = grad0[1];
  fanUpdate->T1 = T1;
  fanUpdate->grad1[0] = grad1[0];
  fanUpdate->grad1[1] = grad1[1];
  fanUpdate->nIndCrTop = nIndCrTop;
  fanUpdate->indCrTop = indCrTop;
  fanUpdate->paramsCrTop = paramsCrTop;
  fanUpdate->nIndStTop = nIndStTop;
  fanUpdate->indStTop = indStTop;
  fanUpdate->paramsStTop = paramsStTop;
  fanUpdate->THat = THat;
  fanUpdate->grads = grads;
  fanUpdate->path = path;
  fanUpdate->gradHat = gradHat;
}


void fanUpdate_initPreOpti(fanUpdateS *fanUpdate, triangleFanS *triFan, double T0,
		    double grad0[2], double T1, double grad1[2]) {
  // Init a triangle fan before using this information in the optimization
  fanUpdate->triFan = triFan;
  fanUpdate->T0 = T0;
  fanUpdate->grad0[0] = grad0[0];
  fanUpdate->grad0[1] = grad0[1];
  fanUpdate->T1 = T1;
  fanUpdate->grad1[0] = grad1[0];
  fanUpdate->grad1[1] = grad1[1];
}

void eik_grid_initFromFile(eik_gridS *eik_g, size_t *start, size_t nStart, char const *pathPoints, char const *pathFaces,
			    char const *pathEdges, char const *pathEdgesInFace,
			    char const *pathNeighbors, char const *pathIncidentFaces,
			    char const *pathIndices, char const *pathBoundary) {
  // the only difference between this and the previous method is that in here we do need to initialize the Mesh structure
  mesh2S *mesh2;
  mesh2_init_from_meshpy(mesh2, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			 pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);
  // and then we can use the previous method
  eik_grid_init( eik_g, start, nStart, mesh2); // voila
}

void printGeneralInfo(eik_gridS *eik_g) {
  printf("\n\n\n\n     GENERAL INFORMATION ON THIS EIKONAL STRUCT     \n\n");
  printf("Number of starting points: %zu \n", eik_g->nStart);
  printf("The starting points are: \n");
  for(int i = 0; i<eik_g->nStart; i++) {
    printf("|   %zu   |", eik_g->start[i]);
  }
  printGeneralInfoMesh(eik_g->mesh2);
  printf("Current state of priority queue: \n");
  printeik_queue(eik_g->p_queueG);
  printf("\nCurrent Eikonal values: \n");
  for(int i = 0; i<eik_g->mesh2->nPoints; i++){
    double x[2];
    x[0] = eik_g->mesh2->points[i][0];
    x[1] = eik_g->mesh2->points[i][1];
    printf("Index   %d    ||  Coordinates:   (%fl, %fl)    ||  Eikonal:   %fl     ||  Current state:   %zu \n", i, x[0], x[1]  , eik_g->eik_vals[i] , eik_g->current_states[i]);
  }
}

void printInfoFanUpdate(eik_gridS *eik_g, size_t k) {
  // prints all the information from a triangle fan update
  triangleFanS *triFan;
  triFan = eik_g->fanUpdate[k].triFan;
  printf("\n\n\n  PRINTING INFORMATION FROM THE %zu-th TRIANGLE FAN UPDATE\n", k);
  printEverythingTriFan(triFan);
  int i;
  size_t nReg = triFan->nRegions;
  printf("\nParams: \n");
  for( i = 0; i<(2*nReg + 1); i++){
    printf(" %fl  ", eik_g->fanUpdate[k].params[i]);
  }

  printf("\nT0: %fl", eik_g->fanUpdate[k].T0);
  printf("\ngrad0:  %fl  |  %fl", eik_g->fanUpdate[k].grad0[0], eik_g->fanUpdate[k].grad0[1]);
  printf("\nT1: %fl", eik_g->fanUpdate[k].T1);
  printf("\ngrad1:  %fl  |  %fl", eik_g->fanUpdate[k].grad1[0], eik_g->fanUpdate[k].grad1[1]);
  
}

void printAllInfoMesh(eik_gridS *eik_g){
  printEverythingInMesh(eik_g->mesh2);
}



/* void initializePointsNear(eik_gridS *eik_g, double rBall) { */
/*   // THIS JUST SETS THE CURRENT STATE TO VALID AND ADDS THE TRUE EIKONAL */
/*   // given a ball of radius rBall around the initial points we initialize all the points inside those balls with the */
/*   // true value of the eikonal (i.e. the distance times the index of refraction). We are going to assume that */
/*   // all those balls belong to the same regions of indices of refraction */
/*   double xMinxStart[2], xStart[2], xCurrent[2], normCurrent, initialIndexRefraction; */
/*   int indexStart; */
/*   for(int j = 0; j<eik_g->nStart; j++){ */
/*     deleteRoot(eik_g->p_queueG); */
/*     indexStart = eik_g->start[j]; */
/*     xStart[0] = eik_g->mesh2->points[ indexStart ][0]; */
/*     xStart[1] = eik_g->mesh2->points[ indexStart ][1]; */
/*     for(int i = 0; i<eik_g->mesh2->nPoints; i++){ */
/*       xCurrent[0] = eik_g->mesh2->points[i][0]; */
/*       xCurrent[1] = eik_g->mesh2->points[i][1]; */
/*       vec2_subtraction( xCurrent, xStart, xMinxStart ); */
/*       normCurrent = l2norm(xMinxStart); */
/*       if(normCurrent < rBall ){ */
/*         if( eik_g->current_states[i] == 1 ){ */
/*           // if it was previously considered as trial we need to delete this from the queue directly */
/*           delete_findIndex(eik_g->p_queueG, i); */
/*         } */
/*         // if this happens, this point is "close enough" to the starting point so that we can initialize it */
/* 	initialIndexRefraction = minEtaFromTwoPoints(eik_g->mesh2, i, indexStart); */
/*         eik_g->current_states[i] = 2; // we initialized it directly */
/*         eik_g->eik_vals[i] = initialIndexRefraction*normCurrent; // add their true Eikonal value */
/*         eik_g->eik_grad[i][0] = initialIndexRefraction*xMinxStart[0]/normCurrent; // update its gradient */
/*         eik_g->eik_grad[i][1] = initialIndexRefraction*xMinxStart[1]/normCurrent; */
/*         addNeighbors_fromAccepted(eik_g, i); // we add its neighbors */
/*       } */
/*     } */
/*   } */

/* } */


void approximateEikonalGradient(double xA[2], double xB[2], double xHat[2], double parameterization, double indexRefraction, double grad[2]) {
  // given two points in a base with its corresponding parameterization we approximate the gradient of the Eikonal
  // depending on the type of update performed on xHat xA is x0 and xB could be either x1 or x2
  double MinParametrization, xBasePart1[2], xBasePart2[2], xBase[2], normBase, direction[2];
  MinParametrization = 1 - parameterization;
  scalar_times_2vec( MinParametrization, xA, xBasePart1 );
  scalar_times_2vec( parameterization, xB, xBasePart2 );
  vec2_addition( xBasePart1, xBasePart2, xBase );
  vec2_subtraction(xHat, xBase, direction);
  normBase = l2norm(direction);
  grad[0] = indexRefraction*direction[0]/normBase;
  grad[1] = indexRefraction*direction[1]/normBase;
}

void findEdgesOnValidFront(eik_gridS *eik_g, size_t index0, int indices1[2], int firstTriangles[2]) {
  // given a newly accepted node with index index0 go around its neighbors
  // and find (at most 2) indices index1 such that the edge x0 x1 is on the valid fron
  // notice that indices1 is an int because it's going to have a -1 if we just have one valid edge from x0
  assert( eik_g->current_states[index0] == 2); // this HAS to be a valid node
  int nNeis, j, i;
  size_t thisNeighbor, previousNeighbor, thisTriangle;
  size_t possibleTriangles[2], possibleThirdVertices[2];
  j = 0;
  i = 1;
  indices1[0] = -1;
  indices1[1] = -1; // no neighbors found yet
  firstTriangles[0] = -1;
  firstTriangles[1] = -1; // no triangles found yet
  nNeis = eik_g->mesh2->neighbors[index0].len; // get the number of neighbors of index0
  // just start with any neighbor (doesn't really matter which one)
  thisNeighbor = (size_t)eik_g->mesh2->neighbors[index0].neis_i[0];
  previousNeighbor = index0;
  // go around
  while( j < 2 & i < nNeis ) {
    // either we find two valid edges or we go around all the possible neighbors
    twoTrianglesFromEdge(eik_g->mesh2, index0, thisNeighbor, possibleTriangles, possibleThirdVertices );
    if( possibleThirdVertices[0] != previousNeighbor ) {
      // the next neighbor is on index 0
      previousNeighbor = thisNeighbor;
      thisNeighbor = possibleThirdVertices[0];
      thisTriangle = possibleTriangles[0];
    }
    else{
      // the next neighbor is on index 1
      previousNeighbor = thisNeighbor;
      thisNeighbor = possibleThirdVertices[1];
      thisTriangle = possibleTriangles[1];
    }
    // after we determined which is the neighbor that we need to cycle next
    if( eik_g->current_states[previousNeighbor] == 2 & eik_g->current_states[thisNeighbor] != 2 ){
      // thisNeighbor is part of the valid front
      indices1[j] = (int)previousNeighbor;
      firstTriangles[j] = (int)thisTriangle;
      j ++;
    }
    i++;
  }
  // figure out if we actually have at least one neighbor
  assert( indices1[0] != -1 | indices1[1] != -1);
}




void initTriFan(eik_gridS *eik_g, triangleFanS *triFan,
		size_t index0, size_t index1, size_t index2,
		size_t indexHat, size_t firstTriangle, double angleMax) {
  // after running twoTrianglesFromEdge select one and initialize a triangle fan like this
  // angle max is the biggest angle on xk1 x0 xk inside the trianlge fan
  // this is going to be useful to know if this triangle fan is feasible or not
  // given two valid indices index0, index1 and a direction index2
  // we set up a triangle fan triFan that updates from the edge index0 index1
  // in the direction of the firstTriangle to update indexHat
  double pi;
  pi = acos(-1.0);
  size_t nRegions, *listFaces, *listEdges, *listIndicesNodes;
  double x0[2], x1[2], xHat[2];
  x0[0] = eik_g->mesh2->points[index0][0];
  x0[1] = eik_g->mesh2->points[index0][1];
  x1[0] = eik_g->mesh2->points[index1][0];
  x1[1] = eik_g->mesh2->points[index1][1];
  xHat[0] = eik_g->mesh2->points[indexHat][0];
  xHat[1] = eik_g->mesh2->points[indexHat][1];
  //////////////////////////
  // first we need to count the number of regions
  double etakM1, etak, xk[2], xk1[2], thisAngle, angleRegion;
  size_t possibleTriangles[2], possibleThirdVertices[2], indexk, indexk1;
  size_t thisTriangle, prevTriangle;
  nRegions = 1;
  indexk = index1;
  indexk1 = index2;
  xk[0] = x1[0];
  xk[1] = x1[1];
  xk1[0] =  eik_g->mesh2->points[index2][0];
  xk1[1] =  eik_g->mesh2->points[index2][1];
  angleMax = angleThreePoints(xk, x0, xk1); // first angle
  thisAngle = angleMax;
  angleRegion = thisAngle;
  thisTriangle = firstTriangle;
  etak = eik_g->mesh2->eta[firstTriangle];
  while( indexk1 != indexHat ){
    // circle around and see what we get
    twoTrianglesFromEdge(eik_g->mesh2, index0, indexk1, possibleTriangles, possibleThirdVertices);
    if( thisTriangle != possibleTriangles[0]){
      // the next triangle is possibleTriangles[0]
      indexk = indexk1;
      indexk1 = possibleThirdVertices[0];
      prevTriangle = thisTriangle;
      thisTriangle = possibleTriangles[0];
    }
    else {
      // the next triangle is possibleTriangles[1]
      indexk = indexk1;
      indexk1 = possibleThirdVertices[1];
      prevTriangle = thisTriangle;
      thisTriangle = possibleTriangles[1];
    }
    // update
    xk[0] = xk1[0];
    xk[1] = xk1[1];
    xk1[0] = eik_g->mesh2->points[indexk1][0];
    xk1[1] = eik_g->mesh2->points[indexk1][1];
    etakM1 = etak;
    etak = eik_g->mesh2->eta[thisTriangle];
    thisAngle = angleThreePoints(xk, x0, xk1);
    // angle from one region to another one
    if( etakM1 != etak ){
      angleRegion = 0;
    }
    else{
      angleRegion = angleRegion + thisAngle;
    }
    // but we are looking for the maximum angle of change in regions
    if( angleRegion > angleMax ){
      angleMax = thisAngle;
    }
    if( angleMax > pi ){
      // we can't  update here
      return;
    }
    nRegions ++; // add a new region
  }
  //////////////////////////
  // with the number of regions set we can go around and malloc everything
  // notice that in mesh2D.c we have a function called
  // triangleFan_initFromIndices and for this we only need the list of indices
  listIndicesNodes = malloc((nRegions + 2)*sizeof(size_t));
  int i = 2;
  indexk = index1;
  indexk1 = index2;
  thisTriangle = firstTriangle;
  listIndicesNodes[0] = indexk;
  listIndicesNodes[1] =  indexk1;
  while( indexk1 != indexHat ) {
    // circle around and see what we get
    twoTrianglesFromEdge(eik_g->mesh2, index0, indexk1, possibleTriangles, possibleThirdVertices);
    if( thisTriangle != possibleTriangles[0] ) {
      // the next triangle is possibleTriangles[0]
      indexk = indexk1;
      indexk1 = possibleThirdVertices[0];
      prevTriangle = thisTriangle;
      thisTriangle = possibleTriangles[0];
    }
    else {
      // the next triangle is possibleTriangles[1]
      indexk = indexk1;
      indexk1 = possibleThirdVertices[1];
      prevTriangle = thisTriangle;
      thisTriangle = possibleTriangles[1];
    }
    listIndicesNodes[i] = indexk1;
    i ++;
    
  }
  //////////////////////////
  // now we have the list of indices, we can initialize the triangle
  triangleFan_initFromIndices(triFan, eik_g->mesh2, nRegions, index0,
			      index1, indexHat, listIndicesNodes);
}


void createJSONinput(fanUpdateS *fanUpdate, char *input_json) {
  // creates a JSON object to represent input data from a triangle fan
  char temp_buffer[2500]; // temporary buffer
  char num[20]; // temporary char where we are going to store our numbers (hopefully they'll fit)
  char comma[5], finArr[5], extra[20], sqIn[5], sqFin[5]; // extra things like commas, [, ] ....
  int i;
  comma = ", ";
  finArr = "], ";
  sqIn = "[";
  sqFin = "]";
  // add x0
  temp_buffer = "{\"x0\": [";
  strcat(temp_buffer, extra);
  snprintf(num, 50, "%f", fanUpdate->triFan->x0[0]); // chartify x0[0]
  strcat(temp_buffer, num); // concatenate x0[0] to the temporary buffer
  strcat(temp_buffer, comma);
  snprintf(num, 50, "%f", fanUpdate->triFan->x0[1]); // chartify x0[1]
  strcat(temp_buffer, num); // concatenate x0[1] to the temporary buffer
  strcat(temp_buffer, finArr); // finish "x0" : [x0[0], x0[1]]

  // add T0
  extra = "\"T0\": ";
  strcat(temp_buffer, extra);
  snprintf(num, 50, "%f", funUpdate->T0); // chartify T0
  strcat(temp_buffer, num);
  strcat(temp_buffer, comma);

  // add grad0
  extra = "\"grad0\": [";
  strcat(temp_buffer, extra);
  snprintf(num, 50, "%f", fanUpdate->grad0[0]); 
  strcat(temp_buffer, num); 
  strcat(temp_buffer, comma);
  snprintf(num, 50, "%f", fanUpdate->grad0[1]);
  strcat(temp_buffer, num);
  strcat(temp_buffer, finArr);
  
  // add x1
  extra = "\"x1\": [";
  strcat(temp_buffer, extra);
  snprintf(num, 50, "%f", fanUpdate->triFan->x1[0]); 
  strcat(temp_buffer, num); 
  strcat(temp_buffer, comma);
  snprintf(num, 50, "%f", fanUpdate->triFan->x1[1]);
  strcat(temp_buffer, num);
  strcat(temp_buffer, finArr);

  // add T1
  extra = "\"T1\": ";
  strcat(temp_buffer, extra);
  snprintf(num, 50, "%f", funUpdate->T1); // chartify T1
  strcat(temp_buffer, num);
  strcat(temp_buffer, comma);

  // add grad1
  extra = "\"grad0\": [";
  strcat(temp_buffer, extra);
  snprintf(num, 50, "%f", fanUpdate->grad1[0]); 
  strcat(temp_buffer, num); 
  strcat(temp_buffer, comma);
  snprintf(num, 50, "%f", fanUpdate->grad1[1]);
  strcat(temp_buffer, num);
  strcat(temp_buffer, finArr);

  // add xHat
  extra = "\"xHat\": [";
  strcat(temp_buffer, extra);
  snprintf(num, 50, "%f", fanUpdate->triFan->xHat[0]);
  strcat(temp_buffer, num);
  strcat(temp_buffer, comma);
  snprintf(num, 50, "%f", fanUpdate->triFan->xHat[1]); 
  strcat(temp_buffer, num);
  strcat(temp_buffer, finArr);

  // use a for loop to add the list of indices
  extra = "\"listIndices\": [";
  strcat(temp_buffer, extra);
  for( i = 0; i<2*fanUpdate->triFan->nRegions + 1; i++){
    snprintf(num, 50, "%f", fanUpdate->triFan->listIndices[i]); // chartify the i-th index
    strcat(temp_buffer, num);
    if( i < 2*fanUpdate->triFan->nRegions){
      strcat(temp_buffer, comma); // add the comma
    }
  }
  strcat(temp_buffer, finArr);


  // use a for loop to add the list of xk
  extra = "\"listxk\": [";
  strcat(temp_buffer, extra);
  for( i = 0; i<fanUpdate->triFan->nRegions + 2; i++){
    strcat(temp_buffer, sqIn); // add [
    snprintf(num, 50, "%f", fanUpdate->triFan->listxk[i][0]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, comma);
    snprintf(num, 50, "%f", fanUpdate->triFan->listxk[i][1]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, sqFin); // add ]
    if( i <fanUpdate->triFan->nRegions + 1){
      strcat(temp_buffer, comma); // add the comma
    }
  }
  strcat(temp_buffer, finArr);


  // use a for loop to add the list of B0k
  extra = "\"listB0k\": [";
  strcat(temp_buffer, extra);
  for( i = 0; i<fanUpdate->triFan->nRegions + 1; i++){
    strcat(temp_buffer, sqIn); // add [
    snprintf(num, 50, "%f", fanUpdate->triFan->listB0k[i][0]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, comma);
    snprintf(num, 50, "%f", fanUpdate->triFan->listB0k[i][1]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, sqFin); // add ]
    if( i <fanUpdate->triFan->nRegions ){
      strcat(temp_buffer, comma); // add the comma
    }
  }
  strcat(temp_buffer, finArr);


  // use a for loop to add the list of Bk
  extra = "\"listBk\": [";
  strcat(temp_buffer, extra);
  for( i = 0; i<fanUpdate->triFan->nRegions + 1; i++){
    strcat(temp_buffer, sqIn); // add [
    snprintf(num, 50, "%f", fanUpdate->triFan->listBk[i][0]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, comma);
    snprintf(num, 50, "%f", fanUpdate->triFan->listBk[i][1]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, sqFin); // add ]
    if( i <fanUpdate->triFan->nRegions ){
      strcat(temp_buffer, comma); // add the comma
    }
  }
  strcat(temp_buffer, finArr);


  // use a for loop to add the list of BkBk1
  extra = "\"listBkBk1\": [";
  strcat(temp_buffer, extra);
  for( i = 0; i<2*fanUpdate->triFan->nRegions ; i++){
    strcat(temp_buffer, sqIn); // add [
    snprintf(num, 50, "%f", fanUpdate->triFan->listBkBk1[i][0]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, comma);
    snprintf(num, 50, "%f", fanUpdate->triFan->listBkBk1[i][1]); 
    strcat(temp_buffer, num);
    strcat(temp_buffer, sqFin); // add ]
    if( i <2*fanUpdate->triFan->nRegions -1 ){
      strcat(temp_buffer, comma); // add the comma
    }
  }
  strcat(temp_buffer, finArr);


  // OPTIONS FOR PLOTTING CHANGE THIS ACCORDINGLY
  extra = "\"plotBefore\": 1, \"plotAfter\": 1, \"plotOpti\": 1 }";
  strcat(temp_buffer, extra);
  
  printf(temp_buffer); // just to see what we get

  // our point should point to temp_buffer[0]
  input_json = &temp_buffer[0];

  
}



void optimizeTriangleFan_wPython(fanUpdateS *fanUpdate){
  // using pipes and a lot of fancy methos call python from here
  // and optimize the triangle fan in fanUpdate, save all info
  int pipefd[2]; // pipe
  pid_t pid; // proces id

  if (pipe(pipefd) == -1) {
    printf("\nProblem when opening pipe");
    exit(EXIT_FAILURE);
  }

  pid = fork(); // in the child process we execute python, in the parent process we write and read
  if (pid == -1) {
    perror("fork");
    exit(EXIT_FAILURE);
  }

  if (pid == 0) { // child process
    if (dup2(pipefd[1], STDOUT_FILENO) == -1) {
      printf("\nProblem with child process");
      exit(EXIT_FAILURE);
    }

    close(pipefd[0]); // CLOSE 0

    // execute Python script
    execlp("python3", "python3", "stepWithPython.py", (char *) NULL);
  }
  else { // parent process
    close(pipefd[1]); // CLOSE 1

    // create JSON object to represent input data
    char* input_json; //
    createJSONinput(fanUpdate, input_json); // put everything in the json order we want

    // write input data to pipe
    if (write(pipefd[0], input_json, strlen(input_json)) == -1) {
      printf("\nProblem when writting triInfo to run Python\n");
      exit(EXIT_FAILURE);
    }

    close(pipefd[0]); // CLOSE 0, we already wrote the input to the python optimizer

    // wait for child process to complete
    wait(NULL);

    // read output data from pipe
    int fd = open("output.json", O_RDONLY);
    if (fd == -1) {
      printf("\nProblem when opening the output json from Python\n");
      exit(EXIT_FAILURE);
    }

    char buffer[1024];
    ssize_t num_read = read(fd, buffer, sizeof(buffer));
    if (num_read == -1) {
      perror("read");
      exit(EXIT_FAILURE);
    }

    close(fd); // CLOSE PIPE

    // PROCESS ALL THE INFORMATION FROM THE JSON FILE AND PARSE IT ACCORDINGLY

      // deserialize output data
      json_object *output_obj = json_tokener_parse(buffer);
      if (output_obj == NULL) {
          fprintf(stderr, "Error parsing JSON: %s\n", json_tokener_error_desc(json_tokener_get_error(json_tokener_new())));
          exit(EXIT_FAILURE);
      }

      // extract array of doubles from output object
      json_object *output_array = json_object_object_get(output_obj, "result");
      if (output_array == NULL || !json_object_is_type(output_array, json_type_array)) {
          fprintf(stderr, "Error extracting array from JSON output\n");
          exit(EXIT_FAILURE);
      }
      int num_elements = json_object_array_length(output_array);
      double output_values[num_elements];
      for (int i = 0; i < num_elements; i++) {
          json_object *value_obj = json_object_array_get_idx(output_array, i);
          if (value_obj == NULL || !json_object_is_type(value_obj, json_type_double)) {
              fprintf(stderr, "Error extracting double value from JSON array\n");
              exit(EXIT_FAILURE);
          }
          output_values[i] = json_object_get_double(value_obj);
      }

}



void addNeighbors_fromAccepted(eik_gridS *eik_g, size_t minIndex) {
  // given a recently accepted node with index minIndex we update its neighbors
  // using triangle fans and up to two valid edges
}





/* void popAddNeighbors(eik_gridS *eik_g) { */
/*   // int nNeighs; */
/*   int minIndex = indexRoot(eik_g->p_queueG); */
/*   deleteRoot(eik_g->p_queueG); // delete the root from the priority queue */
/*   eik_g->current_states[minIndex] = 2; // set the newly accepted index to valid */
/*   addNeighbors_fromAccepted(eik_g, minIndex); // add neighbors from the recently accepted index */
/* } */

int currentMinIndex(eik_gridS *eik_g) {
  return indexRoot(eik_g->p_queueG);
}

int nStillInQueue(eik_gridS *eik_g) {
  return getSize(eik_g->p_queueG);
}

void saveComputedValues(eik_gridS *eik_g, const char *pathFile) {
  FILE *fp;
  fp = fopen(pathFile, "wb");
  fwrite(eik_g->eik_vals, sizeof(double), eik_g->mesh2->nPoints, fp);
  fclose(fp);
}

void saveComputedGradients(eik_gridS *eik_g, const char *pathFile) {
    FILE *fp;
    fp = fopen(pathFile, "wb");
    
    for (int i=0; i<eik_g->mesh2->nPoints; ++i){
      fwrite(eik_g->eik_grad[i], sizeof(double), 2, fp);
    }

    fclose(fp);
}




