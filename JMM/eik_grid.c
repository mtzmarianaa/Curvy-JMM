/* EIKONAL GRID

This is the Eikonal grid with different specifications

*/

#include "eik_grid.h"
#include "priority_queue.h"
// #include "opti_method.h" // currently using Python for the optimization
#include "linAlg.h"
#include "files_methods.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

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
  //eik_g->mesh2 = mesh2;

  // we first set all the current eik_vals to infinity, set all the current_states to 0 (far)
  double *eik_vals;
  double (*grads)[2]; // this is a pointer to a list of the gradients of the eikonal
  size_t *current_states;
  
  eik_vals = malloc(mesh2->nPoints*sizeof(double)); 
  current_states = malloc(mesh2->nPoints*sizeof(size_t));
  grads = malloc(2*mesh2->nPoints*sizeof(double));

  
  for(int i = 0; i<mesh2->nPoints; i++){
    eik_vals[i] = INFINITY; // set them all to infinity
    current_states[i] = 0; // set them all to far
    grads[i][0] = 0;
    grads[i][1] = 0;
  }
  
  eik_g->eik_vals = eik_vals;
  eik_g->current_states = current_states;
  eik_g->grads = grads;

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
  fanUpdate->gradHat[0] = gradHat[0];
  fanUpdate->gradHat[1] = gradHat[1];
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
  mesh2_alloc(&mesh2);
  mesh2_init_from_meshpy(mesh2, pathPoints, pathFaces, pathEdges, pathEdgesInFace,
			 pathNeighbors, pathIncidentFaces, pathIndices, pathBoundary);
  // and then we can use the previous method
  eik_g->mesh2 = mesh2;
  printf("Points 5: %f %f\n", eik_g->mesh2->points[5][0], eik_g->mesh2->points[5][1]);
  printf("Points 5: %f %f\n", mesh2->points[5][0], mesh2->points[5][1]);
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

void printInfoFanUpdate(fanUpdateS *fanUpdate) {
  int i;
  size_t nRegions = fanUpdate->triFan->nRegions;
  size_t nIndCrTop, nIndStTop, nPoints;
  nIndCrTop = fanUpdate->nIndCrTop;
  nIndStTop = fanUpdate->nIndStTop;
  nPoints = 2*nRegions + 1 + 2*nIndCrTop + 2*nIndStTop;
  printf("\n\nPRINTING INFORMATION TRIANGLE FAN UPDATE\n");
  printf("\nTHat: %fl\n", fanUpdate->THat);
  printf("\nGradHat:   %fl    %fl\n", fanUpdate->gradHat[0], fanUpdate->gradHat[1]);
  printf("Triangle fan:\n");
  printEverythingTriFan(fanUpdate->triFan);
  printf("\nParams: \n");
  for(i = 0; i<2*nRegions; i++){
    printf("  %fl  ", fanUpdate->params[i]);
  }
  printf("\nT0: \n");
  printf("%fl", fanUpdate->T0);
  printf("\ngrad0: \n");
  printf("%fl  %fl  ", fanUpdate->grad0[0], fanUpdate->grad0[1]);
  printf("\nT1: \n");
  printf("%fl", fanUpdate->T1);
  printf("\ngrad1: \n");
  printf("%fl  %fl  ", fanUpdate->grad1[0], fanUpdate->grad1[1]);

  printf("\nnIndCrTop:\n");
  printf("%zu", nIndCrTop);
  printf("\nindCrTop: \n");
  for(i = 0; i<nIndCrTop; i++){
    printf("  %zu  ", fanUpdate->indCrTop[i]);
  }
  printf("\nparamsCrTop: \n");
  for(i = 0; i<2*nIndCrTop; i++){
    printf("  %fl  ", fanUpdate->paramsCrTop[i]);
  }

  printf("\nnIndStTop:\n");
  printf("%zu", nIndStTop);
  printf("\nindStTop: \n");
  for(i = 0; i<nIndStTop; i++){
    printf("  %zu  ", fanUpdate->indStTop[i]);
  }
  printf("\nparamsStTop: \n");
  for(i = 0; i<2*nIndStTop; i++){
    printf("  %fl  ", fanUpdate->paramsStTop[i]);
  }

  printf("\nGrads:\n");
  for(i = 0; i<nPoints; i++){
    printf("  |%fl    %fl|  ", fanUpdate->grads[i][0], fanUpdate->grads[i][1]);
  }

  printf("\nPath:\n");
  for(i = 0; i<nPoints; i++){
    printf("  |%fl    %fl|  ", fanUpdate->path[i][0], fanUpdate->path[i][1]);
  }

  
}

void printAllInfoMesh(eik_gridS *eik_g){
  printEverythingInMesh(eik_g->mesh2);
}

void findEdgesOnValidFront(eik_gridS *eik_g, size_t index0, int indices1[2], int indices2[2], int firstTriangles[2]) {
  // given a newly accepted node with index index0 go around its neighbors
  // and find (at most 2) indices index1 such that the edge x0 x1 is on the valid fron
  // notice that indices1 is an int because it's going to have a -1 if we just have one valid edge from x0
  assert( eik_g->current_states[index0] == 2); // this HAS to be a valid node
  int nNeis, j, i;
  size_t thisNeighbor, previousNeighbor, thisTriangle, firstNeighbor;
  size_t possibleTriangles[2], possibleThirdVertices[2];
  j = 0;
  i = 1;
  indices1[0] = -1;
  indices1[1] = -1; // no neighbors found yet
  indices2[0] = -1;
  indices2[1] = -1;
  firstTriangles[0] = -1;
  firstTriangles[1] = -1; // no triangles found yet
  nNeis = eik_g->mesh2->neighbors[index0].len; // get the number of neighbors of index0
  // just start with any neighbor (doesn't really matter which one)
  firstNeighbor = (size_t)eik_g->mesh2->neighbors[index0].neis_i[0];
  thisNeighbor = firstNeighbor;
  previousNeighbor = index0;
  // go around
  printf("\n\n\nFinding valid edge from index0: %zu\n", index0);
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
    printf("PreviousNeighbor: %zu  with current state: %zu\n This Neighbor: %zu  with current state: %zu\n",
	   previousNeighbor, eik_g->current_states[previousNeighbor], thisNeighbor, eik_g->current_states[thisNeighbor]);
    // after we determined which is the neighbor that we need to cycle next
    if( eik_g->current_states[previousNeighbor] == 2 & eik_g->current_states[thisNeighbor] != 2 ){
      // thisNeighbor is part of the valid front
      indices1[j] = (int)previousNeighbor;
      indices2[j] = (int)thisNeighbor;
      firstTriangles[j] = (int)thisTriangle;
      j ++;
    }
    else if( eik_g->current_states[thisNeighbor] == 2 & eik_g->current_states[previousNeighbor] != 2){
      indices1[j] = (int)thisNeighbor;
      indices2[j] = (int)previousNeighbor;
      firstTriangles[j] = (int)thisTriangle;
      j ++;
    }
    i++;
  }
  if( eik_g->current_states[firstNeighbor] == 2 & eik_g->current_states[thisNeighbor] != 2){
    indices1[j] = (int)firstNeighbor;
    indices2[j] = (int)thisNeighbor;
    firstTriangles[j] = (int)faceBetween3Points(eik_g->mesh2, index0, thisNeighbor, firstNeighbor);
  }
  else if( eik_g->current_states[firstNeighbor] != 2 & eik_g->current_states[thisNeighbor] == 2){
    indices1[j] = (int)thisNeighbor;
    indices2[j] = (int)firstNeighbor;
    firstTriangles[j] = (int)faceBetween3Points(eik_g->mesh2, index0, thisNeighbor, firstNeighbor);
  }
  // what if the first and the last ones are the ones we need to consider?
  
  // figure out if we actually have at least one neighbor
  printf("Current eik status of index0: %zu\n", eik_g->current_states[index0]);
  printf("Indices1 found: indices1[0]: %d,   indices1[1]: %d\n", indices1[0], indices1[1]);
  printf("Indices2 found: indices2[0]: %d,   indices2[1]: %d\n", indices2[0], indices2[1]);
  printf("firstTriangles found: firstTriangles[0]: %d,   firstTriangles[1]: %d\n",
	 firstTriangles[0], firstTriangles[1]);
  assert( indices1[0] != -1 | indices1[1] != -1);
}




void initTriFan(eik_gridS *eik_g, triangleFanS *triFan,
		size_t index0, size_t index1, size_t index2,
		size_t indexHat, size_t firstTriangle, double *angleMax) {
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
  angleMax[0] = angleThreePoints(xk, x0, xk1); // first angle
  thisAngle = angleMax[0];
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
    if( angleRegion > angleMax[0] ){
      angleMax[0] = angleRegion;
    }
    printf("angleRegion %f, current angleMax: %f\n", angleRegion, angleMax[0]);
    if( angleMax[0] > pi ){
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
  int i = 3;
  indexk = index1;
  indexk1 = index2;
  thisTriangle = firstTriangle;
  listIndicesNodes[0] = index0;
  listIndicesNodes[1] =  index1;
  listIndicesNodes[2] =  index2;
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
  angleMax[0] = angleMax[0]/pi; // for easy comparison
  //////////////////////////
  // now we have the list of indices, we can initialize the triangle
  printf("\n Trying to initialize a triangleFan with  index0: %zu,  index1:  %zu,  indexHat:  %zu\n",
	 index0, index1, indexHat);
  printf("For this init we have %zu Regions\n", nRegions);
  triangleFan_initFromIndices(triFan, eik_g->mesh2, nRegions, index0,
			      index1, indexHat, listIndicesNodes);
}

void createJSONFile(fanUpdateS *fanUpdate, char const *path) {
  int i;
  FILE *fp = fopen(path, "w");

  fprintf(fp, "{");

  // write x0
  fprintf(fp, "\"x0\": [%g, %g], ",
	  fanUpdate->triFan->x0[0], fanUpdate->triFan->x0[1]);

  // write T0
  fprintf(fp, "\"T0\": %g, ", fanUpdate->T0 );

  // write grad0
  fprintf(fp, "\"grad0\": [%g, %g], ",
	  fanUpdate->grad0[0], fanUpdate->grad0[1]);


  // write x1
  fprintf(fp, "\"x1\": [%g, %g], ",
	  fanUpdate->triFan->x1[0], fanUpdate->triFan->x1[1]);

  // write T1
  fprintf(fp, "\"T1\": %g, ", fanUpdate->T1 );

  // write grad1
  fprintf(fp, "\"grad1\": [%g, %g], ",
	  fanUpdate->grad1[0], fanUpdate->grad1[1]);


  // write xHat
  fprintf(fp, "\"xHat\": [%g, %g], ",
	  fanUpdate->triFan->xHat[0], fanUpdate->triFan->xHat[1]);


  // use a for loop to add the list of indices
  fprintf(fp, "\"listIndices\": [");
  for(i = 0; i<2*fanUpdate->triFan->nRegions + 1; i++){
    fprintf(fp, "%g", fanUpdate->triFan->listIndices[i]);
    if( i < 2*fanUpdate->triFan->nRegions ){
      fprintf(fp, ",");
    }
  }
  fprintf(fp, "],");



  // use a loop to add the list of xk
  fprintf(fp, "\"listxk\": [");
  for(i = 0; i<fanUpdate->triFan->nRegions + 2; i++){
    fprintf(fp, "[%g, %g]", fanUpdate->triFan->listxk[i][0], fanUpdate->triFan->listxk[i][1]);
    if( i < fanUpdate->triFan->nRegions + 1){
      fprintf(fp, ",");
    }
  }
  fprintf(fp, "],");


  // use a loop to add the list of B0k
  fprintf(fp, "\"listB0k\": [");
  for(i = 0; i<fanUpdate->triFan->nRegions + 1; i++){
    fprintf(fp, "[%g, %g]", fanUpdate->triFan->listB0k[i][0], fanUpdate->triFan->listB0k[i][1]);
    if( i < fanUpdate->triFan->nRegions ){
      fprintf(fp, ",");
    }
  }
  fprintf(fp, "],");


  // use a loop to add the list of B0k
  fprintf(fp, "\"listBk\": [");
  for(i = 0; i<fanUpdate->triFan->nRegions + 1; i++){
    fprintf(fp, "[%g, %g]", fanUpdate->triFan->listBk[i][0], fanUpdate->triFan->listBk[i][1]);
    if( i < fanUpdate->triFan->nRegions ){
      fprintf(fp, ",");
    }
  }
  fprintf(fp, "],");


  // use a loop to add the list of BkBk1
  fprintf(fp, "\"listBkBk1\": [");
  for(i = 0; i<2*fanUpdate->triFan->nRegions; i++){
    fprintf(fp, "[%g, %g]", fanUpdate->triFan->listBkBk1[i][0], fanUpdate->triFan->listBkBk1[i][1]);
    if( i < 2*fanUpdate->triFan->nRegions-1 ){
      fprintf(fp, ",");
    }
  }
  fprintf(fp, "],");


  // OPTIONS FOR PLOTTING CHANGE THIS ACCORDINGLY
  fprintf(fp, "\"plotBefore\": 1, \"plotAfter\": 1, \"plotOpti\": 1");
  
  

  fprintf(fp, "}");

  fclose(fp);
}



void deserializeJSONoutput(fanUpdateS *fanUpdate, json_object *output_obj) {
  //given a json_object output from the Python optimizer we deserialize it
  // and save it to fanUpdate
  size_t nRegions = fanUpdate->triFan->nRegions; // useful to have
  // malloc what we can malloc
  fanUpdate->params = malloc((2*nRegions + 1)*sizeof(double));
  // START READING
  assert(output_obj != NULL); // we need a json object, not a null
  // the things we need
  size_t nIndCrTop, nIndStTop;
  double THat;

  // get nIndCrTop, nIndStTop, THat
  nIndCrTop = (size_t)json_object_get_double(json_object_object_get(output_obj, "nIndCrTop"));
  //printf("\nnIndCrTop %zu\n", nIndCrTop);
  nIndStTop = (size_t)json_object_get_double(json_object_object_get(output_obj, "nIndStTop"));
  //printf("\nnIndStTop %zu\n", nIndStTop);
  THat = json_object_get_double(json_object_object_get(output_obj, "THat"));
  //printf("\nTHat %g\n", THat);

  // now for the lists
  json_object *output_list, *value_arr;
  int i, k;

  // get CrTop - extract lists of size_t from output object
  size_t *indCrTop;
  double *paramsCrTop;
  if( nIndCrTop == 0 ){
    indCrTop = NULL;
    paramsCrTop = NULL;
  }
  else{
    size_t indCrTop_vec[nIndCrTop];
    double paramsCrTop_vec[2*nIndCrTop];
    // start with indCrTop
    output_list = json_object_object_get(output_obj, "indCrTop");
    if(output_list == NULL || !json_object_is_type(output_list,  json_type_array)) {
      printf("\nProblem when opening the json file for indCrTop\n");
      exit(EXIT_FAILURE);
    }
    for(i = 0; i<nIndCrTop; i++){
      value_arr = json_object_array_get_idx(output_list, i); // get the value stored on output_list at index i
      if (value_arr == NULL || !json_object_is_type(value_arr, json_type_double)) {
	printf("Error extracting double value from indCrTop\n");
	exit(EXIT_FAILURE);
      }
      indCrTop_vec[i] = (size_t) json_object_get_double(value_arr);
    }
    // now with paramsCrTop
    output_list = json_object_object_get(output_obj, "paramsCrTop");
    if(output_list == NULL || !json_object_is_type(output_list,  json_type_array)) {
      printf("\nProblem when opening the json file for paramsCrTop\n");
      exit(EXIT_FAILURE);
    }
    for(i = 0; i<2*nIndCrTop; i++){
      value_arr = json_object_array_get_idx(output_list, i); // get the value stored on output_list at index i
      if (value_arr == NULL || !json_object_is_type(value_arr, json_type_double)) {
	printf("Error extracting double value from paramsCrTop\n");
	exit(EXIT_FAILURE);
      }
      paramsCrTop_vec[i] = json_object_get_double(value_arr);
    }
    // pointer to the right places
    indCrTop = malloc(nIndCrTop*sizeof(size_t));
    paramsCrTop = malloc(2*nIndCrTop*sizeof(double));
    indCrTop = &indCrTop_vec[0];
    paramsCrTop = &paramsCrTop_vec[0];
  }


  // get StTop - extract lists of size_t from output object
  size_t *indStTop;
  double *paramsStTop;
  if( nIndStTop == 0 ){
    indStTop = NULL;
    paramsStTop = NULL;
  }
  else{
    size_t indStTop_vec[nIndStTop];
    double paramsStTop_vec[2*nIndStTop];
    // start with indStTop
    output_list = json_object_object_get(output_obj, "indStTop");
    if(output_list == NULL || !json_object_is_type(output_list,  json_type_array)) {
      printf("\nProblem when opening the json file for indStTop\n");
      exit(EXIT_FAILURE);
    }
    for(i = 0; i<nIndStTop; i++){
      value_arr = json_object_array_get_idx(output_list, i); // get the value stored on output_list at index i
      if (value_arr == NULL || !json_object_is_type(value_arr, json_type_double)) {
	printf("\nError extracting double value from indStTop\n");
	exit(EXIT_FAILURE);
      }
      indStTop_vec[i] = (size_t) json_object_get_double(value_arr);
    }
    // now with paramsStTop
    output_list = json_object_object_get(output_obj, "paramsStTop");
    if(output_list == NULL || !json_object_is_type(output_list,  json_type_array)) {
      printf("\nProblem when opening the json file for paramsStTop\n");
      exit(EXIT_FAILURE);
    }
    for(i = 0; i<2*nIndStTop; i++){
      value_arr = json_object_array_get_idx(output_list, i); // get the value stored on output_list at index i
      if (value_arr == NULL || !json_object_is_type(value_arr, json_type_double)) {
	printf("\nError extracting double value from paramsStTop\n");
	exit(EXIT_FAILURE);
      }
      paramsStTop_vec[i] = json_object_get_double(value_arr);
    }
    // pointer to the right places
    indStTop = malloc(nIndStTop*sizeof(size_t));
    paramsStTop = malloc(2*nIndStTop*sizeof(double));
    indStTop = &indStTop_vec[0];
    paramsStTop = &paramsStTop_vec[0];
  }


  // get grads
  double (*grads)[2];
  size_t nPath = 2*nRegions + 1 + 2*nIndCrTop + 2*nIndStTop;
  grads = malloc(2*nPath*sizeof(double));
  output_list = json_object_object_get(output_obj, "grads");
  if(output_list == NULL || !json_object_is_type(output_list, json_type_array)){
    printf("\nProblem when opening the json file and reading grads\n");
    exit(EXIT_FAILURE);
  }
  for(i = 0; i<2*nPath; i++){
    value_arr = json_object_array_get_idx(output_list, i);
    if(value_arr == NULL || !json_object_is_type(value_arr, json_type_double)) {
      printf("\nError extracting double value from grads\n");
      exit(EXIT_FAILURE);
    }
    if( (i % 2 == 0) ){
      grads[i][0] = json_object_get_double(value_arr);
    }
    else{
      grads[i][1] = json_object_get_double(value_arr);
    }
  }


  // get path
  double (*path)[2];
  path = malloc(2*nPath*sizeof(double));
  output_list = json_object_object_get(output_obj, "path");
  if(output_list == NULL || !json_object_is_type(output_list, json_type_array)){
    printf("\nProblem when opening the json file and reading path\n");
    exit(EXIT_FAILURE);
  }
  for(i = 0; i<2*nPath; i++){
    value_arr = json_object_array_get_idx(output_list, i);
    if(value_arr == NULL || !json_object_is_type(value_arr, json_type_double)) {
      printf("\nError extracting double value from path\n");
      exit(EXIT_FAILURE);
    }
    if( (i % 2 == 0) ){
      path[i][0] = json_object_get_double(value_arr);
    }
    else{
      path[i][1] = json_object_get_double(value_arr);
    }
  }


  // get gradHat
  double gradHat[2];
  output_list = json_object_object_get(output_obj, "gradHat");
  if(output_list == NULL || !json_object_is_type(output_list, json_type_array)){
    printf("\nProblem when opening the json file and reading gradHat\n");
    exit(EXIT_FAILURE);
  }
  for(i = 0; i<2; i++){
    value_arr = json_object_array_get_idx(output_list, i);
    if(value_arr == NULL || !json_object_is_type(value_arr, json_type_double)) {
      printf("\nError extracting double value from gradHat\n");
      exit(EXIT_FAILURE);
    }
    gradHat[i] = json_object_get_double(value_arr);
  }


  // save everything to fanUpdateS
  fanUpdate->nIndCrTop = nIndCrTop;
  fanUpdate->indCrTop = indCrTop;
  fanUpdate->paramsCrTop = paramsCrTop;
  fanUpdate->nIndStTop = nIndStTop;
  fanUpdate->indStTop = indStTop;
  fanUpdate->paramsStTop = paramsStTop;
  fanUpdate->THat = THat;
  fanUpdate->grads = grads;
  fanUpdate->path = path;
  fanUpdate->gradHat[0] = gradHat[0];
  fanUpdate->gradHat[1] = gradHat[1];

  // print the fan update just to see what happened here
  printInfoFanUpdate(fanUpdate);

}



void optimizeTriangleFan_wPython(fanUpdateS *fanUpdate) {
  // using pipes and a lot of fancy methos call python from here
  // and optimize the triangle fan in fanUpdate, save all info
  int pipefd[2]; // pipe
  int pid; // proces id


  // create JSON object to represent input data
  char* input_json; //
  createJSONFile(fanUpdate, "/Users/marianamartinez/Documents/Curvy-JMM/JMM/update.json");
  //createJSONinput(fanUpdate, &input_json); // put everything in the json order we want
  //printf("created json file for input\n");
  /* if( pipe(pipefd) == -1){ */
  /*   printf("\nProblem opening pipe\n"); */
  /*   exit(EXIT_FAILURE); */
  /* } */
  /* pid = fork(); */
  /* if(pid == 0){ */
  /*   // we are in the child process */
  /*   execlp("python", "python", "stepWithPython.py", "update.json", NULL); */
  /*   printf("executed python\n"); */
  /* } */
  /* else{ */
  /*   wait(NULL); */
  /*   char buffer[5000]; */
  /*   int fd = open("/Users/marianamartinez/Documents/Curvy-JMM/JMM/updateSolve.json", O_RDONLY); */
  /*   if( fd == -1){ */
  /*     printf("\nProblem when opening the output json from Python\n"); */
  /*     exit(EXIT_FAILURE); */
  /*   } */
  /*   ssize_t num_read = read(fd, buffer, sizeof(buffer)); */
  /*   json_object *output_obj = json_tokener_parse(buffer); */
  /*   deserializeJSONoutput(fanUpdate, output_obj); */
  /* } */
  
  //execlp("python", "python", "stepWithPython.py", NULL);
  FILE *fp = popen("python3 /Users/marianamartinez/Documents/Curvy-JMM/JMM/stepWithPython.py /Users/marianamartinez/Documents/Curvy-JMM/JMM/update.json", "r");
  printf("executed python\n");
  assert(fp != NULL);
  char buffer[5000];
  // print what it gets
  char buf[5000 + 1];
  printf("reading from pipe:\n");
  while (fgets(buf, 5000, fp) != NULL){
    printf("%s", buf);
  }
  // try to separate this
  double output[3];
  separateARowDb(&buf[0], 3, output);
  printf("Row in double: %f   %f   %f\n", output[0], output[1], output[2]);
  // add this information to fanUpdate
  fanUpdate->THat = output[0];
  fanUpdate->gradHat[0] = output[1];
  fanUpdate->gradHat[1] = output[2];
  /* json_object *output_obj = json_tokener_parse(buf); */
  /* deserializeJSONoutput(fanUpdate, output_obj); */
  /* pclose(fp); */
  /* ssize_t num_read = read(fd, buffer, sizeof(buffer)); */
  /* json_object *output_obj = json_tokener_parse(buffer); */
  /* deserializeJSONoutput(fanUpdate, output_obj); */

}


void updateOneWay(eik_gridS *eik_g, size_t index0, size_t index1, size_t index2,
		  int indexStop_int, size_t firstTriangle) {
  printf("\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n");
  // given a starting edge on the valid front update a far/trial
  // set of triangles until indexStop is reached (going on the direction of firstTriangle)
  size_t possibleNextTriangles[2], possibleNextIndices[2];
  int i;
  double T0, grad0[2], T1, grad1[2], *angleMax;
  angleMax = malloc(sizeof(double));
  T0 = eik_g->eik_vals[index0];
  grad0[0] = eik_g->grads[index0][0];
  grad0[1] = eik_g->grads[index0][1];
  T1 = eik_g->eik_vals[index1];
  grad1[0] = eik_g->grads[index1][0];
  grad1[1] = eik_g->grads[index1][1];
  printf("Using index0: %zu  with T0: %f,  grad0: %f  %f,  index1: %zu  with T1: %f,  grad1: %f  %f",
	 index0, T0, grad0[0], grad0[1], index1, T1, grad1[0], grad1[1]);
  size_t indexStop;
  if(indexStop < 0){
    // meaning that we dont have an index stop, this is an artificial value
    indexStop = index0;
  }
  else{
    indexStop = (size_t)indexStop_int;
  }
  // for the first triangle
  triangleFanS *currentTriangleFan;
  triangleFan_alloc(&currentTriangleFan);
  fanUpdateS *currentTriangleFanUpdate;
  fanUpdate_alloc(&currentTriangleFanUpdate);
  size_t prevTriang, thisTriang;
  prevTriang = firstTriangle;
  thisTriang = firstTriangle;
  size_t *allPossibleIndicesNodes, nRegions; // in here we are going to store all possible
  // nodes for the triangle fans, the key is to use nRegions to initialize correctly
  // a specific triangle fan
  allPossibleIndicesNodes = malloc(eik_g->mesh2->neighbors->len*sizeof(size_t));
  ///////////////////////////////////////////////////////////
  // FOR THE FIRST ONE
  size_t indexHat, prevUpdated;
  indexHat = index2;
  prevUpdated = index2;
  allPossibleIndicesNodes[0] = index0;
  allPossibleIndicesNodes[1] = index1;
  allPossibleIndicesNodes[2] = index2;
  i = 3;
  nRegions = 1; // because we have one triangle to begin with
  while(  indexHat != indexStop & eik_g->current_states[indexHat] != 2 &
	  i<eik_g->mesh2->neighbors->len ){
    // while we are not on the other valid edge and while the current node we are trying
    // to update is not valid (yet)
    // initialize currentTriangleFan
    printf("\n   Trying to update neighbor: %zu\n", indexHat);
    printf("\n    Index0: %zu     Index1: %zu      Index2: %zu\n", index0, index1, index2);
    initTriFan(eik_g, currentTriangleFan, index0, index1, index2,
	       indexHat, firstTriangle, angleMax);
    printf("Angle in this triangleFan %f\n", angleMax[0]);
    if( fabs(angleMax[0]) > 1){
      printf("\n\nangleMax is greater than pi\n");
      return;
    }
    // initialize the triangle update
    fanUpdate_initPreOpti(currentTriangleFanUpdate, currentTriangleFan,
			  T0, grad0, T1, grad1);
    currentTriangleFanUpdate->THat = eik_g->eik_vals[indexHat];
    //////// OPTIMIZE!
    optimizeTriangleFan_wPython(currentTriangleFanUpdate);
    printf("THat found: %f\n", currentTriangleFanUpdate->THat);
    // if we found a good update, update the values
    if( eik_g->eik_vals[indexHat] > currentTriangleFanUpdate->THat ){
      assert( eik_g->current_states[indexHat] != 2); // no update on valid nodes
      // meaning we found a better update
      eik_g->eik_vals[indexHat] = currentTriangleFanUpdate->THat;
      eik_g->grads[indexHat][0] = currentTriangleFanUpdate->gradHat[0];
      eik_g->grads[indexHat][1] = currentTriangleFanUpdate->gradHat[1];
      if( eik_g->current_states[indexHat] == 1 ){
	// current state of xHat = trial
	update(eik_g->p_queueG, currentTriangleFanUpdate->THat, indexHat); // update THat in the queue
      }
      else{
	// current state of xHat = far
	eik_g->current_states[indexHat] = 1; // set current state of xHat to trial
	insert(eik_g->p_queueG, currentTriangleFanUpdate->THat, indexHat); // insert THat to the queue
      }
    }
    //////// MARCH TO THE OTHER TRIANGLE
    prevTriang = thisTriang;
    twoTrianglesFromEdge(eik_g->mesh2, index0, indexHat, possibleNextTriangles, possibleNextIndices);
    if( possibleNextTriangles[0] == prevTriang ){
      // we should march in the possibleNextTriangles[1] direction
      thisTriang = possibleNextTriangles[1];
      indexHat = possibleNextIndices[1];
    }
    else{
      thisTriang = possibleNextTriangles[0];
      indexHat = possibleNextIndices[0];
    }
    allPossibleIndicesNodes[i] = indexHat;
    nRegions ++;
    i ++;
  }
  triangleFan_dalloc(&currentTriangleFan);
  fanUpdate_dalloc(&currentTriangleFanUpdate);
}



void addNeighbors_fromAccepted(eik_gridS *eik_g, size_t minIndex) {
  // given a recently accepted node with index minIndex we update its neighbors
  // using triangle fans and up to two valid edges
  // update from edges on valid front
  int indices1[2], indices2[2], firstTriangles[2];
  findEdgesOnValidFront(eik_g, minIndex, indices1, indices2, firstTriangles);
  //printf("Edges from valid, indices1: %d %d\n indices2: %d %d\n firsttriangles: %d %d",
  //	 indices1[0], indices1[1], indices2[0], indices2[1], firstTriangles[0], firstTriangles[1]);
  // start with indices2[0]
  size_t index1, index2, firstTriangle;
  int indexStop_int;
  index1 = (size_t) indices1[0];
  index2 = (size_t) indices2[0];
  firstTriangle = (size_t) firstTriangles[0];
  indexStop_int = indices1[1];
  updateOneWay(eik_g, minIndex, index1, index2, indexStop_int, firstTriangle);
  // now with indices2[1] if doable
  if( indices2[1] > 0){
    index1 = (size_t) indices1[1];
    index2 = (size_t) indices2[1];
    firstTriangle = (size_t) firstTriangles[1];
    indexStop_int = indices1[0];
    updateOneWay(eik_g, minIndex, index1, index2, indexStop_int, firstTriangle);
  }
}


/* void triangleFanUpdate_pointNear(eik_gridS *eik_g, size_t index0, size_t indexStart, double eta) { */
/*   // update the triangleFanUpdate for those point who were intialized at the */
/*   // begining of the method, at initializePoitnsNear */
/*   double x0MinxStart[2], xStart[2], x0[2]; */
/*   xStart[0] = eik_g->mesh2->points[indexStart][0]; */
/*   xStart[1] = eik_g->mesh2->points[indexStart][1]; */
/*   x0[0] = eik_g->mesh2->points[index0][0]; */
/*   x0[1] = eik_g->mesh2->points[index0][1]; */
/*   vec2_subtraction(x0, xStart, x0MinxStart); // for both T0 and grad0 */
/*   double grad0[2], T0; */
/*   T0 = eta*l2norm(x0MinxStart); */
/*   printf("T0: %f\n", T0); */
/*   grad0[0] = eta*x0MinxStart[0]; */
/*   grad0[1] = eta*x0MinxStart[1]; */
/*   // initialize the triangle Fan for the triangle fan update */
/*   triangleFanS *triFan; */
/*   triangleFan_alloc(&triFan); */
/*   triangleFan_init(triFan, 0, xStart, xStart, x0, NULL, NULL, NULL, NULL, NULL, NULL, NULL); */
/*   fanUpdate_init(&eik_g->fanUpdate[index0], triFan, NULL, */
/* 		 0, NULL, 0, NULL, */
/* 		 0, NULL, NULL, */
/* 		 0, NULL, NULL, */
/* 		 T0, NULL, NULL, grad0); */
/*   triangleFan_dalloc(&triFan); */
/* } */



void initializePointsNear(eik_gridS *eik_g, double rBall) {
  // THIS JUST SETS THE CURRENT STATE TO VALID AND ADDS THE TRUE EIKONAL
  // given a ball of radius rBall around the initial points we initialize all the points inside those balls with the
  // true value of the eikonal (i.e. the distance times the index of refraction). We are going to assume that
  // all those balls belong to the same regions of indices of refraction
  double initialEta, xCurrent[2], xStart[2], xMinxStart[2], normCurrent;
  size_t indexStart, nei;
  int j, i;
  // initialize the initial index of refraction
  printf("\n\n\nStarting to initialize points near\n\n\n\n");
  for(j = 0; j<eik_g->nStart; j++){
    deleteRoot(eik_g->p_queueG);
    indexStart = eik_g->start[j];
    nei = eik_g->mesh2->neighbors[indexStart].neis_i[0];
    //printf("\nNumber of neighbors of indexStart: %zu", nei);
    initialEta = minEtaFromTwoPoints(eik_g->mesh2, indexStart, nei);
    //printf("\nStarting with initial point: %zu", indexStart);
    xStart[0] = eik_g->mesh2->points[ indexStart ][0];
    xStart[1] = eik_g->mesh2->points[ indexStart ][1];
    //printf("\nCoordinates of xStart: %f   %f\n", xStart[0], xStart[1]);
    for(i = 0; i<eik_g->mesh2->nPoints; i++){
      //printf("\nConsidering the point: %d\n", i);
      xCurrent[0] = eik_g->mesh2->points[i][0];
      xCurrent[1] = eik_g->mesh2->points[i][1];
      //printf("Current points coordinates: %f  %f\n", xCurrent[0], xCurrent[1]);
      vec2_subtraction( xCurrent, xStart, xMinxStart );
      //printf("xMinxStart: %f  %f \n", xMinxStart[0], xMinxStart[1]);
      normCurrent = l2norm(xMinxStart);
      //printf("%f\n", normCurrent);
      //printf("\nAfter assignint normCurrent\n");
      if(normCurrent < rBall ){
	printf("Point inside ball\n");
        if( eik_g->current_states[i] == 1 ){
          // if it was previously considered as trial we need to delete this from the queue directly
	  printf("Point was previously considered trial\n");
          delete_findIndex(eik_g->p_queueG, i);
        }
	printf("Initialize the current state as valid, eikVal and so on\n");
        // if this happens, this point is "close enough" to the starting point so that we can initialize it
        eik_g->current_states[i] = 2; // we initialized it directly
        eik_g->eik_vals[i] = initialEta*normCurrent; // add their true Eikonal value
	printf("Update point near\n");
	eik_g->grads[i][0] = initialEta*xMinxStart[0];
	eik_g->grads[i][1] = initialEta*xMinxStart[1];
      }
    }
  }
  printGeneralInfo(eik_g);
  // now we initialize the trial nodes
  size_t nNeis, flag;
  for( i = 0; i<eik_g->mesh2->nPoints; i++) {
    if(eik_g->current_states[i] == 2){
      // pick a point that is set to valid
      flag = 0;
      nNeis = eik_g->mesh2->neighbors[i].len;
      for(j = 0; j<nNeis; j++){
	nei = eik_g->mesh2->neighbors[i].neis_i[j];
	if(eik_g->current_states[nei] != 2){
	  flag = 1;
	}
      }
      if( flag == 1 ){
	// this means that this node has at least one neighbor that is not valid
	printf("Adding the neighbors from: %d\n", i);
	addNeighbors_fromAccepted(eik_g, i);
      }
    }
  }

}

void popAddNeighbors(eik_gridS *eik_g) {
  // int nNeighs;
  int minIndex = indexRoot(eik_g->p_queueG);
  deleteRoot(eik_g->p_queueG); // delete the root from the priority queue
  eik_g->current_states[minIndex] = 2; // set the newly accepted index to valid
  addNeighbors_fromAccepted(eik_g, minIndex); // add neighbors from the recently accepted index
}

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
      fwrite(eik_g->grads[i], sizeof(double), 2, fp);
    }

    fclose(fp);
}




