/* DIFFERENT TYPES OF UPDATES FOR SEMI LAGRANGIAN MARCHER IN 2D

These are the different types of updates fro the semi lagrangian
marcher given an unstructured triangle mesh in 2D which
also includes information about the gradient of points that
are on the boundary of the geometry.

*/

#include "updates_2D.h"
#include "linAlg.h"
#include "opti_method.h"
#include "SoSFunction.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

void info_update_alloc(info_updateS **info_update) {
  *info_update = malloc(sizeof(info_update));
  assert(*info_update != NULL);
}

void info_update_dealloc(info_updateS **info_update) {
  free(*info_update);
  *info_update = NULL;
}

void info_update_init(info_updateS *info_update, int indexAccepted, int x1_ind, int x2_ind, int xHat_ind, double T0, double grad0[2], double T1, double grad1[2], double indexRef_01, double indexRef_02) {
  // constructor for all the possible information needed for any update
  info_update->indexAccepted = indexAccepted;
  info_update->x1_ind = x1_ind;
  info_update->x2_ind = x2_ind;
  info_update->xHat_ind = xHat_ind;
  info_update->T0 = T0;
  info_update->grad0[0] = grad0[0];
  info_update->grad0[1] = grad0[1];
  info_update->T1 = T1;
  info_update->grad1[0] = grad1[0];
  info_update->grad1[1] = grad1[1];
  info_update->THat = -1; // initial value when starting, eikonal can't be negative useful for asserts
  info_update->indexRef_01 = indexRef_01;
  info_update->indexRef_02 = indexRef_02;
  info_update->lambda = 0; // initial value when starting
  info_update->mu = 0; // initial value when starting
}

void info_update_initCr(info_updateS *info_update, int indexAccepted, int xHat_ind, double T0, double indexMin){
  // constructor for the necessary information for a creeping ray update
  double index, grad0[2], grad1[2];
  grad0[0] = 0;
  grad0[1] = 0;
  grad1[0] = 0;
  grad1[1] = 0;
  info_update_init(info_update, indexAccepted, -1, -1, xHat_ind, T0, grad0, -1, grad1, indexMin, indexMin);
}

void info_update_initTwo(info_updateS *info_update, int indexAccepted, int x1_ind, int xHat_ind, double T0, double grad0[2], double T1, double grad1[2], double indexRef){
  // constructor for the necessary information for a two point update in free space
  info_update_init(info_update, indexAccepted, x1_ind, -1, xHat_ind, T0, grad0, T1, grad1, indexRef, indexRef);
}

void print_info_update(info_updateS *info_update){
  printf("Printing information about the recent update\n");
  printf("Index accepted, x0: %d\n", info_update->indexAccepted);
  printf("Index x1: %d\n", info_update->x1_ind);
  printf("Index x2: %d\n", info_update->x2_ind);
  printf("Index xHat: %d\n", info_update->xHat_ind);
  printf("T0: %lf\n", info_update->T0);
  printf("Grad T0:  %lf   |   %lf \n", info_update->grad0[0], info_update->grad0[1]);
  printf("T1: %lf\n", info_update->T1);
  printf("Grad T1:  %lf   |   %lf \n", info_update->grad1[0], info_update->grad1[1]);
  printf("THat: %lf\n", info_update->THat);
  printf("Index 01 %lf\n", info_update->indexRef_01);
  printf("Index 02 %lf\n", info_update->indexRef_02);
  printf("Lambda %lf\n", info_update->lambda);
  printf("Mu %lf\n", info_update->mu);
}

// First type of update: "creeping ray" type of update

void creepingUpdate(triMesh_2Ds *triM_2D, info_updateS *info_update) {
  // this is when x0, x1, and xHat are all on the boundary
  // second asserts, all of these points are on the boundary
  assert(triM_2D->boundary_tan[info_update->indexAccepted][0] != 0 || triM_2D->boundary_tan[info_update->indexAccepted][1] != 0 );
  assert(triM_2D->boundary_tan[info_update->xHat_ind][0] != 0 || triM_2D->boundary_tan[info_update->xHat_ind][1] != 0 );
  // if all of these are true then we can proceed with a creeping update
  // THat in this case is just the smallest index of refraction times
  // the arc length of the boundary's Hermite interpolation
  double x0[2], xHat[2], B0[2], BHat[2];
  x0[0] = triM_2D->points->x[info_update->indexAccepted];
  x0[1] = triM_2D->points->y[info_update->indexAccepted];
  xHat[0] = triM_2D->points->x[info_update->xHat_ind];
  xHat[1] = triM_2D->points->y[info_update->xHat_ind];
  B0[0] = triM_2D->boundary_tan[info_update->indexAccepted][0];
  B0[1] = triM_2D->boundary_tan[info_update->indexAccepted][1];
  printf("B0 : %lf | %lf \n", B0[0], B0[1]);
  printf("x0: %lf | %lf \n", x0[0], x0[1]);
  BHat[0] = triM_2D->boundary_tan[info_update->xHat_ind][0];
  BHat[1] = triM_2D->boundary_tan[info_update->xHat_ind][1];
  printf("BHat : %lf | %lf \n", BHat[0], BHat[1]);
  printf("xHat: %lf  | %lf \n", xHat[0], xHat[1]);
  // For Simpsons rule que need norm(B0) + norm(g(0.5)) + norm(BHat)
  // for g(0.5):
  double constx0[2], constB0[2], constxHat[2], constBHat[2], auxSum1[2], auxSum2[2], auxSum3[2], gHalves;
  // multiply times the necessary scalars for g(0.5)
  scalar_times_2vec(-1.5, x0, constx0);
  scalar_times_2vec(-0.25, B0, constB0);
  scalar_times_2vec(1.5, xHat, constxHat);
  scalar_times_2vec(-0.25, BHat, constBHat);
  vec2_addition(constx0, constB0, auxSum1);
  vec2_addition(constxHat, constBHat, auxSum2);
  vec2_addition(auxSum1, auxSum2, auxSum3);
  gHalves = l2norm(auxSum3); // norm(g(0.5))
  double normB0, normBHat;
  normB0 = l2norm(B0);
  normBHat = l2norm(BHat);
  // finally we can compute THat using Simpson's rule
  info_update->THat = info_update->T0 + info_update->indexRef_01/6*(normB0 + 4*gHalves + normBHat);
}

void simple_TwoPointUpdate(triMesh_2Ds *triM_2D, info_updateS *info_update){
  // this is a simple two point update, meaning that there are no changes in the region
  // being considered (i.e. no change in the index of refraction). Then we can linearly
  // interpolate x0 to x1, parametrizing it with lambda. T(xLam) is approximated using
  // a cubic hermite polynomial

  // first we need to find the optimal lambda to define xLam
  double tol, lambda0, lambda1, T0, T1, grad0[2], grad1[2], x0[2], x1[2], xHat[2], indexRef;
  int maxIter;
  tol = 0.0000000001;
  maxIter = 50;
  lambda0 = 0;
  lambda1 = 1;
  T0 = info_update->T0;
  T1 = info_update->T1;
  grad0[0] = info_update->grad0[0];
  grad0[1] = info_update->grad0[1];
  grad1[0] = info_update->grad1[0];
  grad1[1] = info_update->grad1[1];
  x0[0] = triM_2D->points->x[info_update->indexAccepted];
  x0[1] = triM_2D->points->y[info_update->indexAccepted];
  x1[0] = triM_2D->points->x[info_update->x1_ind];
  x1[1] = triM_2D->points->y[info_update->x1_ind];
  xHat[0] = triM_2D->points->x[info_update->xHat_ind];
  xHat[1] = triM_2D->points->y[info_update->xHat_ind];
  indexRef = info_update->indexRef_01;

  info_update->lambda = secant_freeSpace(lambda0, lambda1, T0, T1, grad0, grad1, x0, x1, xHat, tol, maxIter, indexRef); // optimal lambda found
  info_update->THat = eikApprox_freeSpace(T0, T1, grad0, grad1, info_update->lambda, x0, x1, xHat, indexRef);
}

void fromBoundary_TwoPointUpdate(triMesh_2Ds *triM_2D, info_updateS *info_update) {
  // two point update but the segment x0x1 is on the boundary
  double tol, T0, grad0[2], B0[2], T1, grad1[2], B1[2], x0[2], x1[2], xHat[2], indexRef;
  int maxIter;
  tol = 0.0000000001;
  maxIter = 50;
  T0 = info_update->T0;
  T1 = info_update->T1;
  grad0[0] = info_update->grad0[0];
  grad0[1] = info_update->grad0[1];
  B0[0] = triM_2D->boundary_tan[info_update->indexAccepted][0];
  B0[1] = triM_2D->boundary_tan[info_update->indexAccepted][1];
  grad1[0] = info_update->grad1[0];
  grad1[1] = info_update->grad1[1];
  B1[0] = triM_2D->boundary_tan[info_update->x1_ind][0];
  B1[1] = triM_2D->boundary_tan[info_update->x1_ind][1];
  x0[0] = triM_2D->points->x[info_update->indexAccepted];
  x0[1] = triM_2D->points->y[info_update->indexAccepted];
  x1[0] = triM_2D->points->x[info_update->x1_ind];
  x1[1] = triM_2D->points->y[info_update->x1_ind];
  xHat[0] = triM_2D->points->x[info_update->xHat_ind];
  xHat[1] = triM_2D->points->y[info_update->xHat_ind];
  indexRef = info_update->indexRef_01;
  // Calculate the optimal lambda
  info_update->lambda = projectedGradient_fromEdge(0, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, tol, maxIter, indexRef);
  // Then calculate the optimal THat
  info_update->THat = fobjective_fromEdge(info_update->lambda, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
}






