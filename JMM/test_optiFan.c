// testing the optimization method for the triangle fan point of view

#include "optiFan.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

int main(){
  optiFanS *optiFan;
  int nRegions;
  double x0[2], T0, x1[2], T1, xHat[2], (*points_fan)[2], (*B_x0)[2], (*B_xk)[2], *indicesRef;

  T0 = 1;
  T1 = 1.4;
  nRegions = 3;

  points_fan = malloc(4*2*sizeof(double));
  B_x0 = malloc(2*4*sizeof(double));
  B_xk = malloc(4*2*sizeof(double));
  indicesRef = malloc(5*sizeof(double));

  x0[0] = 0;
  x0[1] = 0;
  x1[0] = 2;
  x1[1] = -0.2;
  points_fan[0][0] = 0.0;
  points_fan[0][1] = 0.0;
  points_fan[1][0] = 2.0;
  points_fan[1][1] = -0.2;
  points_fan[2][0] = 1.5;
  points_fan[2][1] = 0.8;
  points_fan[3][0] = 0.2;
  points_fan[3][1] = 1.2;
  points_fan[4][0] = -0.8;
  points_fan[4][1] = 0.2;
  xHat[0] = -0.8;
  xHat[1] = 0.2;

  B_x0[0][0] = 0.91036648;
  B_x0[0][1] = 0.41380294;
  B_x0[1][0] = 0.5547002;
  B_x0[1][1] = 0.83205029;
  B_x0[2][0] = 0.09950372;
  B_x0[2][1] = 0.99503719;
  B_x0[3][0] = -0.92847669;
  B_x0[3][1] = 0.37139068;

  B_xk[0][0] = 0.85749293;
  B_xk[0][1] = -0.51449576;
  B_xk[1][0] = 0.99503719;
  B_xk[1][1] = -0.09950372;
  B_xk[2][0] = 0.447213;
  B_xk[2][1] = 0.89442719;
  B_xk[3][0] = -0.9284766;
  B_xk[3][1] = -0.37139068;

  indicesRef[0] = 0.8;
  indicesRef[1] = 1;
  indicesRef[2] = 1.001;
  indicesRef[3] = 1.005;
  indicesRef[4] = 1.5;

  optiFan_alloc(&optiFan);
  optiFan_init(optiFan, nRegions, x0, T0, x1, T1, xHat, indicesRef, points_fan, B_x0, B_xk);

  
}

