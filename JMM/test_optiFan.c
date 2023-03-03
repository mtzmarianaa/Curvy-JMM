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


  printf("\n\n\n\nTest the projections back on to the feasible set\n");

  /* double lambda2min, B0[2], ykPrime[2], x2[2], B1[2], Bk_mu[2], B2[2], xlam[2]; */
  /* B0[0] = B_x0[1][0]; */
  /* B0[1] = B_x0[1][1]; */
  /* B1[0] = B_xk[0][0]; */
  /* B1[1] = B_xk[0][1]; */
  /* x2[0] = points_fan[2][0]; */
  /* x2[1] = points_fan[2][1]; */
  /* B2[0] = B_xk[1][0]; */
  /* B2[1] = B_xk[1][1]; */
  /* hermite_interpolationSpatial(0.85, x0, x1, B0, B1, ykPrime); */
  /* grad_hermite_interpolationSpatial(0.85, x0, x1, B0, B1, Bk_mu); */
  /* printf("ykPrime: %lf %lf\n", ykPrime[0], ykPrime[1]); */
  /* printf("Bkmu: %lf %lf\n", Bk_mu[0], Bk_mu[1]); */
  /* lambda2min = lambda_fromt1(0.3, x0, B0, ykPrime, Bk_mu, x2, B2, 0.000001, 100); */

  /* printf("\n\nThe minimum possible value for lambda2 is: %lf\n\n\n", lambda2min); */

  /* printf("T1 with this lamdba2: %lf\n", t1_ofLam(lambda2min, x0, B0, ykPrime, Bk_mu, x2, B2)); */

  /* hermite_interpolationSpatial(lambda2min, x0, x2, B0, B2, xlam); */
  /* printf("The interpolated value with lambdamin:   %lf  %lf\n", xlam[0], xlam[1]); */

  /* printf("\n\n\n\nTesting finding lamda4min for a type 4 triangle\n\n"); */

  /* double y3Prime[2], B03[2], B3[2], B0Hat[2], BHat[2], x3[2], B3_mu[2]; */
  /* B03[0] = B_x0[2][0]; */
  /* B03[1] = B_x0[2][1]; */
  /* B3[0] = B_xk[2][0]; */
  /* B3[1] = B_xk[2][1]; */
  /* B0Hat[0] = B_x0[3][0]; */
  /* B0Hat[1] = B_x0[3][1]; */
  /* BHat[0] = B_xk[3][0]; */
  /* BHat[1] = B_xk[3][1]; */
  /* x3[0] = points_fan[3][0]; */
  /* x3[1] = points_fan[3][1]; */

  /* hermite_interpolationSpatial(0.15, x0, x3, B03, B3, y3Prime); */
  /* grad_hermite_interpolationSpatial(0.15, x0, x3, B03, B3, B3_mu); */
  /* printf("The coordinates of y3Prime:  %lf   %lf\n\n", y3Prime[0], y3Prime[1]); */

  /* double lambda4Min, lambda4Max, xlam4Min[2], xlam4Max[2]; */

  /* lambda4Min = lambda_fromt2(0, x0, B0Hat, y3Prime, xHat, BHat, 0.000001, 50); */

  /* printf("Minimum lambda4: %lf\n", lambda4Min); */

  /* hermite_interpolationSpatial(lambda4Min, x0, xHat, B0Hat, BHat, xlam4Min); */

  /* printf("Coordinates of minimum %lf  %lf\n", xlam4Min[0], xlam4Min[1]); */

  /* lambda4Max = lambda_fromt1(1, x0, B0Hat, y3Prime, B3_mu, xHat, BHat, 0.000001, 50); */

  /* printf("Maximum lambda4: %lf\n", lambda4Max); */

  /* hermite_interpolationSpatial(lambda4Max, x0, xHat, B0Hat, BHat, xlam4Max); */

  /* printf("Coordinates of minimum %lf  %lf\n", xlam4Max[0], xlam4Max[1]); */
  
}

