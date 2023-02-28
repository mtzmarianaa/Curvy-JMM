/*
Optimization method for the 2D JMM
The approach is based on the triangle fan
 */

#include "optiFan.h"
#include "linAlg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void optiFan_alloc(optiFanS **optiFan) {
  *optiFan = malloc(sizeof(optiFanS));
  assert(*optiFan != NULL);
}

void optiFan_dealloc(optiFanS **optiFan) {
  free(*optiFan);
  *optiFan = NULL;
  assert(*optiFan == NULL);
}

void optiFan_init(optiFanS *optiFan, int nRegions, double x0[2], double T0; double x1[2], double T1, double xHat[2], double *(points_fan)[2], double (*boundary_x0)[2], double *indicesRef) {
  optiFan->nRegions = nRegions;
  optiFan->x0 = x0;
  optiFan->T0 = T0;
  optiFan->x1 = x1;
  optiFan->T1 = T1;
  optiFan->xHat = xHat;
  optiFan->points_fan = points_fan;
  optiFan->boundary_x0 = boundary_x0;
  optiFan->indicesRef = indicesRef;
}

// starting the section with auxiliary functions

void hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double interpolation[2]) {
  // general hermite interpolation evaluated at param for the boundary
  double param2, param3;
  double coef_0[2], coef_1[2], coef_grad0[2], coef_grad1[2];
  param2 = param*param;
  param3 = param2*param;
  scalar_times_2vec(2*param3 - 3*param2 + 1, from, coef_0);
  scalar_times_2vec(param3 - 2*param2 + param, grad_from, coef_grad0);
  scalar_times_2vec(-2*param3 + 3*param2, to, coef_1);
  scalar_times_2vec(param3 - param2, grad_to, coef_grad1);
  double sumpoints[2], sumgrads[2];
  vec2_addition(coef_0, coef_1, sumpoints);
  vec2_addition(coef_grad0, coef_grad1, sumgrads);
  vec2_addition(sumpoints, sumgrads, interpolation);
}

void grad_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double gradient[2]) {
  // gradient of the general hermite interpolation evaluated at param for the boundary (i.e. gradBmu)
  double param2;
  param2 = param*param;
  double twofrom[2], twoto[2], threefrom[2], twogradfrom[2], threeto[2];
  scalar_times_2vec(2, from, twofrom); // 2x0
  scalar_times_2vec(-2, to, twoto); // -2x1
  scalar_times_2vec(-3, from, threefrom); // -3x0
  scalar_times_2vec(-2, grad_from, twogradfrom); // -2B0
  scalar_times_2vec(3, to, threeto); // 3x1
  double sum1[2], sum2[2], sum3[2], sum4[2];
  vec2_addition(twofrom, grad_from, sum1); // 2x0 + B0
  vec2_addition(twoto, grad_to, sum2); // -2x1 + B1
  vec2_addition(threefrom, twogradfrom, sum3); // -3x0 -2B0
  vec2_subtraction(threeto, grad_to, sum4); // -3x0 - B0
  double coefparam2[2], coefparam[2]; 
  vec2_addition(sum1, sum2, coefparam2);
  vec2_addition(sum3, sum4, coefparam);
  double par1[2], par2[2], par3[2];
  scalar_times_2vec(3*param2, coefparam2, par1);
  scalar_times_2vec(2*param, coefparam, par2);
  vec2_addition(par2, grad_from, par3);
  vec2_addition(par1, par3, gradient);
}

void secondDer_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double secondDer[2]) {
  // second derivative with respect to param of the general hermite interpolation evaluated at param for the boundary
  double twofrom[2], twoto[2], threefrom[2], twogradfrom[2], threeto[2];
  scalar_times_2vec(2, from, twofrom);
  scalar_times_2vec(-2, to, twoto);
  scalar_times_2vec(-3, from, threefrom);
  scalar_times_2vec(-2, grad_from, twogradfrom);
  scalar_times_2vec(3, to, threeto);
  double sum1[2], sum2[2], sum3[2], sum4[2];
  vec2_addition(twofrom, grad_from, sum1);
  vec2_addition(twoto, grad_to, sum2);
  vec2_addition(threefrom, twogradfrom, sum3);
  vec2_subtraction(threeto, grad_to, sum4);
  double coefparam2[2], coefparam[2];
  vec2_addition(sum1, sum2, coefparam2);
  vec2_addition(sum3, sum4, coefparam);
  double par1[2], par2[2], par3[2];
  scalar_times_2vec(6*param, coefparam2, par1);
  vec2_addition(coefparam, par1, secondDer);
}

double arclength_hermiteSimpson(double a, double b, double from[2], double to[2], double grad_from[2], double grad_to[2]) {
  // using simpsons rule to calculate the arclength for the parameter being from a to b
  // uses grad_hermite_interpolationSpatial of course
  double gradient_a[2], gradient_mid[2], gradient_b[2], norm_a, norm_mid, norm_b;
  grad_hermite_interpolationSpatial(a, from, to, grad_from, grad_to, gradient_a);
  grad_hermite_interpolationSpatial((a+b)/2, from, to, grad_from, grad_to, gradient_mid);
  grad_hermite_interpolationSpatial(b, from, to, grad_from, grad_to, gradient_b);
  norm_a = l2norm(gradient_a);
  norm_mid = l2norm(gradient_mid);
  norm_b = l2norm(gradient_b);
  return (norm_a + 4*norm_mid + norm_b)*(b-a)/6;
}

double hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]) {
  // hermite interpolation for the value of the eikonal
  double param2, param3;
  param2 = param*param;
  param3 = param2*param;
  double coef_gradA[2], coef_gradB[2], sumgrads[2], xBminxA[2];
  scalar_times_2vec(param3 - 2*param2 + param, gradA, coef_gradA);
  scalar_times_2vec(param3 - param2, gradB, coef_gradB);
  vec2_addition(coef_gradA, coef_gradB, sumgrads);
  vec2_subtraction(xB, xA, xBminxA);
  return (2*param3 - 3*param2 + 1)*TA + (-2*param3 + 3*param2)*TB + dotProd(xBminxA, sumgrads);
}

double der_hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]) {
  // derivative with respect of param for the hermite interpolation for the value of the eikonal
  double param2;
  param2 = param*param;
  double coef0[2], coef1[2], sumGrads[2];
  scalar_times_2vec(3*param2 - 4*param + 1, gradA, coef0);
  scalar_times_2vec(3*param2 - 2*param, gradB, coef1);
  vec2_addition(coef0, coef1, sumGrads);
  double xBminxA[2];
  vec2_subtraction(xB, xA, xBminxA);
  return (6*param2 - 6*param)*TA + (-6*param2 + 6*param)*TB + dotProd(xBminxA,sumGrads); 
}

double t_ofMu(double mu, double x0[2], double ykPrime[2], double Bk_mu[2], double x_k1[2], double B_k1[2]){
  // line through yk' with slope Bk(mu)
  double 
}


