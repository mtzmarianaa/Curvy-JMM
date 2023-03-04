/*
Optimization method for the 2D JMM
The approach is based on the triangle fan
 */

#include "optiFan.h"
#include "linAlg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>

void optiFan_alloc(optiFanS **optiFan) {
  *optiFan = malloc(sizeof(optiFanS));
  assert(*optiFan != NULL);
}

void optiFan_dealloc(optiFanS **optiFan) {
  free(*optiFan);
  *optiFan = NULL;
  assert(*optiFan == NULL);
}

void optiFan_init(optiFanS *optiFan, int nRegions, double x0[2], double T0, double x1[2], double T1, double xHat[2], double *indicesRef, double (*points_fan)[2], double (*B_x0)[2], double (*B_xk)[2]) {
  double gradTest_x0[2], xkMinx0[2], xk[2], B0k1[2], B0k[2], Bk1[2], Bk[2], xk1[2], dotB0k1B0k, dotBk1Bk;
  int i;
  optiFan->nRegions = nRegions;
  optiFan->x0[0] = x0[0];
  optiFan->x0[1] = x0[1];
  optiFan->T0 = T0;
  optiFan->x1[0] = x1[0];
  optiFan->x1[1] = x1[1];
  optiFan->T1 = T1;
  optiFan->xHat[0] = xHat[0];
  optiFan->xHat[1] = xHat[1];
  optiFan->points_fan = points_fan;
  optiFan->indicesRef = indicesRef;
  optiFan->types = malloc(nRegions*sizeof(int));
  optiFan->B_x0 = malloc((nRegions+1)*2*sizeof(double));
  optiFan->B_xk = malloc((nRegions+1)*2*sizeof(double));
  for (i = 0; i<(nRegions+1); i++){
    // we need to know if the gradients go from x0 to xk
    xk[0] = points_fan[i][0];
    xk[1] = points_fan[i][1];
    vec2_subtraction(xk, x0, xkMinx0);
    printf("For i= %d  xk: %lf %lf  and xk-x0: %lf  %lf\n", i, xk[0], xk[1], xkMinx0[0], xkMinx0[1]);
    gradTest_x0[0] = B_x0[i][0];
    gradTest_x0[1] = B_x0[i][1];
    printf("The initial B0k given is: %lf  %lf\n", gradTest_x0[0], gradTest_x0[1]);
    printf("The initial Bk given is: %lf  %lf\n", B_xk[i][0], B_xk[i][1]);
    if( dotProd(gradTest_x0, xkMinx0)>=0 ){
      // direction we want x0->xk
      printf("   We don't have a problem with the original direction of the tangents\n");
      optiFan->B_x0[i][0] = gradTest_x0[0];
      optiFan->B_x0[i][1] = gradTest_x0[1];
      optiFan->B_xk[i][0] = B_xk[i][0];
      optiFan->B_xk[i][1] = B_xk[i][1];
    }
    else{
      // direction we don't want xk->x0
      printf("   We have a problem\n");
      optiFan->B_x0[i][0] = -gradTest_x0[0];
      optiFan->B_x0[i][1] = -gradTest_x0[1];
      optiFan->B_xk[i][0] = -B_xk[i][0];
      optiFan->B_xk[i][1] = -B_xk[i][1];
    }
  }
  // now we compute the type of triangle we have
  printf("\n\n");
  for (i = 0; i<nRegions; i++){
    Bk[0] = optiFan->B_xk[i][0];
    Bk[1] = optiFan->B_xk[i][1];
    Bk1[0] = optiFan->B_xk[i+1][0];
    Bk1[1] = optiFan->B_xk[i+1][1];
    B0k[0] = optiFan->B_x0[i][0];
    B0k[1] = optiFan->B_x0[i][1];
    B0k1[0] = optiFan->B_x0[i+1][0];
    B0k1[1] = optiFan->B_x0[i+1][1];
    printf("For i=%d  Bk= %lf  %lf   B(k+1)=   %lf  %lf\n", i, Bk[0], Bk[1], Bk1[0], Bk1[1]);
    printf("          B0k= %lf  %lf   B0(k+1)=   %lf  %lf\n", B0k[0], B0k[1], B0k1[0], B0k1[1]);
    dotB0k1B0k = dotProd(B0k1, B0k);
    dotBk1Bk = dotProd(Bk1, Bk);
    printf("Dot product B0k1 and B0k  %lf\n", dotB0k1B0k);
    printf("Dot product Bk1 and Bk  %lf\n", dotBk1Bk);
    if( dotB0k1B0k <=0 & dotBk1Bk >=0 ){
      optiFan->types[i] = 3;
    }
    else if(dotB0k1B0k>0 & dotBk1Bk<0){
      optiFan->types[i] = 4;
    }
    else if(dotB0k1B0k >=0 & dotBk1Bk >=0){
      xk1[0] = points_fan[i+1][0];
      xk1[1] = points_fan[i+1][1];
      vec2_subtraction(xk1, x0, xkMinx0);
      if( dotProd(xkMinx0, B0k) >= 0){
	optiFan->types[i] = 1;
      }
      else{
	optiFan->types[i] = 2;
      }
    }
    else{
      // this case should not happen, it means there is an inflection point
      assert(false); // SHOULD NEVER HAPPEN!
      printf("For the point %d  with Bk  %lf  %lf   and B0k  %lf %lf: \n", i, Bk[0], Bk[1], B0k[0], B0k[1]);
      printf("                   and Bk1  %lf  %lf   and B0k1  %lf %lf: \n", Bk1[0], Bk1[1], B0k1[0], B0k1[1]);
    }
    printf(" For the internal triangle %d  the type of curvy triangle is: %d\n\n", i, optiFan->types[i]);
  }
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

// functions to find lambdaMin and lambdaMax for types 1,2, and 4

double t1_ofLam(double lambda, double x0[2], double B0_k1[2], double ykPrime[2], double Bk_mu[2], double x_k1[2], double B_k1[2]) {
  // line through yk' with slope Bk(mu)
  double y_k1[2];
  hermite_interpolationSpatial(lambda, x0, x_k1, B0_k1, B_k1, y_k1);
  return Bk_mu[0]*(y_k1[1] - ykPrime[1]) - Bk_mu[1]*(y_k1[0] - ykPrime[0]);
}

double t1Prime_ofLam(double lambda, double x0[2], double B0_k1[2], double Bk_mu[2], double x_k1[2], double B_k1[2]) {
  // derivative with respect to lambda of the line through yk' with slope Bk(mu)
  double B_lam[2];
  grad_hermite_interpolationSpatial(lambda, x0, x_k1, B0_k1, B_k1, B_lam);
  return Bk_mu[0]*B_lam[1] - Bk_mu[1]*B_lam[0];
}

double backTr_t1(double alpha0, double d, double lambda, double x0[2], double B0_k1[2], double ykPrime[2], double Bk_mu[2], double x_k1[2], double B_k1[2]) {
  // backtracking to find suitable step size for t1
  double t1_prev, t1_cur, alpha;
  int i = 0;
  alpha = alpha0;
  // evaluating objective function
  t1_prev = t1_ofLam(lambda, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1);
  t1_cur = t1_ofLam(lambda - alpha*d, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1);
  while( fabs(t1_prev) <= fabs(t1_cur) & i < 10){
    alpha = alpha*0.5;
    t1_cur = t1_ofLam(lambda - alpha*d, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1);
    i++;
  }
  if( fabs(t1_prev) <= fabs(t1_cur) ){
    alpha = 0;
  }
  return alpha;
}

double lambda_fromt1(double lambda0, double x0[2], double B0_k1[2], double ykPrime[2], double Bk_mu[2], double x_k1[2], double B_k1[2], double tol, int maxIter) {
  // newtons method to find lambda such that t1(lambda) = 0
  double lambda, t1_cur, t1Prime, alpha, d;
  int i = 0;
  lambda = lambda0;
  t1_cur = t1_ofLam(lambda, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1);
  while( fabs(t1_cur) > tol & i < maxIter & alpha > 0){
    t1Prime = t1Prime_ofLam(lambda, x0, B0_k1, Bk_mu, x_k1, B_k1);
    d = t1_cur/t1Prime;
    alpha = backTr_t1(0.1, d, lambda, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1);
    lambda = lambda - alpha*d;
    t1_cur = t1_ofLam(lambda, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1);
    i++;
  }
  return lambda;
}

double t2_ofLam(double lambda, double x0[2], double B0_k1[2], double ykPrime[2], double x_k1[2], double B_k1[2]) {
  // line through yk' with slope B(k+1)(lambda)
  double y_k1[2], B_lam[2];
  hermite_interpolationSpatial(lambda, x0, x_k1, B0_k1, B_k1, y_k1);
  grad_hermite_interpolationSpatial(lambda, x0, x_k1, B0_k1, B_k1, B_lam);
  return B_lam[0]*(y_k1[1] - ykPrime[1]) - B_lam[1]*(y_k1[0] - ykPrime[0]);
}

double t2Prime_ofLam(double lambda, double x0[2], double B0_k1[2], double ykPrime[2], double x_k1[2], double B_k1[2]) {
  // derivative with respect to lambda of the line through yk' with slope B(k+1)(lambda)
  double y_k1[2], B_lam[2], derB_lam[2];
  hermite_interpolationSpatial(lambda, x0, x_k1, B0_k1, B_k1, y_k1);
  grad_hermite_interpolationSpatial(lambda, x0, x_k1, B0_k1, B_k1, B_lam);
  secondDer_hermite_interpolationSpatial(lambda, x0, x_k1, B0_k1, B_k1, derB_lam);
  return derB_lam[0]*(y_k1[1] - ykPrime[1]) - derB_lam[1]*(y_k1[0] - ykPrime[0]);
}

double backTr_t2(double alpha0, double d, double lambda, double x0[2], double B0_k1[2], double ykPrime[2], double x_k1[2], double B_k1[2]) {
  // backtracking to find suitable step size for t2
  double t2_prev, t2_cur, alpha;
  int i = 0;
  alpha = alpha0;
  // evaluating objective function
  t2_prev = t2_ofLam(lambda, x0, B0_k1, ykPrime, x_k1, B_k1);
  t2_cur = t2_ofLam(lambda - alpha*d, x0, B0_k1, ykPrime, x_k1, B_k1);
  while( fabs(t2_prev) < fabs(t2_cur) & i < 10){
    alpha = alpha*0.5;
    t2_cur = t2_ofLam(lambda - alpha*d, x0, B0_k1, ykPrime, x_k1, B_k1);
    i++;
  }
  if( fabs(t2_prev) < fabs(t2_cur) ){
    alpha = 0;
  }
  return alpha;
}

double lambda_fromt2(double lambda0, double x0[2], double B0_k1[2], double ykPrime[2], double x_k1[2], double B_k1[2], double tol, int maxIter) {
  // Newtons method to find lambda such that t2(lambda) = 0
  double lambda, t2_cur, t2Prime, alpha, d;
  int i = 0;
  lambda = lambda0;
  t2_cur = t2_ofLam(lambda, x0, B0_k1, ykPrime, x_k1, B_k1);
  while( fabs(t2_cur) > tol & i < maxIter & alpha > 0){
    t2Prime = t2Prime_ofLam(lambda, x0, B0_k1, ykPrime, x_k1, B_k1);
    d = t2_cur/t2Prime;
    alpha = backTr_t2(1, d, lambda, x0, B0_k1, ykPrime, x_k1, B_k1);
    lambda = lambda - alpha*d;
    t2_cur = t2_ofLam(lambda, x0, B0_k1, ykPrime, x_k1, B_k1);
    i++;
  }
  return lambda;
}

// these are the functions used to project back on the feasible set
// but they depend on the type of curvy triangle we are dealing with

void projectBack_type1(double lambdak1, double yk1[2], double ykPrime[2], double x0[2], double B0_k1[2], double Bk_mu[2], double Bk_mu_perp[2], double x_k[2], double x_k1[2], double B_k1[2], double tol, int maxIter) {
  // given all this information project back yk1. Here we assume that this triangle is
  // a type 1 curvy triangle
  if(lambdak1 > 1){
    lambdak1 = 1;
    yk1[0] = x_k1[0];
    yk1[1] = x_k1[1];
  }
  else if(lambdak1 < 0){
    lambdak1 = 0;
    yk1[0] = x0[0];
    yk1[1] = x0[1];
  }
  // we need to know if yk1 is in the feasible set computed from ykPrime
  double yk1MinykPrime[2], dotTest, lambdat1;
  vec2_subtraction(yk1,ykPrime, yk1MinykPrime);
  dotTest = dotProd(yk1MinykPrime, Bk_mu_perp);
  if(dotTest < 0){
    // then we need to find lambdaMax/lambdaMin and project back
    lambdat1 = lambda_fromt1(lambdak1, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1, tol, maxIter);
    printf("Lambda from t1 found: %lf\n", lambdat1);
    printf("Parameters used in lambda_fromt1:  \n");
    printf("    B0_k1:  %lf   %lf\n", B0_k1[0], B0_k1[1]);
    printf("    ykPrime:  %lf   %lf\n", ykPrime[0], ykPrime[1]);
    printf("    Bk_mu:  %lf   %lf\n", Bk_mu[0], Bk_mu[1]);
    printf("    x_k1:  %lf   %lf\n", x_k1[0], x_k1[1]);
    printf("    B_k1:  %lf   %lf\n", B_k1[0], B_k1[1]);
    if(lambdat1 < 0){
      lambdat1 = 0;
    }
    else if(lambdat1 > 1){
      lambdat1 = 1;
    }
    // compute yk1 from this newly found lambda
    hermite_interpolationSpatial(lambdat1, x0, x_k1, B0_k1, B_k1, yk1);
  }
}

void projectBack_type2(double lambdak1, double yk1[2], double ykPrime[2], double x0[2], double B0_k1[2], double Bk_mu[2], double Bk1_lam_perp[2], double x_k[2], double x_k1[2], double B_k1[2], double tol, int maxIter) {
  // given all this information project back yk1. Here we are going to assume that this
  // is a type 2 curvy triangle
  if(lambdak1 > 1){
    lambdak1 = 1;
    yk1[0] = x_k1[0];
    yk1[1] = x_k1[1];
  }
  else if(lambdak1 < 0){
    lambdak1 = 0;
    yk1[0] = x0[0];
    yk1[1] = x0[1];
  }
  // we need to know if yk1 is in the feasible set computed from ykPrime
  double yk1MinykPrime[2], dotTest, lambdat2;
  vec2_subtraction(yk1, ykPrime, yk1MinykPrime);
  dotTest = dotProd(yk1MinykPrime, Bk1_lam_perp);
  if( dotTest<0){
    // if means that we need to compute lambdaMax/lambdaMin from t2 and then project back
    lambdat2 = lambda_fromt2(lambdak1, x0, B0_k1, ykPrime, x_k1, B_k1, tol, maxIter);
    if(lambdat2 < 0){
      lambdat2 = 0;
    }
    else if(lambdat2 > 1){
      lambdat2 = 1;
    }
    // compute yk1 from this newly found lambda
    hermite_interpolationSpatial(lambdat2, x0, x_k1, B0_k1, B_k1, yk1);
  }
}

void projectBack_type4(double lambdak1, double yk1[2], double ykPrime[2], double x0[2], double B0_k1[2], double Bk_mu[2], double B_k_mu_perp[2], double B_k1_lam_perp[2], double x_k[2], double x_k1[2], double B_k1[2], double tol, double maxIter) {
  // given all of this information project back yk1. Here we are going to assume
  // that this is a type 4 curvy triangle
  if(lambdak1 > 1){
    lambdak1 = 1;
    yk1[0] = x_k1[0];
    yk1[1] = x_k1[1];
  }
  else if(lambdak1 < 0){
    lambdak1 = 0;
    yk1[0] = x0[0];
    yk1[1] = x0[1];
  }
  // then we need to find lambdaMin with t1 and lambdaMax with t2 (if necessary)
  double dotTestMin, dotTestMax, yk1MinykPrime[2], lambdat1, lambdat2;
  vec2_subtraction(yk1, ykPrime, yk1MinykPrime);
  dotTestMin = dotProd(yk1MinykPrime, B_k_mu_perp);
  dotTestMax = dotProd(yk1MinykPrime, B_k1_lam_perp);
  if( dotTestMin < 0){
    // meaning we need to find lambdaMin
    lambdat1 = lambda_fromt1(lambdak1, x0, B0_k1, ykPrime, Bk_mu, x_k1, B_k1, tol, maxIter);
    if(lambdat1 < 0){
      lambdat1 = 0;
    }
    else if(lambdat1 > 1){
      lambdat1 = 1;
    }
    // compute yk1 from this newly found lambda
    hermite_interpolationSpatial(lambdat1, x0, x_k1, B0_k1, B_k1, yk1);
  }
  else if( dotTestMax){
    // if means that we need to compute lambdaMax/lambdaMin from t2 and then project back
    lambdat2 = lambda_fromt2(lambdak1, x0, B0_k1, ykPrime, x_k1, B_k1, tol, maxIter);
    if(lambdat2 < 0){
      lambdat2 = 0;
    }
    else if(lambdat2 > 1){
      lambdat2 = 1;
    }
    // compute yk1 from this newly found lambda
    hermite_interpolationSpatial(lambdat2, x0, x_k1, B0_k1, B_k1, yk1);
  }
}




