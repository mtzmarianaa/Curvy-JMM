
/*
Optimization methods for the 2D FMM
*/

#include "opti_method.h"
#include "linAlg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void linearInterpolation(double param, double from[2], double to[2], double interpolation[2]){
  // useful for xlambda = (1-lambda)*x0 + lambda*x1
  double coeffrom[2], coefto[2];
  scalar_times_2vec(1-param, from, coeffrom);
  scalar_times_2vec(param, to, coefto);
  vec2_addition(coeffrom, coefto, interpolation);
}

void der_linearInterpolation(double param, double from[2], double to[2], double der_interpolation[2]){
  // derivative with respect of lambda of (1-lambda)*x0 + lambda*x1
  vec2_subtraction(to, from, der_interpolation);
}

void hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double interpolation[2]){
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

void grad_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double gradient[2]){
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

void secondDer_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double secondDer[2]){
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

double arclength_hermiteSimpson(double a, double b, double from[2], double to[2], double grad_from[2], double grad_to[2]){
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

double hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]){
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

double der_hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]){
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



double der_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef) {
  // derivative with respect of lambda from the function to minimize for the update in which the
  // segment x0x1 is on the boundary but xHat if fully contained in a region with index indexRef
  // segments x0xHat and x1xHat are fully contained in a region 
  double xLam[2], gradxLam[2], der_Tlam, xHatminxLam[2], normxHatminxLam;
  hermite_interpolationSpatial(lambda, x0, x1, B0, B1, xLam);
  grad_hermite_interpolationSpatial(lambda, x0, x1, B0, B1, gradxLam);
  der_Tlam = der_hermite_interpolationT(lambda, x0, x1, T0, T1, grad0, grad1);
  vec2_subtraction(xLam, xHat, xHatminxLam);
  normxHatminxLam = l2norm(xHatminxLam);
  // Putting everything together
  return der_Tlam + indexRef*dotProd(xHatminxLam, gradxLam)/normxHatminxLam;
}

double backTr_fromEdge(double alpha0, double d, double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef){
  // backtracking method for projected gradient descent from updates from the edge of the domain
  double f_prev, f_cur, alpha;
  int i = 0;
  alpha = alpha0;
  // EVALUATING THE OBJECTIVE FUNCTION
  f_prev = fobjective_fromEdge(lambda, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  f_cur = fobjective_fromEdge(lambda - alpha*d, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  while(f_prev <= f_cur & i < 10){
    alpha = alpha*0.5;
    f_cur = fobjective_fromEdge(lambda - alpha*d, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
    i ++;
  }
  return alpha;
}

double fobjective_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef){

  // T(xlambda) + indexRef*norm(xHat - xLam)

  double Tlam, xLam[2];

  Tlam = hermite_interpolationT(lambda, x0, x1, T0, T1, grad0, grad1);
  hermite_interpolationSpatial(lambda, x0, x1, B0, B1, xLam);

  // put everything together
  double xHatminxLam[2], normxHatminxLam;
  vec2_subtraction(xHat, xLam, xHatminxLam);
  normxHatminxLam = l2norm(xHatminxLam);

  return Tlam + indexRef*normxHatminxLam;
}


double projectedGradient_fromEdge(double lambda0, double lambdaMin, double lambdaMax, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double tol, int maxIter, double indexRef){
  // projected gradient descent for an update in which the segment x0x1 is on the boundary but xHat
  // if fully contained in a region with index indexRef
  double grad_cur, grad_prev, step, alpha, lam_prev, lam_cur, test;
  int i;
  i = 1;
  alpha = 0.25;
  lam_prev = lambda0;
  grad_cur = der_fromEdge(lambda0, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  grad_prev = der_fromEdge(lambda0, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  if(fabs(grad_cur) > tol){
    test = lam_prev - alpha*grad_cur;
  }
  else{
    test = lam_prev;
  }
  if(test>lambdaMax){
    lam_cur = lambdaMax;
  }
  else if(test<lambdaMin){
    lam_cur = lambdaMin;
  }
  else{
    lam_cur = test;
  }
  grad_cur = der_fromEdge(lam_cur, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);

  while(i<maxIter & fabs(grad_cur)>tol & fabs(lam_cur - lam_prev)>0) {
    alpha = backTr_fromEdge(0.25, grad_cur, lam_cur, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
    test = lam_cur - alpha*grad_cur;
    if(test<lambdaMin){
      test = lambdaMin;
    }
    if(test>lambdaMax){
      test = lambdaMax;
    }
    grad_prev = grad_cur;
    lam_prev = lam_cur;
    lam_cur = test;
    grad_cur = der_fromEdge(lam_cur, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
    i ++;
  }
  
  return lam_cur;
}

double der_freeSpace(double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef){
  double xLam[2], gradxLam[2], der_Tlam, xHatminxLam[2], normxHatminxLam;
  linearInterpolation(lambda, xA, xB, xLam);
  der_linearInterpolation(lambda, xA, xB, gradxLam);
  der_Tlam = der_hermite_interpolationT(lambda, xA, xB, TA, TB, gradA, gradB);
  vec2_subtraction(xLam, xHat, xHatminxLam);
  normxHatminxLam = l2norm(xHatminxLam);
  return der_Tlam + indexRef*dotProd(xHatminxLam, gradxLam)/normxHatminxLam;
}

double backTr_freeSpace(double alpha0, double d, double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef){
  double f_prev, f_cur, alpha;
  int i = 0;
  alpha = alpha0;
  // EVALUATING THE OBJECTIVE FUNCTION
  f_prev = fobjective_freeSpace(lambda, TA, gradA, TB, gradB, xA, xB, xHat, indexRef);
  f_cur = fobjective_freeSpace(lambda - alpha*d, TA, gradA, TB, gradB, xA, xB, xHat, indexRef);
  while(f_prev <= f_cur & i < 10){
    alpha = alpha*0.5;
    f_cur = fobjective_freeSpace(lambda - alpha*d, TA, gradA, TB, gradB, xA, xB, xHat, indexRef);
    i ++;
  }
  return alpha;
}

double fobjective_freeSpace(double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef){
  double Tlam, xLam[2];
  Tlam = hermite_interpolationT(lambda, xA, xB, TA, TB, gradA, gradB);
  linearInterpolation(lambda, xA, xB, xLam);
  // put everything together
  double xHatminxLam[2], normxHatminxLam;
  vec2_subtraction(xHat, xLam, xHatminxLam);
  normxHatminxLam = l2norm(xHatminxLam);

  return Tlam + indexRef*normxHatminxLam;
}

double projectedGradient_freeSpace(double lambda0, double lambdaMin, double lambdaMax, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double tol, int maxIter, double indexRef){
  // two point optimization problem. xA and xHat are on the boundary and xB is fully contained in one region.
  // With usual notation xA could be either x0 or x1 and xB could be either x0 or x1 (this is more abstract)
  // THIS IS A PROJECTED GRADIENT DESCENT METHOD
  double grad_cur, grad_prev, step, alpha, lam_prev, lam_cur, test;
  int i =1;
  alpha = 1;
  lam_prev = lambda0;
  grad_cur = der_freeSpace(lambda0, TA, gradA, TB, gradB, xA, xB, xHat, indexRef);
  grad_prev = grad_cur;
  if(fabs(grad_cur) > tol){
    test = lam_prev - alpha*grad_cur;
  }
  else{
    test = lam_prev;
  }
  if(test>lambdaMax){
    test = lambdaMax;
  }
  if(test<lambdaMin){
    test = lambdaMin;
  }
  lam_cur = test;
  grad_cur = der_freeSpace(lam_cur, TA, gradA, TB, gradB, xA, xB, xHat, indexRef);

  while( i<maxIter & fabs(grad_cur)>tol & fabs(lam_cur - lam_prev)>0){
    alpha = backTr_freeSpace(1, grad_cur, lam_cur, TA, gradA, TB, gradB, xA, xB, xHat, indexRef);
    //printf("\n\nIteration %d\n", i);
    //printf("Step size %lf   with direction  %lf, hence change is %lf\n", alpha, -grad_cur, -alpha*grad_cur);
    test = lam_cur - alpha*grad_cur;
    if(test > lambdaMax){
      test = lambdaMax;
    }
    if(test < lambdaMin){
      test = lambdaMin;
    }

    grad_prev = grad_cur;
    lam_prev = lam_cur;
    lam_cur = test;
    grad_cur = der_freeSpace(lam_cur, TA, gradA, TB, gradB, xA, xB, xHat, indexRef);
    //printf("Iteration %d   with lam_prev %lf  and lam_cur %lf   with objective value: %lf   and derivative  %lf \n\n\n", i, lam_prev, lam_cur, fobjective_freeSpace(lam_cur, TA, gradA, TB, gradB, xA, xB, xHat, indexRef), grad_cur);
    
    i++;
  }

  return lam_cur;
}

void grad_twoStep(double gradient[2], double lambda, double mu, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02) {
  // gradient for the function to minimize for a two step update (when there is a change in region)
  double mu2, mu3, lambda2;
  mu2 = mu*mu;
  mu3 = mu2*mu;
  lambda2 = lambda*lambda;
  /////// partial with respect to lambda
  // first time gathering terms
  double x1Minx0[2], x0InxLamMinxMu[2], x1inxLamMinxMu[2], B0inxLamMinxMu[2], x2InxLamMinxMu[2], B2inxLamMinxMu[2];
  vec2_subtraction(x1, x0, x1Minx0);
  scalar_times_2vec(-lambda-2*mu3 + 3*mu2, x0, x0InxLamMinxMu);
  scalar_times_2vec(lambda, x1, x1inxLamMinxMu);
  scalar_times_2vec(mu3 - 2*mu2 + mu, B0, B0inxLamMinxMu);
  scalar_times_2vec(-2*mu3 + 3*mu2, x2, x2InxLamMinxMu);
  scalar_times_2vec(mu3 - mu2, B2, B2inxLamMinxMu);
  // second time gatherieng terms
  double sum1[2], sum2[2], sum3[2], xLamMinxMu[2], t, normxLamMinxMu, coefgrad0[2], coefgrad1[2], sumgrads[2], dot1;
  vec2_addition(x0InxLamMinxMu, x1inxLamMinxMu, sum1);
  vec2_addition(B0inxLamMinxMu, x2InxLamMinxMu, sum2);
  vec2_addition(sum2, B2inxLamMinxMu, sum3);
  vec2_subtraction(sum1, sum3, xLamMinxMu);
  t = dotProd(x1Minx0, xLamMinxMu);
  normxLamMinxMu = l2norm(xLamMinxMu);
  scalar_times_2vec(3*lambda2 - 4*lambda + 1, grad0, coefgrad0);
  scalar_times_2vec(3*lambda2 - 2*lambda, grad1, coefgrad1);
  vec2_addition(coefgrad0, coefgrad1, sumgrads);
  dot1 = dotProd(x1Minx0, sumgrads);
  // then the partial is
  gradient[0] =(6*lambda2 - 6*lambda)*T0 + (-6*lambda2 + 6*lambda)*T1 + dot1 + indexRef_01*t/normxLamMinxMu;
  /////// partial with respect to mu
  // first time gathering terms
  double der_coefx0[2], der_coefB0[2], der_coefx2[2], der_coefB2[2];
  scalar_times_2vec(-6*mu2 + 6*mu, x0, der_coefx0);
  scalar_times_2vec(-3*mu2 + 4*mu - 1, B0, der_coefB0);
  scalar_times_2vec(6*mu2 - 6*mu, x2, der_coefx2);
  scalar_times_2vec(-3*mu2 + 2*mu, B2,  der_coefB2);
  // second time gathering terms
  double coefx0[2], coefB0[2], coefx2[2], coefB2[2];
  scalar_times_2vec(-2*mu3 + 3*mu2 - 1, x0, coefx0);
  scalar_times_2vec(-mu3 + 2*mu2 - mu, B0, coefB0);
  scalar_times_2vec(2*mu3 - 3*mu2, x2, coefx2);
  scalar_times_2vec(-mu3 + mu2, B2, coefB2);
  // third time gathering terms
  double der_sum1[2], der_sum2[2], derMuxMu[2], sum4[2], sum5[2], sum6[2], xLamMinxHat[2], normxLamMinxHat;
  vec2_addition(der_coefx0, der_coefB0, der_sum1);
  vec2_addition(der_coefx2,  der_coefB2, der_sum2);
  vec2_addition(der_sum1, der_sum2, derMuxMu);
  vec2_addition(xHat, coefx0, sum4);
  vec2_addition(coefB0, coefx2, sum5);
  vec2_addition(sum5, coefB2, sum6);
  vec2_addition(sum4, sum6, xLamMinxHat);
  normxLamMinxHat = l2norm(xLamMinxHat);
  // fourth time gathering terms
  double t1, t2;
  t1 = dotProd( derMuxMu, xLamMinxMu);
  t2 = dotProd( derMuxMu, xLamMinxHat);
  // everything together
  gradient[1] = indexRef_01*t1/normxLamMinxMu + indexRef_02*t2/normxLamMinxHat;
}


double backTr_TwoStep(double alpha0, double d[2], double lambda, double mu, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02) {
  // backtracking to compute a step length in the two step update
  double f_prev, f_cur, alpha;
  int i = 0;
  alpha = alpha0;
  // EVALUATING THE OBJECTIVE FUNCTION
  f_prev = fobjective_TwoStep(lambda, mu, T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02);
  //printf("Objective function before %lf  with lambda %lf  and mu  %lf\n", f_prev, lambda, mu);
  f_cur = fobjective_TwoStep(lambda - alpha*d[0], mu - alpha*d[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02);
  while(f_prev <= f_cur & i < 10 ){
    alpha = alpha*0.5;
    f_cur = fobjective_TwoStep(lambda - alpha*d[0], mu - alpha*d[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02);
    i ++;
  }
  if (f_prev <= f_cur){
    alpha = 0;
  }
  //printf("Objective function adter back tracking  %lf\n", fobjective_TwoStep(lambda - alpha*d[0], mu - alpha*d[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02) );
  //printf("with lambda %lf  and mu %lf\n", lambda-alpha*d[0], lambda-alpha*d[1]);
  return alpha;
}

double fobjective_TwoStep(double lambda, double mu, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02){
  // objective function for a two step update
  double lambda2, lambda3, mu2, mu3;
  lambda2 = lambda*lambda;
  lambda3 = lambda2*lambda;
  mu2 = mu*mu;
  mu3 = mu2*mu;
  // firt time gathering terms
  double x1Minx0[2], coefgrad0[2], coefgrad1[2], sumgrads[2], dotProd1;
  vec2_subtraction(x1, x0, x1Minx0);
  scalar_times_2vec(lambda3 - 2*lambda2 + lambda, grad0, coefgrad0);
  scalar_times_2vec(lambda3 - lambda2, grad1, coefgrad1);
  vec2_addition(coefgrad0, coefgrad1, sumgrads);
  dotProd1 = dotProd(x1Minx0, sumgrads);
  // second time gathering terms
  double x0InxLamMinxMu[2], x1inxLamMinxMu[2], B0inxLamMinxMu[2], x2InxLamMinxMu[2], B2inxLamMinxMu[2];
  scalar_times_2vec(-lambda-2*mu3 + 3*mu2, x0, x0InxLamMinxMu);
  scalar_times_2vec(lambda, x1, x1inxLamMinxMu);
  scalar_times_2vec(-mu3 + 2*mu2 - mu, B0, B0inxLamMinxMu);
  scalar_times_2vec(2*mu3 - 3*mu2, x2, x2InxLamMinxMu);
  scalar_times_2vec(-mu3 + mu2, B2, B2inxLamMinxMu);
  // third time gatherieng terms
  double sum1[2], sum2[2], sum3[2], xLamMinxMu[2], t, normxLamMinxMu, dot1;
  vec2_addition(x0InxLamMinxMu, x1inxLamMinxMu, sum1);
  vec2_addition(B0inxLamMinxMu, x2InxLamMinxMu, sum2);
  vec2_addition(sum2, B2inxLamMinxMu, sum3);
  vec2_addition(sum1, sum3, xLamMinxMu);
  normxLamMinxMu = l2norm(xLamMinxMu);
  // fourth time gathering terms
  double coefx0[2], coefB0[2], coefx2[2], coefB2[2];
  scalar_times_2vec(-2*mu3 + 3*mu2 - 1, x0, coefx0);
  scalar_times_2vec(-mu3 + 2*mu2 - mu, B0, coefB0);
  scalar_times_2vec(2*mu3 - 3*mu2, x2, coefx2);
  scalar_times_2vec(-mu3 + mu2, B2, coefB2);
  // fifth time gathering terms
  double sum4[2], sum5[2], sum6[2], xLamMinxHat[2], normxLamMinxHat;
  vec2_addition(xHat, coefx0, sum4);
  vec2_addition(coefB0, coefx2, sum5);
  vec2_addition(sum5, coefB2, sum6);
  vec2_addition(sum4, sum6, xLamMinxHat);
  normxLamMinxHat = l2norm(xLamMinxHat);
  
  return (2*lambda3 - 3*lambda2 + 1)*T0 + (-2*lambda3 + 3*lambda2)*T1 + dotProd1 + indexRef_01*normxLamMinxMu + indexRef_02*normxLamMinxHat;
  
}


void projectedGradient_TwoStep(double optimizers[2], double lambdaMin, double lambdaMax, double muMin, double muMax, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02, double tol, int maxIter) {
  // projected gradient descent for a two step update (when the index of refraction changes)
  double grad_cur[2], grad_prev[2], step, alpha, optimizers_prev[2], optimizers_cur[2], test[2], difStep[2];
  int i = 1;
  alpha = 0.1;
  optimizers[0] = lambdaMax;
  optimizers[1] = muMax;
  // compute the gradient
  grad_twoStep(grad_cur, optimizers[0], optimizers[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02);
  grad_prev[0] = grad_cur[0];
  grad_prev[1] = grad_cur[1];
  if (l2norm(grad_cur) > tol){
    test[0] = optimizers[0] - alpha*grad_cur[0];
    test[1] = optimizers[1] - alpha*grad_cur[1];
  }
  else{
    test[0] = optimizers[0];
    test[1] = optimizers[1];
  }
  if( test[0] > lambdaMax){
    test[0] = lambdaMax;
  }
  if(test[1] > muMax){
    test[1] = muMax;
  }
  if(test[0] < lambdaMin){
    test[0] = lambdaMin;
  }
  if(test[1] < muMin){
    test[1] = muMin;
  }
  // start the iteration
  optimizers_cur[0] = test[0];
  optimizers_cur[1] = test[1];
  grad_twoStep(grad_cur, optimizers_cur[0], optimizers_cur[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02);
  vec2_subtraction(optimizers_prev, optimizers_cur, difStep);
  while(i < maxIter & l2norm(grad_cur) > tol & l2norm(difStep) > 0){
    //printf("\n\nIteration %d\n", i);
    alpha = backTr_TwoStep(0.1, grad_cur, optimizers_cur[0], optimizers_cur[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02);
    //printf("Step size %lf  with direction  %lf  %lf, hence the change is  %lf   %lf\n", alpha, -grad_cur[0], -grad_cur[1], -alpha*grad_cur[0], - alpha*grad_cur[1]);
    test[0] = optimizers_cur[0] - alpha*grad_cur[0];
    test[1] = optimizers_cur[1] - alpha*grad_cur[1];
    //printf("Values before projecting back  %lf   %lf\n", test[0], test[1]);
    // project back if neccesary
    if( test[0] > lambdaMax){
      test[0] = lambdaMax;
    }
    if(test[1] > muMax){
      test[1] = muMax;
    }
    if(test[0] < lambdaMin){
      test[0] = lambdaMin;
    }
    if(test[1] < muMin){
      test[1] = muMin;
    }
    // if there is no better function value don't update
    if( fobjective_TwoStep(optimizers_cur[0], optimizers_cur[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02) < fobjective_TwoStep(test[0], test[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02) ) {
      test[0] = optimizers_cur[0];
      test[1] = optimizers_cur[1];
    }
    // update
    grad_prev[0] = grad_cur[0];
    grad_prev[1] = grad_cur[1];
    optimizers_prev[0] = optimizers_cur[0];
    optimizers_prev[1] = optimizers_cur[1];
    optimizers_cur[0] = test[0];
    optimizers_cur[1] = test[1];
    vec2_subtraction(optimizers_prev, optimizers_cur, difStep);
    grad_twoStep(grad_cur, optimizers_cur[0], optimizers_cur[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02);
    //printf("Iteration %d  with lambda_prev   %lf   and lambda_cur  %lf\n", i, optimizers_prev[0], optimizers_cur[0]);
    //printf("Iteration %d  with mu_prev   %lf   and mu_cur  %lf\n", i, optimizers_prev[1], optimizers_cur[1]);
    //printf("Objective value   %lf    and gradient %lf  %lf \n", fobjective_TwoStep(optimizers_cur[0], optimizers_cur[1], T0, grad0, T1, grad1, x0, x1, x2, xHat, B0, B2, indexRef_01, indexRef_02), grad_cur[0], grad_cur[1]);
    i ++;
  }
  optimizers[0] = optimizers_cur[0];
  optimizers[1] = optimizers_cur[1];
}


double t0_ofMu(double mu, double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]) {
  // function to define muMin (xB - xMu)T(gradMuPerp)
  double Bmu[2], xMu[2], Bmu_perp[2];
  grad_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, Bmu);
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
  Bmu_perp[0] = Bmu[1];
  Bmu_perp[1] = -Bmu[0];
  double xBminxMu[2];
  vec2_subtraction(xB, xMu, xBminxMu);
  return dotProd(xBminxMu, Bmu_perp);
}

double der_t0_ofMu(double mu, double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]){
  // derivative of t0_ofMu with respect to mu
  double xMu[2], derBmu[2], derBmu_perp[2];
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
  double xBminxMu[2];
  vec2_subtraction(xB, xMu, xBminxMu);
  secondDer_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, derBmu);
  derBmu_perp[0] = derBmu[1];
  derBmu_perp[1] = -derBmu[0];
  return dotProd(xBminxMu, derBmu_perp);
}

double findMumin_shootCr(double mu0, double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double tol, int maxIter){
  // Newtons method to find minimum mu of a shoot and creep type of update
  // note: this only makes sense if xR is neither xA or xB because if
  //       xR = xA or xR = xB then miMin is just 0
  double t0, der_t0, dif_mus, mu, mu_opt;
  int i = 1;
  t0 = t0_ofMu(mu0, xB, xHat, xR, BHat, BR);
  der_t0 = der_t0_ofMu(mu0, xB, xHat, xR, BHat, BR);
  dif_mus = t0/der_t0;
  mu = mu0 - dif_mus;
  while(fabs(dif_mus)>tol & i < maxIter){
    t0 = t0_ofMu(mu, xB, xHat, xR, BHat, BR);
    der_t0 = der_t0_ofMu(mu, xB, xHat, xR, BHat, BR);
    dif_mus = t0/der_t0;
    mu = mu - dif_mus;
    i++;
  }
  mu_opt = mu;
  if(mu>1){
    mu_opt = 1;
  }
  else if(mu<0){
    mu_opt = 0;
  }
  return mu_opt;
}

double t_ofMu(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]){
  // lambda = t(mu)
  double xMu[2], Bmu[2], Bmu_perp[2];
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
  grad_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, Bmu);
  double t, b;
  t = (xMu[0] - xA[0])*Bmu[1] - (xMu[1] - xA[1])*Bmu[0];
  b = (xA[0] - xB[0])*Bmu[1] - (xA[1] - xB[1])*Bmu[0];
  return -t/b;
}

double der_t_ofMu(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]){
  // derivative if t_ofMu with respect to mu
  double xMu[2], Bmu[2], Bmu_perp[2], derBmu[2], derBmu_perp[2];
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
  grad_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, Bmu);
  Bmu_perp[0] = Bmu[1];
  Bmu_perp[1] = -Bmu[0];
  secondDer_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, derBmu);
  derBmu_perp[0] = derBmu[1];
  derBmu_perp[1] = -derBmu[0];
  double xMuminxA[2], xBminxA[2];
  vec2_subtraction(xMu, xA, xMuminxA);
  vec2_subtraction(xB, xA, xBminxA);
  double b;
  b = dotProd(xBminxA, Bmu_perp);
  return (dotProd(xMuminxA, derBmu_perp)*dotProd(xBminxA, Bmu_perp) - dotProd(xBminxA, derBmu_perp)*dotProd(xMuminxA, Bmu_perp))/(b*b);
}

double backTr_find_minMu(double mu, double alpha0, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double tol, double maxIter){
  // backtracking algorithm to find the best step for finding the minimum mu
  //  (0<= mu <= 1) such that the segment x0xMu is parallel to Bmu)
  double mu_new, fEval, derEval, fEval_new, derEval_new, xMu[2], Bmu[2];
  double Bmu_perp[2], xMuminxB[2], derBmu[2], derBmu_perp[2], alpha;
  int i = 1;
  alpha = alpha0;
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
  grad_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, Bmu);
  Bmu_perp[0] = Bmu[1];
  Bmu_perp[1] = -Bmu[0];
  secondDer_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, derBmu);
  derBmu_perp[0] = derBmu[1];
  derBmu_perp[1] = -derBmu[0];
  vec2_subtraction(xMu, xB, xMuminxB);
  fEval = dotProd(xMuminxB, Bmu_perp) ; // (xMu - xB)^T(Bmu_perp) 
  derEval = dotProd(xMuminxB, derBmu_perp); // (xMu - xB)^T(derBmu_perp)
  mu_new = mu - alpha*fEval/derEval; // new step
  hermite_interpolationSpatial(mu_new, xR, xHat, BR, BHat, xMu);
  grad_hermite_interpolationSpatial(mu_new, xR, xHat, BR, BHat, Bmu);
  Bmu_perp[0] = Bmu[1];
  Bmu_perp[1] = -Bmu[0];
  secondDer_hermite_interpolationSpatial(mu_new, xR, xHat, BR, BHat, derBmu);
  derBmu_perp[0] = derBmu[1];
  derBmu_perp[1] = -derBmu[0];
  vec2_subtraction(xMu, xB, xMuminxB);
  fEval_new = dotProd(xMuminxB, Bmu_perp) ; // (xMu - xB)^T(Bmu_perp)
  while(fabs(fEval) < fabs(fEval_new)){
    alpha = alpha/2;
    mu_new = mu - alpha*fEval/derEval; // new step
    hermite_interpolationSpatial(mu_new, xR, xHat, BR, BHat, xMu);
    grad_hermite_interpolationSpatial(mu_new, xR, xHat, BR, BHat, Bmu);
    Bmu_perp[0] = Bmu[1];
    Bmu_perp[1] = -Bmu[0];
    secondDer_hermite_interpolationSpatial(mu_new, xR, xHat, BR, BHat, derBmu);
    derBmu_perp[0] = derBmu[1];
    derBmu_perp[1] = -derBmu[0];
    vec2_subtraction(xMu, xB, xMuminxB);
    fEval_new = dotProd(xMuminxB, Bmu_perp) ; // (xMu - xB)^T(Bmu_perp)
    derEval_new = dotProd(xMuminxB, derBmu_perp);
  }
  return alpha;
}

double find_minMu(double mu0, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double tol, double maxIter){
  // find the minimum mu (0<= mu <= 1) such that the segment x0xMu is parallel to Bmu)
  double mu, fEval, derEval, xMu[2], Bmu[2], Bmu_perp[2], xMuminxB[2], derBmu[2], derBmu_perp[2], change, alpha;
  int i = 1;
  mu = mu0;
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
  grad_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, Bmu);
  Bmu_perp[0] = Bmu[1];
  Bmu_perp[1] = -Bmu[0];
  secondDer_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, derBmu);
  derBmu_perp[0] = derBmu[1];
  derBmu_perp[1] = -derBmu[0];
  vec2_subtraction(xMu, xB, xMuminxB);
  fEval = dotProd(xMuminxB, Bmu_perp) ; // (xMu - xB)^T(Bmu_perp) 
  derEval = dotProd(xMuminxB, derBmu_perp); // (xMu - xB)^T(derBmu_perp)
  change = 10;
  while(fabs(fEval)> tol & i < maxIter & fabs(change)>tol){
    // Newton's method
    alpha = backTr_find_minMu(mu, 1, xA, xB, xHat, xR, BHat, BR, tol, maxIter);
    mu = mu - alpha*fEval/derEval;
    hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
    grad_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, Bmu);
    Bmu_perp[0] = Bmu[1];
    Bmu_perp[1] = -Bmu[0];
    secondDer_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, derBmu);
    derBmu_perp[0] = derBmu[1];
    derBmu_perp[1] = -derBmu[0];
    vec2_subtraction(xMu, xB, xMuminxB);
    fEval = dotProd(xMuminxB, Bmu_perp); // (xMu - xB)^T(Bmu_perp)
    derEval = dotProd(xMuminxB, derBmu_perp); // (xMu - xB)^T(derBmu_perp)
    i ++;
    if(mu<0){
      mu = 0;
    }
    if(mu>1){
      mu = 1;
    }
  }
  return mu;
}

double der_linearShootCr(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef){
  // gradient of fobjective_linearShootCr
  double lambda, xMu[2], xLam[2], Bmu[2], derBmu[2], Tprime, Bmu_halves[2], derBmu_halves[2];
  double tPrime;
  lambda = t_ofMu(mu, xA, xB, xHat, xR, BHat, BR); // because lambda is uniquely defined by mu
  linearInterpolation(lambda, xA, xB, xLam); // xLam
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu); // xMu
  Tprime = der_hermite_interpolationT(lambda, xA, xB, TA, TB, gradA, gradB); // T'(lambda)
  grad_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, Bmu); // Bmu
  secondDer_hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, derBmu); // Bmu'
  grad_hermite_interpolationSpatial((1-mu)/2, xR, xHat, BR, BHat, Bmu_halves); // Bmu of (1-mu)/2
  secondDer_hermite_interpolationSpatial((1-mu)/2, xR, xHat, BR, BHat, derBmu_halves); // Bmu' of (1-mu)/2
  tPrime = der_t_ofMu(mu, xA, xB, xHat, xR, BHat, BR); // t'(mu)
  double xMuminxLam[2], xBminxA[2];
  vec2_subtraction(xMu, xLam, xMuminxLam);
  vec2_subtraction(xB, xA, xBminxA);
  double tPrimexBminxA[2], BmumintPrimexBminxA[2];
  scalar_times_2vec(tPrime, xBminxA, tPrimexBminxA);
  vec2_subtraction(Bmu, tPrimexBminxA, BmumintPrimexBminxA);
  double der2, L;
  // we compute kind of simpson but with the (1-mu)/6 factor in front
  der2 = indexRef*dotProd(BmumintPrimexBminxA, xMuminxLam)/l2norm(xMuminxLam);
  L = (indexRef/6)*(l2norm(Bmu) + 4*l2norm(Bmu_halves) + l2norm(BHat));
  double der3;
  der3 = -L + ((indexRef*(1-mu)/6))*( dotProd(derBmu, Bmu)/l2norm(Bmu) + 2*dotProd(derBmu_halves, Bmu_halves)/l2norm(Bmu_halves));
  return Tprime*tPrime + der2 + der3;
}

double backTr_linearShootCr(double alpha0, double d, double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef){
  double f_prev, f_cur, alpha;
  int i = 0;
  alpha = alpha0;
  // evaluating the objective function
  f_prev = fobjective_linearShootCr(mu, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
  f_cur = fobjective_linearShootCr(mu-alpha*d, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
  // start backtracking if necessary
  while(f_prev<=f_cur & i<10){
    alpha = alpha*0.5;
    f_cur = fobjective_linearShootCr(mu-alpha*d, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
    i++;
  }
  return alpha;
}

double projectedGradient_linearShootCr(double mu0, double muMin, double muMax, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double tol, int maxIter, double indexRef){
  // projected gradient descent for when xHat and xR are on the boundary and we want to update xHat with xA and xB which
  // are not on the boundary but from 0 to muMin we can't join xAxB to xHat with a straight line
  // this is why we shoot to the boundary and then we do a creeping update until we reach xHat
  double der_cur, der_prev, step, alpha, mu_prev, mu_cur, test;
  int i = 1;
  mu_prev = mu0;
  der_cur = der_linearShootCr(mu0, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
  der_prev = der_cur;
  alpha = backTr_linearShootCr(1, der_cur, mu0, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
  if(fabs(der_cur)>tol){
    test = mu_prev - alpha*der_cur;
  }
  else{
    test = mu_prev;
  }
  if(test>muMax){
    test = muMax;
  }
  if(test<muMin){
    test = muMin;
  }
  mu_cur = test;
  der_cur = der_linearShootCr(mu_cur, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
  // start the iteration part
  while(i<maxIter & fabs(der_cur)>tol & fabs(mu_cur - mu_prev)>0){
    alpha = backTr_linearShootCr(1, der_cur, mu_cur, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
    test = mu_cur - alpha*der_cur;
    if(test>muMax){
      test = muMax;
    }
    if(test<muMin){
      test = muMin;
    }
    der_prev = der_cur;
    mu_prev = mu_cur;
    mu_cur = test;
    der_cur = der_linearShootCr(mu_cur, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
    double fCur;
    fCur = fobjective_linearShootCr(mu_cur, xA, xB, xHat, xR, BHat, BR, TA, TB, gradA, gradB, indexRef);
    i++;
    }
    return mu_cur;
}



double fobjective_linearShootCr(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef) {
  // objective function for when xHat and xR are on the boundary and we want to update xHat with xA and xB which
  // are not on the boundary but from muMin to 1 we can't join xAxB to xHat with a straight line
  // this is why we shoot to the boundary and then we do a creeping update until we reach xHat
  double lambda;
  lambda = t_ofMu(mu, xA, xB, xHat, xR, BHat, BR); // because lambda is uniquely defined by mu
  /* printf("\n\nWith mu: %lf\n", mu); */
  /* printf("Value of lambda: %lf\n", lambda); */
  double xLam[2], xMu[2], L, Tlam;
  linearInterpolation(lambda, xA, xB, xLam);
  hermite_interpolationSpatial(mu, xR, xHat, BR, BHat, xMu);
  Tlam = hermite_interpolationT(lambda, xA, xB, TA, TB, gradA, gradB);
  double xMuminxLam[2];
  vec2_subtraction(xMu, xLam, xMuminxLam);
  double arcL;
  arcL = arclength_hermiteSimpson(mu, 1, xR, xHat, BR, BHat);
  //printf("The computed L: %lf\n", arcL);
  return Tlam + indexRef*l2norm(xMuminxLam) + indexRef*arcL;
}

double tN_ofLamMu(double lambda, double mu, double xA[2], double xB[2], double BA[2], double BB[2], double xHat[2], double BHat[2]){
  // t for a nonlinear shoot and creep update depends on BOTH lambda and mu
  double xLam[2], Bmu[2], Bmu_perp[2];
  hermite_interpolationSpatial(lambda, xA, xB, BA, BB, xLam); // xLam
  grad_hermite_interpolationSpatial(mu, xA, xHat, BA, BHat, Bmu); // Bmu
  Bmu_perp[0] = Bmu[1];
  Bmu_perp[1] = -Bmu[0];
  // then compute the subtraction and the dot product
  double xHatminxLam[2];
  vec2_subtraction(xHat, xLam, xHatminxLam);
  return dotProd(xHatminxLam, Bmu_perp);
}

double find_MufromLam_tN(double mu0, double mu1, double lambda, double xA[2], double xB[2], double BA[2], double BB[2], double xHat[2], double BHat[2], double tol, int maxIter) {
  // find mu from lambda for tN_ofLamMu = 0 (i.e. <xHat - xLam, Bmu_perp> = 0)
  // uses the secant method
  double tN_cur, tN_prev, mu_cur, mu_prev, aux_mu;
  int i = 1;
  mu_cur = mu1;
  mu_prev = mu0;
  tN_cur = tN_ofLamMu(lambda, mu_cur, xA, xB, BA, BB, xHat, BHat);
  tN_prev = tN_ofLamMu(lambda, mu_prev, xA, xB, BA, BB, xHat, BHat);
  while( fabs(tN_cur) > tol & i < maxIter){
    aux_mu = mu_cur;
    mu_cur = (mu_prev*tN_cur - mu_cur*tN_prev)/(tN_cur - tN_prev);
    mu_prev = aux_mu;
    tN_cur = tN_ofLamMu(lambda, mu_cur, xA, xB, BA, BB, xHat, BHat);
    tN_prev = tN_ofLamMu(lambda, mu_prev, xA, xB, BA, BB, xHat, BHat);
    i++;
  }
  return mu_cur;
}

double find_LamfromMu_tN(double lam0, double lam1, double mu, double xA[2], double xB[2], double BA[2], double BB[2], double xHat[2], double BHat[2], double tol, int maxIter) {
  // find lambda from mu for tN_ofLamMu = 0 (i.e. <xHat - xLam, Bmu_perp> = 0)
  // uses the secant method
  double tN_cur, tN_prev, lam_cur, lam_prev, aux_lam;
  int i = 1;
  lam_cur = lam1;
  lam_prev = lam0;
  tN_cur = tN_ofLamMu(lam_cur, mu, xA, xB, BA, BB, xHat, BHat);
  tN_prev = tN_ofLamMu(lam_prev, mu, xA, xB, BA, BB, xHat, BHat);
  while( fabs(tN_cur) > tol & i < maxIter) {
    aux_lam = lam_cur;
    lam_cur = (lam_prev*tN_cur - lam_cur*tN_prev)/(tN_cur - tN_prev);
    lam_prev = aux_lam;
    tN_cur = tN_ofLamMu(lam_cur, mu, xA, xB, BA, BB, xHat, BHat);
    tN_prev = tN_ofLamMu(lam_prev, mu, xA, xB, BA, BB, xHat, BHat);
    i++;
  }
  return lam_cur;
}


void grad_NonLinShootCr(double lambda, double mu, double xA[2], double xB[2], double BA[2], double BB[2], double xHat[2], double BHat[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef, double grad[2]) {
  double Tprime, xMu[2], xLam[2], Bmu[2], derBmu[2], Blam[2], Bmu_halves[2], derBmu_halves[2];
  hermite_interpolationSpatial(lambda, xA, xB, BA, BB, xLam); // xLam
  hermite_interpolationSpatial(mu, xA, xHat, BA, BHat, xMu); // xMu
  Tprime = der_hermite_interpolationT(lambda, xA, xB, TA, TB, gradA, gradB); // T'(lambda)
  grad_hermite_interpolationSpatial(lambda, xA, xB, BA, BB, Blam); // Blam
  grad_hermite_interpolationSpatial(mu, xA, xHat, BA, BHat, Bmu); // Bmu
  grad_hermite_interpolationSpatial((1+mu)/2, xA, xHat, BA, BHat, Bmu_halves); // Bmu
  secondDer_hermite_interpolationSpatial(mu, xA, xHat, BA, BHat, derBmu); // Bmu'
  secondDer_hermite_interpolationSpatial((1+mu)/2, xA, xHat, BA, BHat, derBmu_halves); // Bmu_halves'
  double xMuminxLam[2], normxMuminxLam;
  vec2_subtraction(xMu, xLam, xMuminxLam); // xMu - xLam
  normxMuminxLam = l2norm(xMuminxLam);
  double normBmu, normBmu_halves, normBHat;
  normBmu = l2norm(Bmu);
  normBmu_halves = l2norm(Bmu_halves);
  normBHat = l2norm(BHat);
  // put everything together
  grad[0] = Tprime - indexRef*(dotProd(Blam, xMuminxLam)/normxMuminxLam );
  double grad1, grad2;
  grad1 = (indexRef/6)*( normBmu + 4*normBmu_halves + normBHat);
  grad2 = indexRef*((1-mu)/6)*( dotProd(derBmu, Bmu)/normBmu + 2*dotProd(derBmu_halves, Bmu_halves)/normBmu_halves);
  grad[1] = indexRef*(dotProd(Bmu, xMuminxLam)/normxMuminxLam) - grad1 + grad2;
}

double backTr_NonLinShootCr(double alpha0, double d, double lambda, double mu, double xA[2], double xB[2], double BA[2], double BB[2], double xHat[2], double BHat[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef) {
  // backtracking algorithm for the nonlinear shoot and creep optimization regine
}

double fobjective_NonLinShootCr(double lambda, double mu, double xA[2], double xB[2], double BA[2], double BB[2], double xHat[2], double BHat[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef) {
  // THat for a non linear shoot and creep update
  double xLam[2], xMu[2], Tlam, xMuminxLam[2];
  hermite_interpolationSpatial(lambda, xA, xB, BA, BB, xLam); // xLam
  hermite_interpolationSpatial(mu, xA, xHat, BA, BHat, xMu); // xMu
  vec2_subtraction(xMu, xLam, xMuminxLam); // xMu - xLam
  Tlam = hermite_interpolationT(lambda, xA, xB, TA, TB, gradA, gradB); // T(lambda)
  double arcL;
  arcL = arclength_hermiteSimpson(mu, 1, xA, xHat, BA, BHat);
  // put everything together
  return Tlam + indexRef*l2norm(xMuminxLam) + indexRef*arcL;
  
}



