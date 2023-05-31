
/*
Optimization methods for the 2D JMM
*/

#include "opti_method.h"
#include "linAlg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>




double eikApproxCubic(double T0, double T1, double grad0[2], double grad1[2],
		    double lambda, double x0[2], double x1[2], double xHat[2], double indexRef) {
    // this is the approximation to the eikonal solution using Hermite interpolation (because
   // we also approximate the gradient so we use it

    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xLamMinxHat[2], normAux, min1Lambda, gApprox;
    double lambdaSq, lambdaCub, hermPol, dotXGrad0, dotXGrad1, x1Minx0[2];
      
    xlambda[0] = 0;
    xlambda[1] = 0;
    firstPartxLambda[0] = 0;
    firstPartxLambda[1] = 0;
    secondPartxLambda[0] = 0;
    secondPartxLambda[1] = 0;
    xLamMinxHat[0] = 0;
    xLamMinxHat[1] = 0;
    min1Lambda = 1-lambda;
    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1
    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1
    vec2_subtraction( xlambda, xHat, xLamMinxHat ); // xlambda - xHat
    normAux = l2norm( xLamMinxHat ); //   || xlambda - xHat ||
    // For the part of the Hermite interpolation we need the following auxiliary variables
    lambdaSq = pow(lambda, 2); // lambda^2
    lambdaCub = pow(lambda, 3); // lambda^3
    vec2_subtraction(x1, x0, x1Minx0); // x1-x0
    dotXGrad0 = dotProd(x1Minx0, grad0); // (x1-x0)T(grad0)
    dotXGrad1 = dotProd(x1Minx0, grad1); // (x1-x0)T(grad1)
    hermPol = T0*(2*lambdaCub - 3*lambdaSq + 1) + T1*(-2*lambdaCub + 3*lambdaSq);
    hermPol += dotXGrad0*(lambdaCub - 2*lambdaSq + lambda) + dotXGrad1*(lambdaCub - lambdaSq);
    // now we have the approximation gApprox = HermiteInterpolation + nu*distance
    gApprox = hermPol + indexRef*normAux;

    return gApprox;
}

double gPrimeCubic(double T0, double T1, double grad0[2], double grad1[2], double lambda,
		   double x0[2], double x1[2], double xHat[2], double indexRef){
    // auxiliary function to compute the function gPrime when we are approximation T(xLambda)
   // using Hermite interpolation
    double x1Minx0[2], dotProduct, normAux, gPrim;
    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xLamMinxHat[2], min1Lambda;
    double dotXGrad0, dotXGrad1, lambdaSq, devHerm;

    xlambda[0] = 0;
    xlambda[1] = 0;
    firstPartxLambda[0] = 0;
    firstPartxLambda[1] = 0;
    secondPartxLambda[0] = 0;
    secondPartxLambda[1] = 0;
    xLamMinxHat[0] = 0;
    xLamMinxHat[1] = 0;
    min1Lambda = 1-lambda;

    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1
    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1
    vec2_subtraction( xlambda, xHat, xLamMinxHat ); // xlambda - xHat
    normAux = l2norm( xLamMinxHat ); //   || xlambda - xHat ||

    x1Minx0[0] = 0;
    x1Minx0[1] = 0;

    dotProduct = 0;
    normAux= 0;

    vec2_subtraction(x1, x0, x1Minx0);

    dotProduct = dotProd(x1Minx0, xLamMinxHat);
    normAux = l2norm(xLamMinxHat);

    // For the part of the derivative that corresponds to the hermite polynomial:
    lambdaSq = pow(lambda, 2); // lambda^2
    vec2_subtraction(x1, x0, x1Minx0); // x1-x0
    dotXGrad0 = dotProd(x1Minx0, grad0); // (x1-x0)T(grad0)
    dotXGrad1 = dotProd(x1Minx0, grad1); // (x1-x0)T(grad1)
    devHerm = T0*(6*lambdaSq -6*lambda) + T1*(-6*lambdaSq + 6*lambda);
    devHerm += dotXGrad0*(3*lambdaSq - 4*lambda + 1) + dotXGrad1*(3*lambdaSq - 2*lambda);

    gPrim = devHerm + indexRef*dotProduct/normAux;

    return gPrim;

}

double secantCubic_2D(double lambda0, double lambda1, double T0, double T1, double grad0[2], double grad1[2],
		     double x0[2], double x1[2], double xHat[2], double tol, int maxIter, double indexRef){
    // This method is the implementation of the secant method to optimize the path from xlambda to xhat
   // T(xlambda) is approximated with a cubic Hermite polynomial 
    int k = 1;
    double lam, gPrime0, gPrime1;

    gPrime1 = gPrimeCubic(T0, T1, grad0, grad1, lambda1, x0, x1, xHat, indexRef);
    // printf("Initial, with lambda = %fl   the value of gPrime is: %fl\n", lambda0, gPrime(T0, T1, lambda0, x0, x1, xHat, indexRef));
    // printf("Initial, with lambda = %fl   the value of gPrime is: %fl\n", lambda1, gPrime1);
    
    while(k < maxIter & fabs(gPrime1)>tol){
      gPrime0 = gPrimeCubic(T0, T1, grad0, grad1, lambda0, x0, x1, xHat, indexRef);
      gPrime1 = gPrimeCubic(T0, T1, grad0, grad1, lambda1, x0, x1, xHat, indexRef);
      lam = lambda1 - gPrime1*(lambda1 - lambda0)/( gPrime1 - gPrime0 );
       lambda0 = lambda1;
       lambda1 = lam;
       k ++;
    }
    if(lambda1<0){
        lambda1 = 0;
    }
    if(lambda1 > 1){
        lambda1 = 1;
    }
    return lambda1;
}
