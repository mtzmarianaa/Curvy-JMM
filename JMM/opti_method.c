

/*
Optimization methods for the 2D FMM
*/

#include "opti_method.h"
#include "linAlg.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double eikApproxLin(double T0, double T1, double lambda, double x0[2], double x1[2], double xHat[2], double indexRef) {
    // this is the approximation to the eikonal solution using local approximations everywhere

    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xLamMinxHat[2], normAux, min1Lambda, gApprox;

    xlambda[0] = 0;
    xlambda[1] = 0;
    firstPartxLambda[0] = 0;
    firstPartxLambda[1] = 0;
    secondPartxLambda[0] = 0;
    secondPartxLambda[1] = 0;
    xLamMinxHat[0] = 0;
    xLamMinxHat[1] = 0;
    min1Lambda = 1-lambda;
    // printf("   aux printing eikApproxLin\n");
    // printf("            1-lambda = %fl\n", min1Lambda);
    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    // printf("            (1-lambda)*x0 = (%fl,%fl)\n", firstPartxLambda[0], firstPartxLambda[1]);
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1
    // printf("            lambda*x1 = (%fl,%fl)\n", secondPartxLambda[0], secondPartxLambda[1]);
    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1
    // printf("            xlambda = (%fl,%fl)\n", xlambda[0], xlambda[1]);
    vec2_subtraction( xlambda, xHat, xLamMinxHat ); // xlambda - xHat
    // printf("            xlambda-xHat = (%fl, %fl)\n", xLamMinxHat[0], xLamMinxHat[1]);
    normAux = l2norm( xLamMinxHat ); //   || xlambda - xHat ||
    // printf("            norm(xlambda-xHat) = %fl\n", normAux);

    gApprox = min1Lambda*T0 + lambda*T1 + indexRef*normAux;

    return gApprox;
}

double gPrime(double T0, double T1, double lambda, double x0[2], double x1[2], double xHat[2], double indexRef){
    // auxiliary function to compute the function gPrime
    double x1Minx0[2], dotProduct, normAux, gPrim;
    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xLamMinxHat[2], min1Lambda;

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

    gPrim = T1 - T0 + indexRef*dotProduct/normAux;

    return gPrim;

}

double secant_2D(double lambda0, double lambda1, double T0, double T1, double x0[2], double x1[2], double xHat[2], double tol, int maxIter, double indexRef){
    // This method is the implementation of the secant method for the 2D fmm using the
    // function defined in SoSFunction.h as the speed of sound
    int k = 1;
    double lam, gPrime0, gPrime1;

    gPrime1 = gPrime(T0, T1, lambda1, x0, x1, xHat, indexRef);
    // printf("Initial, with lambda = %fl   the value of gPrime is: %fl\n", lambda0, gPrime(T0, T1, lambda0, x0, x1, xHat, indexRef));
    // printf("Initial, with lambda = %fl   the value of gPrime is: %fl\n", lambda1, gPrime1);
    
    while(k < maxIter & fabs(gPrime1)>tol){
        // u
        gPrime0 = gPrime(T0, T1, lambda0, x0, x1, xHat, indexRef);
        gPrime1 = gPrime(T0, T1, lambda1, x0, x1, xHat, indexRef);
        lam = lambda1 - gPrime1*(lambda1 - lambda0)/( gPrime1 - gPrime0 );
        // printf("In this iteration lambda_i has value of %fl", lam);
        // printf("    |  value of T(xhat): %fl   with %fl as index of regraction\n", eikApproxLin(T0, T1, lam, x0, x1, xHat, indexRef), indexRef );
        // printf("    |   In this iteration the value of gPrime is: %fl \n", gPrime(T0, T1, lam, x0, x1, xHat, indexRef));
        lambda0 = lambda1;
        lambda1 = lam;
        k ++;
    }
    // check that lambda1 stays in the interval [0, 1]
    if(lambda1<0){
        lambda1 = 0;
    }
    if(lambda1 > 1){
        lambda1 = 1;
    }
    // printf("Optimum lambda found %fl\n", lambda1);
    // printf("With gPrime value of: %fl\n", gPrime(T0, T1, lambda1, x0, x1, xHat, indexRef));
    // printf("With the secant method, the optimal lambda found (for a simple update) is: %fl\n", lambda1);
    return lambda1;
}


// for the "two step" optimization part (i.e. the piece wise linear path goes through two different regions with different indices of refraction)

double eikApproxLin_2Regions(double T0, double T1, double lambda, double mu, double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02) {
    // piece wise linear approximation of a ray from x0x1 passing through a change of region defined at x0x2 to get to xHat

    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xMuMinxLam[2], normAux1, min1Lambda, gApprox;
    double xmu[2], firstPartxMu[2], seconPartxMu[2], xHatMinxMu[2], normAux2, min1Mu;

    min1Lambda = 1 - lambda;
    min1Mu = 1 - mu;

    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1

    scalar_times_2vec( min1Mu, x0, firstPartxMu ); // (1-mu)*x0
    scalar_times_2vec(mu, x2, seconPartxMu); // mu*x2

    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1

    vec2_addition( firstPartxMu, seconPartxMu, xmu ); // (1-mu)*x0 + mu*x2

    vec2_subtraction( xmu, xlambda, xMuMinxLam ); // xmu - xlambda

    vec2_subtraction( xHat, xmu, xHatMinxMu ); // xhat - xmu

    normAux1 = l2norm( xMuMinxLam ); //   || xmu - xlambda ||

    normAux2 = l2norm( xHatMinxMu ); //   || xhat - xmu ||

    gApprox = min1Lambda*T0 + lambda*T1 + indexRef_01*normAux1 + indexRef_02*normAux2;

    return gApprox;
}

void gradient_2Regions(double grad[2], double T0, double T1, double lambda, double mu, double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02) {
    // auxiliary function to calculate the gradient of eikApproxLin_2Regions, the first entry corresponds to the partial with respect to lambda
    // and the second entry to the partial with respect to mu
    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xMuMinxLam[2], normAux1, min1Lambda;
    double xmu[2], firstPartxMu[2], seconPartxMu[2], xHatMinxMu[2], normAux2, min1Mu;
    double x1Minx0[2], x2Minx0[2], dotProd1, dotProd2, dotProd3;

    min1Lambda = 1 - lambda;
    min1Mu = 1 - mu;

    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1

    scalar_times_2vec( min1Mu, x0, firstPartxMu ); // (1-mu)*x0
    scalar_times_2vec(mu, x2, seconPartxMu); // mu*x2

    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1

    vec2_addition( firstPartxMu, seconPartxMu, xmu ); // (1-mu)*x0 + mu*x2

    vec2_subtraction( xmu, xlambda, xMuMinxLam ); // xmu - xlambda

    vec2_subtraction( xHat, xmu, xHatMinxMu ); // xhat - xmu

    normAux1 = l2norm( xMuMinxLam ); //   || xmu - xlambda ||

    normAux2 = l2norm( xHatMinxMu ); //   || xhat - xmu ||

    vec2_subtraction(x1, x0, x1Minx0); // x1 - x0

    vec2_subtraction(x2, x0, x2Minx0); // x2 - x0

    dotProd1 = dotProd( xMuMinxLam, x1Minx0 ); // (xmu-xlam)T(x1-x0)
    dotProd2 = dotProd( xMuMinxLam, x2Minx0 ); // (xmu-xlam)T(x2-x0)
    dotProd3 = dotProd(xHatMinxMu, x2Minx0); // (xhat - xmu)T(x2-x0)
    // assemble
    grad[0] = T1 - T0 - indexRef_01*dotProd1/normAux1;
    grad[1] = indexRef_01*dotProd2/normAux1 - indexRef_02*dotProd3/normAux2;
}

void projectedGradientDescent(double optimizers[2], double T0, double T1, double x0[2], double x1[2], double x2[2], double xHat[2], double tol, int maxIter, double indexRef_01, double indexRef_02) {
    // this is the optimization method to find lambda and mu such that xlam is on the segment x0x1, xmu is on the segment x0x2, where
    // x0x2 defined the change in region from index of refraction indexRef_01 to index of refraction indexRef_02
    // optimizers[0] = lambda, optimizers[1] = mu
    // this is the projected gradient descent method with step sizes equal to 0.01
    double grad[2], step[2], yk[2], projectedGradient[2], t;
    int iter;
    iter = 0;
    optimizers[0] = 1.0;
    optimizers[1] = 1.0; // initial value of the optimizers
    // printf("\n\nIteration %d\n", iter);
    // printf("Optimizers at this step: %fl    , %fl\n", optimizers[0], optimizers[1]);
    // we have to calculate the current gradient
    gradient_2Regions(grad, T0, T1, optimizers[0], optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
    // printf("Gradient: (%fl, %fl)\n", grad[0], grad[1]);
    // steepest descent and calculation of the projected gradient
    projection01Cube(grad, projectedGradient);
    // printf("Projected gradient: (%fl, %fl)\n", projectedGradient[0], projectedGradient[1]);

    // iteration part
    while( l2norm(projectedGradient) > tol & iter < maxIter  ){
        iter ++;
        // printf("\n\nIteration %d\n", iter);
        // we have to calculate the current gradient
        // printf("Projected gradient: (%fl, %fl)  with norm %fl\n", projectedGradient[0], projectedGradient[1], l2norm(projectedGradient));
        t = backtracking(optimizers, grad, T0, T1, x0, x1, x2, xHat, indexRef_01, indexRef_02);
        scalar_times_2vec(t, grad, step);
        // printf("This steps %fl with direction: (%fl, %fl)\n", t, step[0], step[1]);
        vec2_subtraction( optimizers, step, yk );
        // printf("Without projecting back to the domain, the next optimizers would be   %fl   %fl\n", yk[0], yk[1]);
        // then we need to project it back to the feasible set
        projection01Cube(yk, optimizers);
        // printf("Optimizers at this step: %fl    , %fl\n", optimizers[0], optimizers[1]);
        // Tcurrent = eikApproxLin_2Regions(T0, T1, optimizers[0], optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
        // printf("Current value of T: %fl\n", Tcurrent);
        gradient_2Regions(grad, T0, T1, optimizers[0], optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
        // printf("Gradient: (%fl, %fl)  with norm %fl\n", grad[0], grad[1], l2norm(grad));
        // steepest descent and calculation of the projected gradient
        projection01Cube(grad, projectedGradient);
        if(optimizers[0] == 0 & optimizers[1] == 0){
            break;
        }
    }

    // printf("With the projected gradient method the value of lambda: %fl    , the value of mu: %fl\n", optimizers[0], optimizers[1]);

}

double backtracking(double optimizers[2], double gradient[2], double T0, double T1, double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02){
    // backtracking algorithm to find the step size of the projected gradient descent method so that the
    // Armijo Wolfe conditions are met
    double t, alpha, beta, Tval_t, Tval, dotGradDirection, grad[2], direction[2];
    direction[0] = -gradient[0];
    direction[1] = -gradient[1];
    t = 0.05;
    alpha = 0.05;
    beta = 0.1;
    Tval_t = eikApproxLin_2Regions(T0, T1, optimizers[0] + t*direction[0], optimizers[1] + t*direction[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
    Tval = eikApproxLin_2Regions(T0, T1, optimizers[0], optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
    gradient_2Regions(grad, T0, T1, optimizers[0], optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
    dotGradDirection = dotProd( grad, direction );
    // backtracking part
    while(  Tval_t > Tval + alpha*t*dotGradDirection  ){
        t = beta*t;
        Tval_t = eikApproxLin_2Regions(T0, T1, optimizers[0] + t*direction[0], optimizers[1] + t*direction[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
    }
    return t;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
///// OPTIMIZATION METHODS - INTERPOLATING T(LAMBDA) WITH CUBIC HERMITE INTERPOLATION

double eikApprox_freeSpace(double T0, double T1, double grad0[2], double grad1[2],
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

double secant_freeSpace(double lambda0, double lambda1, double T0, double T1, double grad0[2], double grad1[2],
		     double x0[2], double x1[2], double xHat[2], double tol, int maxIter, double indexRef){
  // optimization on free space with Hermite interpolation for T
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

double eikApproxCubic_2Regions(double T0, double T1, double grad0[2], double grad1[2],
			       double lambda, double mu, double x0[2], double x1[2],
			       double x2[2], double xHat[2], double indexRef_01, double indexRef_02) {
    // approximation of a ray from x0x1 passing through a change of region defined at x0x2 to get to xHat
    // the approximation of T(xlambda) is computed using cubic hermite polynomial

    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xMuMinxLam[2], normAux1, min1Lambda, gApprox;
    double xmu[2], firstPartxMu[2], seconPartxMu[2], xHatMinxMu[2], normAux2, min1Mu;
    double lambdaSq, lambdaCub, hermPol, dotXGrad0, dotXGrad1, x1Minx0[2];

    min1Lambda = 1 - lambda;
    min1Mu = 1 - mu;

    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1

    scalar_times_2vec( min1Mu, x0, firstPartxMu ); // (1-mu)*x0
    scalar_times_2vec(mu, x2, seconPartxMu); // mu*x2

    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1

    vec2_addition( firstPartxMu, seconPartxMu, xmu ); // (1-mu)*x0 + mu*x2

    vec2_subtraction( xmu, xlambda, xMuMinxLam ); // xmu - xlambda

    vec2_subtraction( xHat, xmu, xHatMinxMu ); // xhat - xmu

    normAux1 = l2norm( xMuMinxLam ); //   || xmu - xlambda ||

    normAux2 = l2norm( xHatMinxMu ); //   || xhat - xmu ||

    // For the part of the Hermite interpolation we need the following auxiliary variables
    lambdaSq = pow(lambda, 2); // lambda^2
    lambdaCub = pow(lambda, 3); // lambda^3
    vec2_subtraction(x1, x0, x1Minx0); // x1-x0
    dotXGrad0 = dotProd(x1Minx0, grad0); // (x1-x0)T(grad0)
    dotXGrad1 = dotProd(x1Minx0, grad1); // (x1-x0)T(grad1)
    hermPol = T0*(2*lambdaCub - 3*lambdaSq + 1) + T1*(-2*lambdaCub + 3*lambdaSq);
    hermPol += dotXGrad0*(lambdaCub - 2*lambdaSq + lambda) + dotXGrad1*(lambdaCub - lambdaSq);

    gApprox = hermPol + indexRef_01*normAux1 + indexRef_02*normAux2;

    return gApprox;
}

void gradientCubic_2Regions(double grad[2], double T0, double T1, double grad0[2], double grad1[2],
			    double lambda, double mu, double x0[2], double x1[2], double x2[2],
			    double xHat[2], double indexRef_01, double indexRef_02) {
    // auxiliary function to calculate the gradient of eikApproxCubic_2Regions,
    //the first entry corresponds to the partial with respect to lambda
    // and the second entry to the partial with respect to mu
    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xMuMinxLam[2], normAux1, min1Lambda;
    double xmu[2], firstPartxMu[2], seconPartxMu[2], xHatMinxMu[2], normAux2, min1Mu;
    double x1Minx0[2], x2Minx0[2], dotProd1, dotProd2, dotProd3;
    double dotXGrad0, dotXGrad1, lambdaSq, devHerm;

    min1Lambda = 1 - lambda;
    min1Mu = 1 - mu;

    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1

    scalar_times_2vec( min1Mu, x0, firstPartxMu ); // (1-mu)*x0
    scalar_times_2vec(mu, x2, seconPartxMu); // mu*x2

    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1

    vec2_addition( firstPartxMu, seconPartxMu, xmu ); // (1-mu)*x0 + mu*x2

    vec2_subtraction( xmu, xlambda, xMuMinxLam ); // xmu - xlambda

    vec2_subtraction( xHat, xmu, xHatMinxMu ); // xhat - xmu

    normAux1 = l2norm( xMuMinxLam ); //   || xmu - xlambda ||

    normAux2 = l2norm( xHatMinxMu ); //   || xhat - xmu ||

    vec2_subtraction(x1, x0, x1Minx0); // x1 - x0

    vec2_subtraction(x2, x0, x2Minx0); // x2 - x0

    dotProd1 = dotProd( xMuMinxLam, x1Minx0 ); // (xmu-xlam)T(x1-x0)
    dotProd2 = dotProd( xMuMinxLam, x2Minx0 ); // (xmu-xlam)T(x2-x0)
    dotProd3 = dotProd(xHatMinxMu, x2Minx0); // (xhat - xmu)T(x2-x0)

    // For the part of the derivative that corresponds to the hermite polynomial in grad[0]:
    lambdaSq = pow(lambda, 2); // lambda^2
    vec2_subtraction(x1, x0, x1Minx0); // x1-x0
    dotXGrad0 = dotProd(x1Minx0, grad0); // (x1-x0)T(grad0)
    dotXGrad1 = dotProd(x1Minx0, grad1); // (x1-x0)T(grad1)
    devHerm = T0*(6*lambdaSq -6*lambda) + T1*(-6*lambdaSq + 6*lambda);
    devHerm += dotXGrad0*(3*lambdaSq - 4*lambda + 1) + dotXGrad1*(3*lambdaSq - 2*lambda);
    
    // assemble
    grad[0] = devHerm - indexRef_01*dotProd1/normAux1;
    grad[1] = indexRef_01*dotProd2/normAux1 - indexRef_02*dotProd3/normAux2;
}

void projectedGradientDescentCubic(double optimizers[2], double T0, double T1, double grad0[2], double grad1[2],
			      double x0[2], double x1[2], double x2[2], double xHat[2], double tol,
			      int maxIter, double indexRef_01, double indexRef_02) {
    // this is the optimization method to find lambda and mu such that xlam is on the segment x0x1,
    //xmu is on the segment x0x2, where
    // x0x2 defined the change in region from index of refraction indexRef_01 to index of refraction indexRef_02
    // optimizers[0] = lambda, optimizers[1] = mu
    // this is the projected gradient descent method with step sizes equal to 0.01
    double grad[2], step[2], yk[2], projectedGradient[2], t;
    int iter;
    iter = 0;
    optimizers[0] = 1.0;
    optimizers[1] = 1.0; // initial value of the optimizers
    // we have to calculate the current gradient
    gradientCubic_2Regions(grad, T0, T1, grad0,grad1,optimizers[0], optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
    // steepest descent and calculation of the projected gradient
    projection01Cube(grad, projectedGradient);

    // iteration part
    while( l2norm(projectedGradient) > tol & iter < maxIter  ){
        iter ++;
        // we have to calculate the current gradient
        t = backtracking(optimizers, grad, T0, T1, x0, x1, x2, xHat, indexRef_01, indexRef_02);
        scalar_times_2vec(t, grad, step);
        vec2_subtraction( optimizers, step, yk );
        // then we need to project it back to the feasible set
        projection01Cube(yk, optimizers);
        gradientCubic_2Regions(grad,T0,T1,grad0,grad1,optimizers[0],optimizers[1], x0, x1, x2, xHat, indexRef_01, indexRef_02);
        // steepest descent and calculation of the projected gradient
        projection01Cube(grad, projectedGradient);
        if(optimizers[0] == 0 & optimizers[1] == 0){
            break;
        }
    }
}

double der_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef) {
  // derivative with respect of lambda from the function to minimize for the update in which the
  // segment x0x1 is on the boundary but xHat if fully contained in a region with index indexRef
  double lambda3, lambda2;
  lambda2 = lambda*lambda;
  lambda3 = lambda2*lambda;
  double x1Minx0[2], B0PlusB1[2], B1Plus2B0[2], twoB0[2], grad0Plusgrad1[2], twoGrad0PlusGrad1[2], twoGrad0[2];
  double lamB0[2], x0MinxHat[2];
  // first time we gather terms
  vec2_subtraction(x1, x0, x1Minx0);
  vec2_addition(B0, B1, B0PlusB1);
  scalar_times_2vec(2, B0, twoB0);
  vec2_addition(B1, twoB0, B1Plus2B0);
  vec2_addition(grad0, grad1, grad0Plusgrad1);
  scalar_times_2vec(2, grad0, twoGrad0);
  vec2_addition(twoGrad0, grad1, twoGrad0PlusGrad1);
  scalar_times_2vec(lambda, B0, lamB0);
  vec2_subtraction(x0, xHat, x0MinxHat);
  // second time gathering terms
  double dotProd1, dotProd2, sixX0minx1[2], threeB0plusB1[2], twoB1Plus2B0[2], twox0Minx1[2], threex0Minx1[2];
  dotProd1 = dotProd(x1Minx0, grad0Plusgrad1);
  dotProd2 = dotProd(x1Minx0, twoGrad0PlusGrad1);
  scalar_times_2vec(-6, x1Minx0, sixX0minx1);
  scalar_times_2vec(3, B0PlusB1, threeB0plusB1);
  scalar_times_2vec(2, B1Plus2B0, twoB1Plus2B0);
  scalar_times_2vec(-2, x1Minx0, twox0Minx1);
  scalar_times_2vec(-3, x1Minx0, threex0Minx1);
  // third time gathering terms
  double coeflam2_1[2], coeflam_1[2], coef_lam3_2[2], coef_lam2_2[2], coef_lam_2[2];
  vec2_addition(sixX0minx1, threeB0plusB1, coeflam2_1);
  vec2_addition(sixX0minx1,  twoB1Plus2B0, coeflam_1);
  vec2_addition(twox0Minx1, B0PlusB1, coef_lam3_2);
  vec2_addition(threex0Minx1, B1Plus2B0, coef_lam2_2);
  scalar_times_2vec(lambda, B0, coef_lam_2);
  // fourth time gathering terms
  double lam2_1[2], lam1_1[2], rest1_1[2], lam3_2[2], lam2_2[2], rest1_2[2], rest2_2[2];
  scalar_times_2vec(lambda2, coeflam2_1, lam2_1);
  scalar_times_2vec(lambda, coeflam_1, lam1_1);
  vec2_subtraction(lam2_1, lam1_1, rest1_1);
  scalar_times_2vec(lambda3, coef_lam3_2, lam3_2);
  scalar_times_2vec(lambda2, coef_lam2_2, lam2_2);
  vec2_addition(lamB0, x0MinxHat, rest2_2);
  vec2_subtraction(lam3_2, lam2_2, rest1_2);
  // fifth time gathering terms
  double derInside[2], xHatMinxLam[2], dotProd3, norm1, boundaryPart, tLamPart, dotProd4;
  vec2_addition(rest1_1, B0, derInside);
  vec2_addition(rest1_2, rest2_2, xHatMinxLam);
  dotProd3 = dotProd(derInside, xHatMinxLam);
  dotProd4 = dotProd(x1Minx0, grad0);
  norm1 = l2norm(xHatMinxLam);
  boundaryPart = indexRef*dotProd3/norm1;
  tLamPart = 6*T0*lambda2 - 6*T1*lambda2 + 3*lambda2*dotProd1 -6*T0*lambda + 6*T1*lambda - 2*lambda*dotProd2 + dotProd4;
  // Putting everything together
  return tLamPart + boundaryPart;
}

double backTr_fromEdge(double alpha0, double d, double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef){
  // backtracking method for projected gradient descent from updates from the edge of the domain
  double f_prev, f_cur, alpha;
  alpha = alpha0;
  // EVALUATING THE OBJECTIVE FUNCTION
  f_prev = fobjective_fromEdge(lambda, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  f_cur = fobjective_fromEdge(lambda - alpha*d, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  while(f_prev < f_cur){
    alpha = alpha*0.5;
    f_cur = fobjective_fromEdge(lambda - alpha*d, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  }
  printf("Objective function currently %lf\n", f_cur);
  return alpha;
}

double fobjective_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef){
  double lambda3, lambda2;
  lambda2 = lambda*lambda;
  lambda3 = lambda2*lambda;
  double x1Minx0[2], B0PlusB1[2], B1Plus2B0[2], twoB0[2], grad0Plusgrad1[2], twoGrad0PlusGrad1[2], twoGrad0[2];
  double lamB0[2], x0MinxHat[2];
  // first time we gather terms
  vec2_subtraction(x1, x0, x1Minx0);
  vec2_addition(B0, B1, B0PlusB1);
  scalar_times_2vec(2, B0, twoB0);
  vec2_addition(B1, twoB0, B1Plus2B0);
  vec2_addition(grad0, grad1, grad0Plusgrad1);
  scalar_times_2vec(2, grad0, twoGrad0);
  vec2_addition(twoGrad0, grad1, twoGrad0PlusGrad1);
  scalar_times_2vec(lambda, B0, lamB0);
  vec2_subtraction(x0, xHat, x0MinxHat);
  // second time gathering terms
  double dotProd1, dotProd2, twox0Minx1[2], threex0Minx1[2];
  dotProd1 = dotProd(x1Minx0, grad0Plusgrad1);
  dotProd2 = dotProd(x1Minx0, twoGrad0PlusGrad1);
  scalar_times_2vec(-2, x1Minx0, twox0Minx1);
  scalar_times_2vec(-3, x1Minx0, threex0Minx1);
  // third time gathering terms
  double coef_lam3_2[2], coef_lam2_2[2], coef_lam_2[2];
  vec2_addition(twox0Minx1, B0PlusB1, coef_lam3_2);
  vec2_addition(threex0Minx1, B1Plus2B0, coef_lam2_2);
  scalar_times_2vec(lambda, B0, coef_lam_2);
  // fourth time gathering terms
  double lam3_2[2], lam2_2[2], rest1_2[2], rest2_2[2];
  scalar_times_2vec(lambda3, coef_lam3_2, lam3_2);
  scalar_times_2vec(lambda2, coef_lam2_2, lam2_2);
  vec2_addition(lamB0, x0MinxHat, rest2_2);
  vec2_subtraction(lam3_2, lam2_2, rest1_2);
  // fifth time gathering terms
  double xHatMinxLam[2], norm1, boundaryPart, tLamPart, dotProd4;
  vec2_addition(rest1_2, rest2_2, xHatMinxLam);
  dotProd4 = dotProd(x1Minx0, grad0);
  norm1 = l2norm(xHatMinxLam);
  tLamPart = 2*T0*lambda3 - 2*T1*lambda3 + lambda3*dotProd1 -3*T0*lambda2 + 3*T1*lambda2 - lambda2*dotProd2 + lambda*dotProd4 + T0;
  return tLamPart + indexRef*norm1;
}


double projectedGradient_fromEdge(double lambda0, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double tol, double maxIter, double indexRef){
  // projected gradient descent for an update in which the segment x0x1 is on the boundary but xHat
  // if fully contained in a region with index indexRef
  double grad_cur, grad_prev, step, alpha, lam_prev, lam_cur, test;
  int i;
  i = 1;
  alpha = 0.001;
  lam_prev = lambda0;
  grad_cur = der_fromEdge(lambda0, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  grad_prev = der_fromEdge(lambda0, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
  if(fabs(grad_cur) > tol){
    test = lam_prev - alpha*grad_cur/fabs(grad_cur);
  }
  else{
    test = lam_prev;
  }
  if(test>1){
    lam_cur = 1;
  }
  if(test<0){
    lam_cur = 0;
  }
  else{
    lam_cur = test;
  }
  grad_cur = der_fromEdge(lam_cur, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);

  while(i<maxIter & fabs(grad_cur)>tol & fabs(lam_cur - lam_prev)>0) {
    printf("\n\nIteration %d\n", i);
    alpha = backTr_fromEdge(0.08, grad_cur, lam_cur, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
    printf("Step size %lf\n", alpha);
    test = lam_prev - alpha*grad_cur/fabs(grad_cur);
    printf("Direction %lf\n", -alpha*grad_cur/grad_cur);
    printf("Goes from %lf  to   %lf\n", lam_cur, test);
    if(test<0){
      test = 0;
    }
    if(test>1){
      test = 1;
    }
    grad_prev = grad_cur;
    lam_prev = lam_cur;
    lam_cur = test;
    grad_cur = der_fromEdge(lam_cur, T0, grad0, B0, T1, grad1, B1, x0, x1, xHat, indexRef);
    printf("Lambda previous %lf \n", lam_prev);
    printf("Lambda current %lf \n", lam_cur);
    printf("Current derivative: %lf \n", grad_cur);
    i ++;
  }
  
  return lam_cur;
}


