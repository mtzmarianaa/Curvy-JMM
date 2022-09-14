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
    vec2_substraction( xlambda, xHat, xLamMinxHat ); // xlambda - xHat
    // printf("            xlambda-xHat = (%fl, %fl)\n", xLamMinxHat[0], xLamMinxHat[1]);
    normAux = l2norm( xLamMinxHat ); //   || xlambda - xHat ||
    // printf("            norm(xlambda-xHat) = %fl\n", normAux);

    gApprox = min1Lambda*T0 + lambda*T1 + indexRef*normAux;

    return gApprox;
}

double gPrime(double T0, double T1, double lambda, double x0[2], double x1[2], double xHat[2], double indexRef){
    // auxiliary function to compute the function gPrime
    double x1Minx0[2], dotProduct, normAux, gPrim;
    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xLamMinxHat[2], min1Lambda, gApprox;

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
    vec2_substraction( xlambda, xHat, xLamMinxHat ); // xlambda - xHat
    normAux = l2norm( xLamMinxHat ); //   || xlambda - xHat ||

    x1Minx0[0] = 0;
    x1Minx0[1] = 0;

    dotProduct = 0;
    normAux= 0;

    vec2_substraction(x1, x0, x1Minx0);

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

    vec2_substraction( xmu, xlambda, xMuMinxLam ); // xmu - xlambda

    vec2_substraction( xHat, xmu, xHatMinxMu ); // xhat - xmu

    normAux1 = l2norm( xMuMinxLam ); //   || xmu - xlambda ||

    normAux2 = l2norm( xHatMinxMu ); //   || xhat - xmu ||

    gApprox = min1Lambda*T0 + lambda*T1 + indexRef_01*normAux1 + indexRef_02*normAux2;

    return gApprox;
}

void gradient_2Regions(double grad[2], double T0, double T1, double lambda, double mu, double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02) {
    // auxiliary function to calculate the gradient of eikApproxLin_2Regions, the first entry corresponds to the partial with respect to lambda
    // and the second entry to the partial with respect to mu
    double xlambda[2], firstPartxLambda[2], secondPartxLambda[2], xMuMinxLam[2], normAux1, min1Lambda, gApprox;
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

    vec2_substraction( xmu, xlambda, xMuMinxLam ); // xmu - xlambda

    vec2_substraction( xHat, xmu, xHatMinxMu ); // xhat - xmu

    normAux1 = l2norm( xMuMinxLam ); //   || xmu - xlambda ||

    normAux2 = l2norm( xHatMinxMu ); //   || xhat - xmu ||

    vec2_substraction(x1, x0, x1Minx0); // x1 - x0

    vec2_substraction(x2, x0, x2Minx0); // x2 - x0

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
    double lambda, mu, grad[2], step[2], yk[2], projectedGradient[2], Tcurrent, t;
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
        vec2_substraction( optimizers, step, yk );
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