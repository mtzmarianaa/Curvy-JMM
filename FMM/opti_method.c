/*
Optimization methods for the 2D FMM
*/

#include "opti_method.h"
#include "linAlg.h"
#include "SoSFunction.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>


double eikApproxLin(double T1, double T0, double lambda, double x0[2], double x1[2], double xHat[2], int regionIndex) {
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
    printf("   aux printing eikApproxLin\n");
    printf("            1-lambda = %fl\n", min1Lambda);
    scalar_times_2vec( min1Lambda, x0, firstPartxLambda ); // (1-lambda)*x0
    printf("            (1-lambda)*x0 = (%fl,%fl)\n", firstPartxLambda[0], firstPartxLambda[1]);
    scalar_times_2vec(lambda, x1, secondPartxLambda); // lambda*x1
    printf("            lambda*x1 = (%fl,%fl)\n", secondPartxLambda[0], secondPartxLambda[1]);
    vec2_addition( firstPartxLambda, secondPartxLambda, xlambda ); // (1-lambda)*x0 + lambda*x1
    printf("            xlambda = (%fl,%fl)\n", xlambda[0], xlambda[1]);
    vec2_substraction( xlambda, xHat, xLamMinxHat ); // xlambda - xHat
    printf("            xlambda-xHat = (%fl, %fl)\n", xLamMinxHat[0], xLamMinxHat[1]);
    normAux = l2norm( xLamMinxHat ); //   || xlambda - xHat ||
    printf("            norm(xlambda-xHat) = %fl\n", normAux);

    gApprox = min1Lambda*T0 + lambda*T1 + s_function_threeSections(xHat, regionIndex)*normAux;

    return gApprox;
}

double gPrime(double T1, double T0, double lambda, double x0[2], double x1[2], double xHat[2], int regionIndex){
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

    gPrim = T1 - T0 + s_function_threeSections(xHat, regionIndex)*dotProduct/normAux;

    return gPrim;

}

double secant_2D(double lambda0, double lambda1, double T0, double T1, double x0[2], double x1[2], double xHat[2], double tol, int maxIter, int regionIndex){
    // This method is the implementation of the secant method for the 2D fmm using the
    // function defined in SoSFunction.h as the speed of sound
    int k = 1;
    double lam, gPrime0, gPrime1;

    gPrime1 = gPrime(T1, T0, lambda1, x0, x1, xHat, regionIndex);
    printf("Initial, with lambda = %fl   the value of gPrime is: %fl\n", lambda0, gPrime(T1, T0, lambda0, x0, x1, xHat, regionIndex));
    printf("Initial, with lambda = %fl   the value of gPrime is: %fl\n", lambda1, gPrime1);
    
    while(k < maxIter & fabs(gPrime1)>tol){
        // u
        gPrime0 = gPrime(T1, T0, lambda0, x0, x1, xHat, regionIndex);
        gPrime1 = gPrime(T1, T0, lambda1, x0, x1, xHat, regionIndex);
        lam = lambda1 - gPrime1*(lambda1 - lambda0)/( gPrime1 - gPrime0 );
        printf("In this iteration lambda_i has value of %fl", lam);
        printf("    |  value of T(xhat): %fl   with %fl as index of regraction\n", eikApproxLin(T1, T0, lam, x0, x1, xHat, regionIndex), s_function_threeSections(xHat, regionIndex) );
        printf("    |   In this iteration the value of gPrime is: %fl \n", gPrime(T1, T0, lam, x0, x1, xHat, regionIndex));
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
    printf("Optimum lambda found %fl\n", lambda1);
    printf("With gPrime value of: %fl\n", gPrime(T1, T0, lambda1, x0, x1, xHat, regionIndex));
    return lambda1;
}

