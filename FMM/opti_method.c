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
    double x0Minx1[2], xHatMinx0[2], lamx0Minx1[2], vecAux[2], normAux, gApprox;
    printf("\nx0 is (%fl, %fl), x1 is (%fl, %fl), and lambda is %fl\n", x0[0], x0[1], x1[0], x1[1], lambda);

    x0Minx1[0] = 0;
    x0Minx1[1] = 0;

    xHatMinx0[0] = 0;
    xHatMinx0[1] = 0;

    lamx0Minx1[0] = 0;
    lamx0Minx1[1] = 0;

    vecAux[0] = 0;
    vecAux[1] = 0;

    normAux= 0;

    vec2_substraction(x0, x1, x0Minx1);
    vec2_substraction(xHat, x0, xHatMinx0);
    scalar_times_2vec(lambda, x0Minx1, lamx0Minx1);
    vec2_addition(xHatMinx0, lamx0Minx1, vecAux);
    normAux = l2norm(vecAux);
    printf("The norm of this path is %fl\n", normAux);
    gApprox = lambda*(T1 - T0) + T0 + s_function_threeSections(xHat, regionIndex)*normAux;

    return gApprox;
}

double gPrime(double T1, double T0, double lambda, double x0[2], double x1[2], double xHat[2], int regionIndex){
    // auxiliary function to compute the function gPrime
    double x0Minx1[2], xHatMinx0[2], lamx0Minx1[2], vecAux[2], dotProduct, normAux, gPrim;

    x0Minx1[0] = 0;
    x0Minx1[1] = 0;

    xHatMinx0[0] = 0;
    xHatMinx0[1] = 0;

    lamx0Minx1[0] = 0;
    lamx0Minx1[1] = 0;

    vecAux[0] = 0;
    vecAux[1] = 0;

    dotProduct = 0;
    normAux= 0;

    vec2_substraction(x0, x1, x0Minx1);
    vec2_substraction(xHat, x0, xHatMinx0);
    scalar_times_2vec(lambda, x0Minx1, lamx0Minx1);
    vec2_addition(xHatMinx0, lamx0Minx1, vecAux);
    dotProduct = dotProd(x0Minx1, vecAux);
    normAux = l2norm(vecAux);

    gPrim = T1 - T0 + s_function_threeSections(xHat, regionIndex)*dotProduct/normAux;

    return gPrim;

}

double secant_2D(double lambda0, double lambda1, double T0, double T1, double x0[2], double x1[2], double xHat[2], double tol, int maxIter, int regionIndex){
    // This method is the implementation of the secant method for the 2D fmm using the
    // function defined in SoSFunction.h as the speed of sound
    int k = 1;
    double lam, gPrime0, gPrime1;

    gPrime1 = gPrime(T1, T0, lambda1, x0, x1, xHat, regionIndex);
    
    while(k < maxIter & gPrime1>tol){
        // u
        gPrime0 = gPrime(T1, T0, lambda0, x0, x1, xHat, regionIndex);
        gPrime1 = gPrime(T1, T0, lambda1, x0, x1, xHat, regionIndex);
        lam = lambda1 - gPrime1*(lambda1 - lambda0)/( gPrime1 - gPrime0 );
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
    //printf("Optimum lambda found %fl\n", lambda1);
    return lambda1;
}

