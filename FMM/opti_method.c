/*
Optimization methods for the 2D FMM
*/

#include "opti_method.h"
#include "linAlg.h"
#include "SoSFunction.h"

#include <math.h>
#include <stdio.h>

void secant_2D(double lambda0, double lambda1, double x0[], double x1[], double xHat[], double tol, int maxIter){
    // This method is the implementation of the secant method for the 2D fmm using the
    // function defined in SoSFunction.h as the speed of sound
    int k;
    k = 1;
    double g_prime, temp1[2], temp2[2], temp3[2], temp4[2], temp5[2], temp6[2], temp7[2], norm1, norm2, lam;
    temp1[0] = 0; // xHat - x1
    temp1[1] = 0;
    temp2[0] = 0; // x1 - x0
    temp2[1] = 0;
    temp3[0] = 0; // lambda1*(x1-x0)
    temp3[1] = 0;
    temp4[0] = 0; // lambda0*(x1-x0)
    temp4[1] = 0;
    temp5[0] = 0; // xHat - x1 + lambda1*(x1-x0) = temp1 + temp3;
    temp5[1] = 0;
    temp6[0] = 0; // xHat - x1 + lambda0*(x1 - x0) = temp1 + temp4;
    temp6[1] = 0;
    temp7[0] = 0;
    temp7[1] = 0;
    vec2_substraction(xHat, x1, temp1);
    vec2_substraction(x1, x0, temp2);
    scalar_times_2vec(lambda1, temp2, temp3);
    scalar_times_2vec(lambda0, temp2, temp4);
    vec2_addition(temp1, temp3, temp5);
    vec2_addition(temp1, temp4, temp6);
    norm1 = l2norm(temp5);
    norm2 = l2norm(temp6);
    g_prime = 1*dotProd(temp2, temp5); // numerator of g_prime is enough to know if g_prime is zero or not
    while(k < maxIter & g_prime!=0){
        vec2_addition(temp1, temp3, temp7);
        lam = lambda1 - s_function(xHat)*dotProd(temp7, temp2)*(lambda1-lambda0)/(norm2*dotProd(temp2, temp5) + 1*norm1*dotProd(temp2, temp6) );
        vec2_substraction(xHat, x1, temp1);
        vec2_substraction(x1, x0, temp2);
        scalar_times_2vec(lambda1, temp2, temp3);
        scalar_times_2vec(lambda0, temp2, temp4);
        vec2_addition(temp1, temp3, temp5);
        vec2_addition(temp1, temp4, temp6);
        norm1 = l2norm(temp5);
        norm2 = l2norm(temp6);
        g_prime = 1*dotProd(temp2, temp5);
        printf("Lamba 1: %f \n", lambda1);
        printf("Lamba 2: %f \n", lambda0);
        printf("Lam : %f \n", lam);
        lambda0 = lambda1;
        lambda1 = lam;
        k ++;
    }
}

