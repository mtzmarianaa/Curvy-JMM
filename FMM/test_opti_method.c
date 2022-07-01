#include "opti_method.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){


    double lambda0, lambda1, x0[2], x1[2], xHat[2], tol, lamb_opt, T0, T1;
    int maxIter;
    lambda0 = 0.1;
    lambda1 = 0.8;

    x0[0] = 1;
    x0[1] = 0;
    T0 = 1;

    x1[0] = 0;
    x1[1] = 1;
    T1 = 1;

    xHat[0] = 1;
    xHat[1] = 1;

    tol = 0.001;
    maxIter = 25 ;
    lamb_opt = secant_2D(lambda0, lambda1, T0, T1, x0, x1, xHat, tol, maxIter);
    printf("The optimum lambda should be 0.5 (middle point of the base of the triangle)\n");
    printf("x0: %f, %f \n", x0[0], x0[1]);
    printf("T(x0): %f \n", T0);
    printf("x1: %f, %f \n", x1[0], x1[1]);
    printf("T(x0): %f \n", T0);
    printf("Lambda 1: %f \n", lambda1);
    printf("Lambda 2: %f \n", lambda0);
    printf("Optimum lambda: %f \n", lamb_opt);

    printf("\n\n\n");

    lambda0 = 0.0;
    lambda1 = 1.0;

    x0[0] = 1;
    x0[1] = 0;
    T0 = 1;

    x1[0] = 2;
    x1[1] = 1;
    T1 = sqrt(5);

    xHat[0] = 2;
    xHat[1] = 0;
    tol = 0.001;
    maxIter = 25 ;
    lamb_opt = secant_2D(lambda0, lambda1, T0, T1, x0, x1, xHat, tol, maxIter);
    printf("The optimum lambda should be 0\n");
    printf("x0: %f, %f \n", x0[0], x0[1]);
    printf("T(x0): %f \n", T0);
    printf("x1: %f, %f \n", x1[0], x1[1]);
    printf("T(x1): %f \n", T1);
    printf("Lambda 1: %f \n", lambda1);
    printf("Lambda 2: %f \n", lambda0);
    printf("Optimum lambda: %f \n", lamb_opt);





}


