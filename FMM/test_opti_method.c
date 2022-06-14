#include "opti_method.h"

#include <stdio.h>

int main(){
    double lambda0, lambda1, x0[2], x1[2], xHat[2], tol;
    int maxIter;
    lambda0 = 0.0;
    lambda1 = 1.0;
    x0[0] = 1;
    x0[1] = 0;
    x1[0] = 0;
    x1[1] = 1;
    xHat[0] = 1;
    xHat[1] = 1;
    tol = 0.001;
    maxIter = 5;
    secant_2D(lambda0, lambda1, x0, x1, xHat, tol, maxIter);
    printf("\f", lambda1);
    printf("\f", lambda0);
}