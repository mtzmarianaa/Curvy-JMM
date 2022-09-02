#include "opti_method.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(){
    void projectedGradientDescent(double optimizers[2], double T0, double T1, double x0[2], double x1[2], double x2[2], double xHat[2], double tol, int maxIter, double indexRef_01, double indexRef_02);

    double optimizers[2], T0, T1, x0[2], x1[2], x2[2], xHat[2], tol, indexRef_01, indexRef_02;
    int maxIter;

    maxIter = 50;
    tol = 0.0000001;
    indexRef_01 = 1.0;
    indexRef_02 = 1.453;

    x0[0] = -2.0;
    x0[1] = 1.0;
    x1[0] = 0.0;
    x1[1] = 1.0;
    x2[0] = 0.0;
    x2[1] =2.0;
    xHat[0] = -1.0;
    xHat[1] = 3.0;

    T0 = sqrt(2);
    T1 = sqrt(2);

    projectedGradientDescent(optimizers, T0, T1, x0, x1, x2, xHat, tol, maxIter, indexRef_01, indexRef_02);





}


