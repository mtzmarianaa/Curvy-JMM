/*
USEFUL LINEAR ALGEBRA COMPUTATIONS
*/
#pragma once
#include <math.h>

double l2norm(double x0, double x1)
{
    return sqrt( pow(x0, 2) + pow(x1, 2)  );
}

double dotProd(double x0, double x1, double y0, double y1)
{
    return x0*y0 + x1*y1;
}

