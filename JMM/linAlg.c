/*
USEFUL LINEAR ALGEBRA COMPUTATIONS
*/
#include "linAlg.h"

#include <math.h>
#include <time.h>

double l2norm(double x[])
{
    double result;
    result = sqrt( pow(x[0], 2) + pow(x[1], 2)  );
    return result;
}

double lInfnorm(double x[])
{
    double result;
    result = fabs(x[0]);
    if(fabs(x[0]) < fabs(x[1]) ){
        result = fabs(x[1]);
    }
    return result;
}

void projection01Cube(double x[], double projectedx[]){
    // projection of x in R2 onto the square 0 1 in R2
    projectedx[0] = x[0];
    projectedx[1] = x[1];
    if( x[0] < 0 ){
        projectedx[0] = 0;
    }
    if( x[1] < 0) {
        projectedx[1] = 0;
    }
    if( x[0] > 1  ){
        projectedx[0] = 1;
    }
    if( x[1] > 1){
        projectedx[1] = 1;
    }
}

double dotProd(double x[], double y[])
{
    return x[0]*y[0] + x[1]*y[1];
}

void scalar_times_2vec(double alpha, double x[], double output[])
{
    output[0] = alpha*x[0];
    output[1] = alpha*x[1];
}

void vec2_addition(double x[], double y[], double output[])
{
    output[0] = x[0] + y[0];
    output[1] = x[1] + y[1];
}

void vec2_subtraction(double x[], double y[], double output[])
{
    output[0] = x[0] - y[0];
    output[1] = x[1] - y[1];
}

double angleThreePoints(double A[], double B[], double C[]) {
    double BA[2], BC[2];
    vec2_subtraction( A, B, BA );
    vec2_subtraction( C, B, BC );
    double normBA, normBC;
    normBA = l2norm(BA);
    normBC = l2norm(BC);
    return acos(dotProd(BA, BC)/(normBA*normBC));
}


double toc() {
  static clock_t t1 = 0;
  clock_t t0 = t1;
  t1 = clock();
  return ((double)t1 - (double)t0)/CLOCKS_PER_SEC;
}

double determinant(matrix_2x2 const A){
    return A[0][0]*A[1][1] - A[0][1]*A[1][0];
}

void inverse2x2(matrix_2x2 const A, matrix_2x2 Ainv){
    double determinantA;
    determinantA = determinant(A);
    Ainv[0][0] = A[1][1]/determinantA;
    Ainv[0][1] = -A[0][1]/determinantA;
    Ainv[1][0] = -A[1][0]/determinantA;
    Ainv[1][1] = A[0][0]/determinantA;
}

void matrixXvec2x2(matrix_2x2 const A, double x[2], double sol[2]){
  sol[0] = A[0][0]*x[0] + A[0][1]*x[1];
  sol[1] = A[1][0]*x[0] + A[1][1]*x[1];
}

double min(double a, double b){
  if(a <= b){
    return a;
  }
  else{
    return b;
  }
}

