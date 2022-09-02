/*
USEFUL LINEAR ALGEBRA COMPUTATIONS
*/
#pragma once

#define MATRIX_2x2(a00, a01, a10, a11) (matrix_2x2) {{a00, a01}, {a10, a11}} // we define a matrix

typedef double matrix_2x1[2]; // this a 2x1 vector
typedef matrix_2x1 matrix_2x2[2]; // these are 2 2x1 vectors to form a 2x2 matrix

double l2norm(double x[]);

double lInfnorm(double x[]);

void projection01Cube(double x[], double projectedx[]);

double dotProd(double x[], double y[]);

void scalar_times_2vec(double alpha, double x[], double output[]);

void vec2_addition(double x[], double y[], double output[]);

void vec2_substraction(double x[], double y[], double output[]);

double angleThreePoints(double A[], double B[], double C[]);

double toc();

double determinant(matrix_2x2 const A);

void inverse2x2(matrix_2x2 const A, matrix_2x2 Ainv);