/*
USEFUL LINEAR ALGEBRA COMPUTATIONS
*/
#pragma once

double l2norm(double x[]);

double dotProd(double x[], double y[]);

void scalar_times_2vec(double alpha, double x[], double output[]);

void vec2_addition(double x[], double y[], double output[]);

void vec2_substraction(double x[], double y[], double output[]);

double toc();
