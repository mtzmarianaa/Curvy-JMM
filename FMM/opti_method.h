#pragma once

double eikApproxLin(double T1, double T0, double lambda, double x0[2], double x1[2], double xHat[2], double indexRef);

double gPrime(double T1, double T0, double lambda, double x0[], double x1[], double xHat[], double indexRef);

double secant_2D(double lambda0, double lambda1, double T0, double T1, double x0[], double x1[], double xHat[], double tol, int maxIter, double indexRef);

double eikApproxLin_2Regions(double T0, double T1, double lambda, double mu, double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02);

void gradient_2Regions(double grad[2], double T0, double T1, double lambda, double mu, double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02);

void projectedGradientDescent(double optimizers[2], double T0, double T1, double x0[2], double x1[2], double x2[2], double xHat[2], double tol, int maxIter, double indexRef_01, double indexRef_02);

double backtracking(double optimizers[2], double direction[2], double T0, double T1, double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02);