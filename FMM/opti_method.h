#pragma once

double eikApproxLin(double T1, double T0, double lambda, double x0[2], double x1[2], double xHat[2], int regionIndex);

double gPrime(double T1, double T0, double lambda, double x0[], double x1[], double xHat[], int regionIndex);

double secant_2D(double lambda0, double lambda1, double T0, double T1, double x0[], double x1[], double xHat[], double tol, int maxIter, int regionIndex);