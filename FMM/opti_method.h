#pragma once

double gPrime(double T1, double T0, double lambda, double x0[], double x1[], double xHat[]);

double secant_2D(double lambda0, double lambda1, double T0, double T1, double x0[], double x1[], double xHat[], double tol, int maxIter);