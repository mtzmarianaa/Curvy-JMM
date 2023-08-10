#pragma once

#include <stdio.h>
#include <stdlib.h>

double eikApproxCubic(double T0, double T1, double grad0[2], double grad1[2],
		      double lambda, double x0[2], double x1[2], double xHat[2], double indexRef);

double gPrimeCubic(double T0, double T1, double grad0[2], double grad1[2], double lambda,
		   double x0[2], double x1[2], double xHat[2], double indexRef);

double stepSize(double d, double T0, double T1, double grad0[2], double grad1[2],
		double lambda, double x0[2], double x1[2], double xHat[2], double indexRef);

double gradDescentCubic_2D(double T0, double T1, double grad0[2], double grad1[2],
			   double x0[2], double x1[2], double xHat[2],
			   double tol, size_t maxIter, double indexRef);

double secantCubic_2D(double T0, double T1, double grad0[2], double grad1[2],
		      double x0[2], double x1[2], double xHat[2], double tol,
		      int maxIter, double indexRef);
