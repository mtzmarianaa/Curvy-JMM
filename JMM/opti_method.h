#pragma once

double eikApproxCubic(double T0, double T1, double grad0[2], double grad1[2],
		      double lambda, double x0[2], double x1[2], double xHat[2], double indexRef);

double gPrimeCubic(double T0, double T1, double grad0[2], double grad1[2], double lambda,
		   double x0[2], double x1[2], double xHat[2], double indexRef);

double secantCubic_2D(double lambda0, double lambda1, double T0, double T1,
		      double grad0[2], double grad1[2],
		     double x0[2], double x1[2], double xHat[2],
		      double tol, int maxIter, double indexRef);
