#pragma once

double der_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef);

double backTr_fromEdge(double alpha0, double d, double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef);

double fobjective_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef);

double projectedGradient_fromEdge(double lambda0, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double tol, double maxIter, double indexRef);

double der_freeSpace(double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef);

double backTr_freeSpace(double alpha0, double d, double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef);

double fobjective_freeSpace(double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef);

double projectedGradient_freeSpace(double lambda0, double lambdaMin, double lambdaMax, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double tol, double maxIter, double indexRef);



  



double eikApprox_freeSpace(double T0, double T1, double grad0[2], double grad1[2],
			   double lambda, double x0[2], double x1[2], double xHat[2], double indexRef);

double secant_freeSpace(double lambda0, double lambda1, double T0, double T1, double grad0[2], double grad1[2],
			double x0[2], double x1[2], double xHat[2], double tol, int maxIter, double indexRef);

double eikApproxLin(double T1, double T0, double lambda, double x0[2], double x1[2],
		    double xHat[2], double indexRef);

double gPrime(double T1, double T0, double lambda, double x0[], double x1[],
	      double xHat[], double indexRef);

double secant_2D(double lambda0, double lambda1, double T0, double T1, double x0[],
		 double x1[], double xHat[], double tol, int maxIter, double indexRef);

double eikApproxLin_2Regions(double T0, double T1, double lambda, double mu,
			     double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02);

void gradient_2Regions(double grad[2], double T0, double T1, double lambda, double mu,
		       double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02);

void projectedGradientDescent(double optimizers[2], double T0, double T1, double x0[2],
			      double x1[2], double x2[2], double xHat[2], double tol,
			      int maxIter, double indexRef_01, double indexRef_02);

double backtracking(double optimizers[2], double direction[2], double T0, double T1,
		    double x0[2], double x1[2], double x2[2], double xHat[2], double indexRef_01, double indexRef_02);


double gPrimeCubic(double T0, double T1, double grad0[2], double grad1[2], double lambda,
		   double x0[2], double x1[2], double xHat[2], double indexRef);

double eikApproxCubic_2Regions(double T0, double T1, double grad0[2], double grad1[2],
			       double lambda, double mu, double x0[2], double x1[2],
			       double x2[2], double xHat[2], double indexRef_01, double indexRef_02);

void gradientCubic_2Regions(double grad[2], double T0, double T1, double grad0[2], double grad1[2],
			    double lambda, double mu, double x0[2], double x1[2], double x2[2],
			    double xHat[2], double indexRef_01, double indexRef_02);

void projectedGradientDescentCubic(double optimizers[2], double T0, double T1, double grad0[2], double grad1[2],
			      double x0[2], double x1[2], double x2[2], double xHat[2], double tol,
			      int maxIter, double indexRef_01, double indexRef_02);




