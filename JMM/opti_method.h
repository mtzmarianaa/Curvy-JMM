#pragma once

void linearInterpolation(double param, double from[2], double to[2], double interpolation[2]);

void der_linearInterpolation(double param, double from[2], double to[2], double der_interpolation[2]);

void hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double interpolation[2]);

void grad_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double gradient[2]);

void secondDer_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double secondDer[2]);

double arclength_hermiteSimpson(double a, double b, double from[2], double to[2], double grad_from[2], double grad_to[2]);

double hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]);

double der_hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]);

double der_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef);

double backTr_fromEdge(double alpha0, double d, double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef);

double fobjective_fromEdge(double lambda, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double indexRef);

double projectedGradient_fromEdge(double lambda0, double lambdaMin, double lambdaMax, double T0, double grad0[2], double B0[2], double T1, double grad1[2], double B1[2], double x0[2], double x1[2], double xHat[2], double tol, int maxIter, double indexRef);

double der_freeSpace(double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef);

double backTr_freeSpace(double alpha0, double d, double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef);

double fobjective_freeSpace(double lambda, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double indexRef);

double projectedGradient_freeSpace(double lambda0, double lambdaMin, double lambdaMax, double TA, double gradA[2], double TB, double gradB[2], double xA[2], double xB[2], double xHat[2], double tol, int maxIter, double indexRef);

void grad_twoStep(double gradient[2], double lambda, double mu, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02);

double backTr_TwoStep(double alpha0, double d[2], double lambda, double mu, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02);

double fobjective_TwoStep(double lambda, double mu, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02);

void projectedGradient_TwoStep(double optimizers[2], double lambdaMin, double lambdaMax, double muMin, double muMax, double T0, double grad0[2], double T1, double grad1[2], double x0[2], double x1[2], double x2[2], double xHat[2], double B0[2], double B2[2], double indexRef_01, double indexRef_02, double tol, int maxIter);

double t0_ofMu(double mu, double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]);

double der_t0_ofMu(double mu, double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]);

double findMumin_shootCr(double mu0, double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double tol, int maxIter);

double t_ofMu(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]);

double der_t_ofMu(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2]);

double backTr_find_minMu(double mu, double alpha0, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double tol, double maxIter);

double find_minMu(double mu0, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double tol, double maxIter);

double der_shootCr(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef);

double backTr_shootCr(double alpha0, double d, double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef);

double projectedGradient_shootCr(double mu0, double muMin, double muMax, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double tol, int maxIter, double indexRef);

double fobjective_shootCr(double mu, double xA[2], double xB[2], double xHat[2], double xR[2], double BHat[2], double BR[2], double TA, double TB, double gradA[2], double gradB[2], double indexRef);

