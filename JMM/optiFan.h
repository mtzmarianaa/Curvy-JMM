#pragma once

typedef struct {
  // all the information needed
  // this info was gather after marching
  // around the triangle fan
  int nRegions; // number of indices of refraction inside the triangle fan
  double x0[2]; // recently accepted point
  double T0; // eikonal value at x0
  double x1[2]; // farthest neighbor of x0 set to valid
  double T1; // eikonal value at x1
  double xHat[2]; // point we want to update
  double (*points_fan)[2]; // coordinates of the points on the fan (total = n)
  double (*B_x0)[2]; // gradients at x0 from the boundary (total = n+1)
  double (*B_xk)[2]; // gradients at xk from the boundary (total = n+1)
  double *indicesRef; // n+2 different indices of refraction
  int *types; // types of curved boundaries inside the triangle fan (type 1, 2, 3, or 4)
} optiFanS;

void optiFan_alloc(optiFanS **optiFan);

void optiFan_dealloc(optiFanS **optiFan);

void optiFan_init(optiFanS *optiFan, int nRegions, double x0[2], double T0; double x1[2], double T1, double xHat[2], double *(points_fan)[2], double (*B_x0)[2], double (*B_xk)[2], double *indicesRef);

//////////// AUXILIARY FUNCTIONS

void hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double interpolation[2]);

void grad_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double gradient[2]);

void secondDer_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double secondDer[2]);

double arclength_hermiteSimpson(double a, double b, double from[2], double to[2], double grad_from[2], double grad_to[2]);

double hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]);

double der_hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]);

// functions to find lambdaMin and lambdaMax for types 1,2, and 4

double t1_ofLam(double lambda, double x0[2], double B0[2], double ykPrime[2], double Bk_mu[2], double x_k1[2], double B_k1[2]);

double t1Prime_ofLam(double lambda, double x0[2], double B0[2], double Bk_mu[2], double x_k1[2], double B_k1[2]);

double backTr_t1(double alpha0, double d, double lambda, double x0[2], double B0[2], double ykPrime[2], double Bk_mu[2], double x_k1[2], double B_k1[2]);

double lambda_fromt1(double lambda0, double x0[2], double B0[2], double ykPrime[2], double Bk_mu[2], double x_k1[2], double B_k1[2], double tol, int maxIter);

double t2_ofLam(double lambda, double x0, double B0[2], double ykPrime[2], double x_k1[2], double B_k1[2]);

double t2Prime_ofLam(double lambda, double x0, double B0[2], double x_k1[2], double B_k1[2]);

double lambda_fromt2(double lambda0, double x0, double B0[2], double ykPrime[2], double x_k1[2], double B_k1[2], double tol, int maxIter);


