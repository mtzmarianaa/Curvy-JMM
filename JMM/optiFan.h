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
  double (*points_fan)[2]; // coordinates of the points on the fan
  double (*boundary_x0)[2]; // divisions inside the triangle fan, length n+2
  double *indicesRef; // n+2 different indices of refraction
} optiFanS;

void optiFan_alloc(optiFanS **optiFan);

void optiFan_dealloc(optiFanS **optiFan);

void optiFan_init(optiFanS *optiFan, int nRegions, double x0[2], double T0; double x1[2], double T1, double xHat[2], double *(points_fan)[2], double (*boundary_x0)[2], double *indicesRef);

//////////// AUXILIARY FUNCTIONS

void hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double interpolation[2]);

void grad_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double gradient[2]);

void secondDer_hermite_interpolationSpatial(double param, double from[2], double to[2], double grad_from[2], double grad_to[2], double secondDer[2]);

double arclength_hermiteSimpson(double a, double b, double from[2], double to[2], double grad_from[2], double grad_to[2]);

double hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]);

double der_hermite_interpolationT(double param, double xA[2], double xB[2], double TA, double TB, double gradA[2], double gradB[2]);


