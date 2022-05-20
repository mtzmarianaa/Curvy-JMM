#pragma once

#include "eik_grid.h"


void FMM_2D( eik_gridS *eikonal_g );

static double speed(double x, double y);

static void ActualSolution(double x_min, double y_min, int start[2], double *trueSolution, int M, int N, double h);

static double onePointUpdate(  eik_gridS *eikonal_g, int coordinate  );

static void generateDataForPlot(double *x, int M, int N, char data_name[10]);